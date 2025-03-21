version 1.0

# Merge and process RDS object with QC, filtering, and normalization, and convert to adata object, integrate and cluster

import "../../../wf-common/wdl/tasks/write_cohort_sample_list.wdl" as WriteCohortSampleList
import "../../integrate_data/integrate_data.wdl" as IntegrateData
import "../../../wf-common/wdl/tasks/upload_final_outputs.wdl" as UploadFinalOutputs

workflow cohort_analysis {
	input {
		String cohort_id
		Array[Array[String]] project_sample_ids
		Array[File] preprocessed_rds_objects

		# If provided, these files will be uploaded to the staging bucket alongside other intermediate files made by this workflow
		Array[String] preprocessing_output_file_paths = []

		# Filtering parameters
		File cell_type_markers_list
		Float min_genes_detected_in_percent_segment

		# Integrate and cluster parameters
		Int n_top_genes
		Int n_comps
		String batch_key
		Float leiden_resolution

		String workflow_name
		String workflow_version
		String workflow_release
		String run_timestamp
		String raw_data_path_prefix
		Array[String] staging_data_buckets
		String billing_project
		String container_registry
		String zones
	}

	String sub_workflow_name = "cohort_analysis"
	String sub_workflow_version = "1.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version, workflow_release]]

	String raw_data_path = "~{raw_data_path_prefix}/~{sub_workflow_name}/~{sub_workflow_version}/~{run_timestamp}"

	call WriteCohortSampleList.write_cohort_sample_list {
		input:
			cohort_id = cohort_id,
			project_sample_ids = project_sample_ids,
			billing_project = billing_project,
			workflow_info = workflow_info,
			raw_data_path = raw_data_path,
			container_registry = container_registry,
			zones = zones
	}

	scatter (preprocessed_rds_object in preprocessed_rds_objects) {
		call process {
			input:
				preprocessed_rds_object = preprocessed_rds_object,
				cell_type_markers_list = cell_type_markers_list,
				min_genes_detected_in_percent_segment = min_genes_detected_in_percent_segment,
				raw_data_path = raw_data_path,
				workflow_info = workflow_info,
				billing_project = billing_project,
				container_registry = container_registry,
				zones = zones
		}

		call rds_to_adata {
			input:
				processed_rds_object = process.processed_rds_object,
				container_registry = container_registry,
				zones = zones
		}
	}

	call merge_and_prep {
		input:
			cohort_id = cohort_id,
			processed_adata_objects = rds_to_adata.processed_adata_object,
			n_top_genes = n_top_genes,
			n_comps = n_comps,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call IntegrateData.integrate_data {
		input:
			cohort_id = cohort_id,
			processed_adata_object = merge_and_prep.merged_adata_object, #!FileCoercion
			n_comps = n_comps,
			batch_key = batch_key,
			leiden_resolution = leiden_resolution,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call UploadFinalOutputs.upload_final_outputs as upload_preprocess_files {
		input:
			output_file_paths = preprocessing_output_file_paths,
			staging_data_buckets = staging_data_buckets,
			staging_data_path = "~{workflow_name}/preprocess",
			billing_project = billing_project,
			zones = zones
	}

	Array[String] cohort_analysis_final_output_paths = flatten([
		[
			write_cohort_sample_list.cohort_sample_list
		],
		process.segment_gene_detection_plot_png,
		process.gene_detection_rate_csv,
		process.q3_negprobe_plot_png,
		process.normalization_plot_png,
		rds_to_adata.processed_adata_object,
		[
			merge_and_prep.merged_adata_object
		],
		[
			integrate_data.umap_cluster_plots_png
		]
	]) #!StringCoercion

	call UploadFinalOutputs.upload_final_outputs as upload_cohort_analysis_files {
		input:
			output_file_paths = cohort_analysis_final_output_paths,
			staging_data_buckets = staging_data_buckets,
			staging_data_path = "~{workflow_name}/~{sub_workflow_name}",
			billing_project = billing_project,
			zones = zones
	}

	output {
		File cohort_sample_list = write_cohort_sample_list.cohort_sample_list #!FileCoercion

		# Processed RDS object and plots
		Array[File] processed_rds_object = process.processed_rds_object
		Array[File] segment_gene_detection_plot_png = process.segment_gene_detection_plot_png #!FileCoercion
		Array[File] gene_detection_rate_csv = process.gene_detection_rate_csv #!FileCoercion
		Array[File] q3_negprobe_plot_png = process.q3_negprobe_plot_png #!FileCoercion
		Array[File] normalization_plot_png = process.normalization_plot_png #!FileCoercion

		# Converted AnnData object
		Array[File] processed_adata_object = rds_to_adata.processed_adata_object

		# Merged and prepped AnnData object
		File merged_adata_object = merge_and_prep.merged_adata_object #!FileCoercion

		# Integrate data outputs
		File integrated_adata_object = integrate_data.integrated_adata_object
		File clustered_adata_object = integrate_data.clustered_adata_object
		File umap_cluster_plots_png = integrate_data.umap_cluster_plots_png

		Array[File] preprocess_manifest_tsvs = upload_preprocess_files.manifests #!FileCoercion
		Array[File] cohort_analysis_manifest_tsvs = upload_cohort_analysis_files.manifests #!FileCoercion
	}
}

task process {
	input {
		File preprocessed_rds_object

		File cell_type_markers_list

		Float min_genes_detected_in_percent_segment

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int mem_gb = ceil(size([preprocessed_rds_object, cell_type_markers_list], "GB") * 2 + 20)
	Int disk_size = ceil(size([preprocessed_rds_object, cell_type_markers_list], "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		sample_id=$(basename ~{preprocessed_rds_object} | cut -d '.' -f 1)

		Rscript /opt/scripts/process.R \
			--sample-id "$sample_id" \
			--input ~{preprocessed_rds_object} \
			--celltype-markers ~{cell_type_markers_list} \
			--min-segment ~{min_genes_detected_in_percent_segment} \
			--output "$sample_id".processed.rds

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "$sample_id.segment_gene_detection_plot.png" \
			-o "$sample_id.gene_detection_rate.csv" \
			-o "$sample_id.q3_negprobe_plot.png" \
			-o "$sample_id.normalization_plot.png"
	>>>

	output {
		File processed_rds_object = "$sample_id.processed.rds"
		String segment_gene_detection_plot_png = "~{raw_data_path}/$sample_id.segment_gene_detection_plot.png"
		String gene_detection_rate_csv = "~{raw_data_path}/$sample_id.gene_detection_rate.csv"
		String q3_negprobe_plot_png = "~{raw_data_path}/$sample_id.q3_negprobe_plot.png"
		String normalization_plot_png = "~{raw_data_path}/$sample_id.normalization_plot.png"
	}

	runtime {
		docker: "~{container_registry}/spatial_r:1.0.0"
		cpu: 2
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 30
		zones: zones
	}
}

task rds_to_adata {
	input {
		File processed_rds_object

		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(processed_rds_object, "GB") * 2 + 20)
	Int disk_size = ceil(size(processed_rds_object, "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		sample_id=$(basename ~{processed_rds_object} | cut -d '.' -f 1)

		Rscript /opt/scripts/rds_to_adata.R \
			--input ~{processed_rds_object} \
			--output-prefix "$sample_id".processed
	>>>

	output {
		File processed_adata_object = "$sample_id.processed.h5ad"
	}

	runtime {
		docker: "~{container_registry}/spatial_r:1.0.0"
		cpu: 2
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 30
		zones: zones
	}
}

task merge_and_prep {
	input {
		String cohort_id
		Array[File] processed_adata_objects

		Int n_top_genes
		Int n_comps

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(processed_adata_objects, "GB") * 2 + 20)
	Int disk_size = ceil(size(processed_adata_objects, "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		python3 /opt/scripts/merge_and_prep_geomx.py \
			--adata-paths-input ~{sep=' ' processed_adata_objects} \
			--n-top-genes ~{n_top_genes} \
			--n-comps ~{n_comps} \
			--plots-prefix ~{cohort_id} \
			--adata-output ~{cohort_id}.merged.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.merged.h5ad"
	>>>

	output {
		String merged_adata_object = "~{raw_data_path}/~{cohort_id}.merged.h5ad"
	}

	runtime {
		docker: "~{container_registry}/spatial_py:1.0.0"
		cpu: 2
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 30
		zones: zones
	}
}
