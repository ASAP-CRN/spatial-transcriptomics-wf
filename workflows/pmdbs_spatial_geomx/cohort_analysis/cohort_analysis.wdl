version 1.0

# Merge and process RDS object with QC, filtering, and normalization, and convert to adata object and cluster

import "../../../wf-common/wdl/tasks/write_cohort_sample_list.wdl" as WriteCohortSampleList
import "../../../wf-common/wdl/tasks/upload_final_outputs.wdl" as UploadFinalOutputs

workflow cohort_analysis {
	input {
		String cohort_id
		Array[Array[String]] project_sample_ids
		Array[File] preprocessed_rds_objects

		# If provided, these files will be uploaded to the staging bucket alongside other intermediate files made by this workflow
		Array[String] preprocessing_output_file_paths = []

		File cell_type_markers_list

		# Filtering parameters
		Int min_genes_detected_in_percent_segment

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

	call merge {
		input:
			cohort_id = cohort_id,
			preprocessed_rds_objects = preprocessed_rds_objects,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call process {
		input:
			cohort_id = cohort_id,
			merged_rds_object = merge.merged_rds_object, #!FileCoercion
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
			cohort_id = cohort_id,
			processed_rds_object = process.processed_rds_object,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call IntegrateData.integrate_data {
		input:
			cohort_id = cohort_id,
			processed_adata_object = process.processed_adata_object, #!FileCoercion
			n_comps = n_comps,
			batch_key = batch_key,
			leiden_resolution = leiden_resolution,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call plot_spatial {
		input:
			cohort_id = cohort_id,
			clustered_adata_object = integrate_data.clustered_adata_object, #!FileCoercion
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call SpatialStatistics.spatial_statistics {
		input:
			cohort_id = cohort_id,
			clustered_adata_object = cluster.umap_cluster_adata_object, #!FileCoercion
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
		[
			merge_and_plot_qc_metrics.merged_adata_object,
			merge_and_plot_qc_metrics.qc_plots_png
		],
		[
			cluster.umap_cluster_plot_png,
		],
		[
			spatial_statistics.moran_top_10_variable_genes_csv,
			spatial_statistics.nhood_enrichment_plot_png,
			spatial_statistics.final_adata_object,
			spatial_statistics.co_occurrence_plot_png
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

		# Merged RDS objects
		File merged_rds_object = merge.merged_rds_object #!FileCoercion

		# Processed RDS object and plots
		File processed_rds_object = process.processed_rds_object
		File segment_gene_detection_plot_png = process.segment_gene_detection_plot_png #!FileCoercion
		File gene_detection_rate_csv = process.gene_detection_rate_csv #!FileCoercion
		File q3_negprobe_plot_png = process.q3_negprobe_plot_png #!FileCoercion
		File normalization_plot_png = process.normalization_plot_png #!FileCoercion

		# Converted AnnData object
		File processed_adata_object = rds_to_adata.processed_adata_object #!FileCoercion

		# Leiden clustered adata object and UMAP and spatial coordinates plots
		File umap_cluster_adata_object = cluster.umap_cluster_adata_object #!FileCoercion
		File umap_cluster_plot_png = cluster.umap_cluster_plot_png #!FileCoercion

		# Spatial statistics outputs
		File moran_adata_object = spatial_statistics.moran_adata_object
		File moran_top_10_variable_genes_csv = spatial_statistics.moran_top_10_variable_genes_csv
		File nhood_enrichment_adata_object = spatial_statistics.nhood_enrichment_adata_object
		File nhood_enrichment_plot_png = spatial_statistics.nhood_enrichment_plot_png
		File final_adata_object = spatial_statistics.final_adata_object
		File co_occurrence_plot_png = spatial_statistics.co_occurrence_plot_png

		Array[File] preprocess_manifest_tsvs = upload_preprocess_files.manifests #!FileCoercion
		Array[File] cohort_analysis_manifest_tsvs = upload_cohort_analysis_files.manifests #!FileCoercion
	}
}

task merge {
	input {
		String cohort_id
		Array[File] preprocessed_rds_objects

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(preprocessed_rds_objects, "GB") * 2 + 20)
	Int disk_size = ceil(size(preprocessed_rds_objects, "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/merge_rds.R \
			--paths-input ~{sep=' ' preprocessed_rds_objects} \
			--output ~{cohort_id}.merged_rds_object.rds

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.merged_rds_object.rds"
	>>>

	output {
		String merged_rds_object = "~{raw_data_path}/~{cohort_id}.merged_rds_object.rds"
	}

	runtime {
		docker: "~{container_registry}/spatial_r:1:0:0"
		cpu: 2
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 30
		zones: zones
	}
}

task process {
	input {
		String cohort_id
		File merged_rds_object

		File cell_type_markers_list

		Int min_genes_detected_in_percent_segment

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int mem_gb = ceil(size([merged_rds_object, cell_type_markers_list], "GB") * 2 + 20)
	Int disk_size = ceil(size([merged_rds_object, cell_type_markers_list], "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/process.R \
			--cohort-id ~{cohort_id} \
			--input ~{merged_rds_object} \
			--celltype-markers ~{cell_type_markers_list} \
			--min-segment ~{min_genes_detected_in_percent_segment} \
			--output ~{cohort_id}.processed.rds

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.segment_gene_detection_plot.png" \
			-o "~{cohort_id}.gene_detection_rate.csv" \
			-o "~{cohort_id}.q3_negprobe_plot.png" \
			-o "~{cohort_id}.normalization_plot.png"
	>>>

	output {
		File processed_rds_object = "~{cohort_id}.processed.rds"
		String segment_gene_detection_plot_png = "~{raw_data_path}/~{cohort_id}.segment_gene_detection_plot.png"
		String gene_detection_rate_csv = "~{raw_data_path}/~{cohort_id}.gene_detection_rate.csv"
		String q3_negprobe_plot_png = "~{raw_data_path}/~{cohort_id}.q3_negprobe_plot.png"
		String normalization_plot_png = "~{raw_data_path}/~{cohort_id}.normalization_plot.png.png"
	}

	runtime {
		docker: "~{container_registry}/spatial_r:1:0:0"
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
		String cohort_id
		File processed_rds_object

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(processed_rds_object, "GB") * 2 + 20)
	Int disk_size = ceil(size(processed_rds_object, "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		Rscript /opt/scripts/rds_to_adata.R \
			--input ~{processed_rds_object} \
			--output ~{cohort_id}.processed.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.processed.h5ad"
	>>>

	output {
		String processed_adata_object = "~{raw_data_path}/~{cohort_id}.processed.h5ad"
	}

	runtime {
		docker: "~{container_registry}/spatial_r:1:0:0"
		cpu: 2
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 30
		zones: zones
	}
}
