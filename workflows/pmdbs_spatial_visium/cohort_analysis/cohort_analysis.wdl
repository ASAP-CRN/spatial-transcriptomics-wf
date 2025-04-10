version 1.0

# Merge and process adata object with QC, filtering, normalization, dimensionality reduction, integration, and clustering

import "../../../wf-common/wdl/tasks/write_cohort_sample_list.wdl" as WriteCohortSampleList
import "../../integrate_data/integrate_data.wdl" as IntegrateData
import "spatial_statistics/spatial_statistics.wdl" as SpatialStatistics
import "../../../wf-common/wdl/tasks/upload_final_outputs.wdl" as UploadFinalOutputs

workflow cohort_analysis {
	input {
		String cohort_id
		Array[Array[String]] project_sample_ids
		Array[File] preprocessed_adata_objects

		# If provided, these files will be uploaded to the staging bucket alongside other intermediate files made by this workflow
		Array[String] preprocessing_output_file_paths = []

		# Processing parameters
		Int filter_cells_min_counts
		Int filter_cells_min_genes
		Int filter_genes_min_cells
		Float filter_mt_max_percent
		Float normalize_target_sum
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

	call merge_and_plot_qc_metrics {
		input:
			cohort_id = cohort_id,
			preprocessed_adata_objects = preprocessed_adata_objects,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call filter_and_normalize {
		input:
			cohort_id = cohort_id,
			merged_adata_object = merge_and_plot_qc_metrics.merged_adata_object, #!FileCoercion
			filter_cells_min_counts = filter_cells_min_counts,
			filter_cells_min_genes = filter_cells_min_genes,
			filter_genes_min_cells = filter_genes_min_cells,
			filter_mt_max_percent = filter_mt_max_percent,
			normalize_target_sum = normalize_target_sum,
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
			processed_adata_object = filter_and_normalize.processed_adata_object, #!FileCoercion
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
			clustered_adata_object = integrate_data.clustered_adata_object, #!FileCoercion
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
			merge_and_plot_qc_metrics.merged_adata_object
		],
		merge_and_plot_qc_metrics.qc_plots_png,
		[
			filter_and_normalize.hvg_plot_png
		],
		[
			integrate_data.clustered_adata_object,
			integrate_data.umap_cluster_plots_png
		],
		[
			plot_spatial.spatial_scatter_plot_png
		],
		[	
			spatial_statistics.final_adata_object,
			spatial_statistics.moran_top_10_variable_genes_csv,
			spatial_statistics.moran_top_3_variable_genes_spatial_scatter_plot_png
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

		# Merged adata objects and QC plots
		File merged_adata_object = merge_and_plot_qc_metrics.merged_adata_object #!FileCoercion
		Array[File] qc_plots_png = merge_and_plot_qc_metrics.qc_plots_png #!FileCoercion

		# Processed outputs
		File processed_adata_object = filter_and_normalize.processed_adata_object
		File hvg_plot_png = filter_and_normalize.hvg_plot_png #!FileCoercion

		# Integrate data outputs
		File integrated_adata_object = integrate_data.integrated_adata_object
		File clustered_adata_object = integrate_data.clustered_adata_object
		File umap_cluster_plots_png = integrate_data.umap_cluster_plots_png

		# Spatial plots
		File spatial_scatter_plot_png = plot_spatial.spatial_scatter_plot_png #!FileCoercion

		# Spatial statistics outputs
		File final_adata_object = spatial_statistics.final_adata_object
		File moran_top_10_variable_genes_csv = spatial_statistics.moran_top_10_variable_genes_csv
		File moran_top_3_variable_genes_spatial_scatter_plot_png = spatial_statistics.moran_top_3_variable_genes_spatial_scatter_plot_png

		Array[File] preprocess_manifest_tsvs = upload_preprocess_files.manifests #!FileCoercion
		Array[File] cohort_analysis_manifest_tsvs = upload_cohort_analysis_files.manifests #!FileCoercion
	}
}

task merge_and_plot_qc_metrics {
	input {
		String cohort_id
		Array[File] preprocessed_adata_objects

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(preprocessed_adata_objects, "GB") * 2 + 20)
	Int disk_size = ceil(size(preprocessed_adata_objects, "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		python3 /opt/scripts/visium_merge_and_plot_qc.py \
			--adata-paths-input ~{write_lines(preprocessed_adata_objects)} \
			--qc-plots-prefix ~{cohort_id} \
			--merged-adata-output ~{cohort_id}.merged_adata_object.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.merged_adata_object.h5ad" \
			-o "~{cohort_id}.qc_violin.png" \
			-o "~{cohort_id}.qc_dist.png"
	>>>

	output {
		String merged_adata_object = "~{raw_data_path}/~{cohort_id}.merged_adata_object.h5ad"
		Array[String] qc_plots_png = [
			"~{raw_data_path}/~{cohort_id}.qc_violin.png",
			"~{raw_data_path}/~{cohort_id}.qc_dist.png"
		]
	}

	runtime {
		docker: "~{container_registry}/spatial_py:1.0.0"
		cpu: 2
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 5
		zones: zones
	}
}

task filter_and_normalize {
	input {
		String cohort_id
		File merged_adata_object

		Int filter_cells_min_counts
		Int filter_cells_min_genes
		Int filter_genes_min_cells
		Float filter_mt_max_percent
		Float normalize_target_sum
		Int n_top_genes
		Int n_comps

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(merged_adata_object, "GB") * 2 + 20)
	Int disk_size = ceil(size(merged_adata_object, "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		python3 /opt/scripts/visium_process.py \
			--adata-input ~{merged_adata_object} \
			--min-counts ~{filter_cells_min_counts} \
			--min-genes ~{filter_cells_min_genes} \
			--min-cells ~{filter_genes_min_cells} \
			--mt-max-percent ~{filter_mt_max_percent} \
			--target-sum ~{normalize_target_sum} \
			--n-top-genes ~{n_top_genes} \
			--n-comps ~{n_comps} \
			--plots-prefix ~{cohort_id} \
			--adata-output ~{cohort_id}.processed.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.hvg_dispersion.png"
	>>>

	output {
		File processed_adata_object = "~{cohort_id}.processed.h5ad"
		String hvg_plot_png = "~{raw_data_path}/~{cohort_id}.hvg_dispersion.png"
	}

	runtime {
		docker: "~{container_registry}/spatial_py:1.0.0"
		cpu: 4
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 5
		zones: zones
	}
}

task plot_spatial {
	input {
		String cohort_id
		File clustered_adata_object

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(clustered_adata_object, "GB") * 2 + 20)
	Int disk_size = ceil(size(clustered_adata_object, "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		python3 /opt/scripts/visium_plot_spatial.py \
			--adata-input ~{clustered_adata_object} \
			--plots-prefix ~{cohort_id}

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.spatial_scatter.png"
	>>>

	output {
		String spatial_scatter_plot_png = "~{raw_data_path}/~{cohort_id}.spatial_scatter.png"
	}

	runtime {
		docker: "~{container_registry}/spatial_py:1.0.0"
		cpu: 2
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 5
		zones: zones
	}
}
