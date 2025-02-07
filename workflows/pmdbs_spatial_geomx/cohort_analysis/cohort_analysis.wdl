version 1.0

# Merge and process adata object with QC, filtering, normalization, and clustering

import "../../../wf-common/wdl/tasks/write_cohort_sample_list.wdl" as WriteCohortSampleList
import "../../spatial_statistics/spatial_statistics.wdl" as SpatialStatistics
import "../../../wf-common/wdl/tasks/upload_final_outputs.wdl" as UploadFinalOutputs

workflow cohort_analysis {
	input {
		String cohort_id
		Array[Array[String]] project_sample_ids
		Array[File] preprocessed_adata_objects

		# If provided, these files will be uploaded to the staging bucket alongside other intermediate files made by this workflow
		Array[String] preprocessing_output_file_paths = []

		# Filter parameters
		Int filter_cells_min_counts
		Int filter_genes_min_cells

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
			filter_genes_min_cells = filter_genes_min_cells,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call cluster {
		input:
			cohort_id = cohort_id,
			filtered_normalized_adata_object = filter_and_normalize.filtered_normalized_adata_object, #!FileCoercion
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
		cluster.umap_and_spatial_coord_plots_png,
		[
			spatial_statistics.nhood_enrichment_plot_png,
			spatial_statistics.co_occurrence_plot_png,
			spatial_statistics.final_adata_object,
			spatial_statistics.moran_top_10_variable_genes_csv
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
		File qc_plots_png = merge_and_plot_qc_metrics.qc_plots_png #!FileCoercion
		Float qc_unassigned_ctrl_probes_percentage = merge_and_plot_qc_metrics.qc_unassigned_ctrl_probes_percentage

		# Filtered and normalized adata object
		File filtered_normalized_adata_object = filter_and_normalize.filtered_normalized_adata_object #!FileCoercion

		# Leiden clustered adata object and UMAP and spatial coordinates plots
		File umap_cluster_adata_object = cluster.umap_cluster_adata_object #!FileCoercion
		Array[File] umap_and_spatial_coord_plots_png = cluster.umap_and_spatial_coord_plots_png #!FileCoercion

		# Spatial statistics outputs
		File nhood_enrichment_adata_object = spatial_statistics.nhood_enrichment_adata_object
		File nhood_enrichment_plot_png = spatial_statistics.nhood_enrichment_plot_png
		File co_occurrence_adata_object = spatial_statistics.co_occurrence_adata_object
		File co_occurrence_plot_png = spatial_statistics.co_occurrence_plot_png
		File final_adata_object = spatial_statistics.final_adata_object
		File moran_top_10_variable_genes_csv = spatial_statistics.moran_top_10_variable_genes_csv

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

		python3 /opt/scripts/merge_and_plot_geomx_qc.py \
			--adata-paths-input ~{sep=' ' preprocessed_adata_objects} \
			--merged-adata-output ~{cohort_id}.merged_adata_object.h5ad \
			--qc-plots-output ~{cohort_id}.qc_hist.png

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.merged_adata_object.h5ad" \
			-o "~{cohort_id}.qc_hist.png"
	>>>

	output {
		String merged_adata_object = "~{raw_data_path}/~{cohort_id}.merged_adata_object.h5ad"
		String qc_plots_png = "~{raw_data_path}/~{cohort_id}.qc_hist.png"
		Float qc_unassigned_ctrl_probes_percentage = read_float("unassigned_ctrl_probes_percentage.txt")
	}

	runtime {
		docker: "~{container_registry}/squidpy:1.6.2_1"
		cpu: 2
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 30
		zones: zones
	}
}

task filter_and_normalize {
	input {
		String cohort_id
		File merged_adata_object

		Int filter_cells_min_counts
		Int filter_genes_min_cells

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

		python3 /opt/scripts/filter_and_normalize.py \
			--adata-input ~{merged_adata_object} \
			--min-counts ~{filter_cells_min_counts} \
			--min-cells ~{filter_genes_min_cells} \
			--adata-output ~{cohort_id}.filtered_normalized.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.filtered_normalized.h5ad"
	>>>

	output {
		String filtered_normalized_adata_object = "~{raw_data_path}/~{cohort_id}.filtered_normalized.h5ad"
	}

	runtime {
		docker: "~{container_registry}/squidpy:1.6.2_1"
		cpu: 2
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 30
		zones: zones
	}
}

task cluster {
	input {
		String cohort_id
		File filtered_normalized_adata_object

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(filtered_normalized_adata_object, "GB") * 2 + 20)
	Int disk_size = ceil(size(filtered_normalized_adata_object, "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		python3 /opt/scripts/cluster.py \
			--cohort-id ~{cohort_id} \
			--adata-input ~{filtered_normalized_adata_object} \
			--adata-output ~{cohort_id}.umap_cluster.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.umap_cluster.h5ad" \
			-o "~{cohort_id}.umap.png" \
			-o "~{cohort_id}.spatial_coord_by_counts.png" \
			-o "~{cohort_id}.spatial_coord_by_clusters.png"
	>>>

	output {
		String umap_cluster_adata_object = "~{raw_data_path}/~{cohort_id}.umap_cluster.h5ad"
		Array[String] umap_and_spatial_coord_plots_png = [
			"~{raw_data_path}/~{cohort_id}.umap.png",
			"~{raw_data_path}/~{cohort_id}.spatial_coord_by_counts.png",
			"~{raw_data_path}/~{cohort_id}.spatial_coord_by_clusters.png"
		]
	}

	runtime {
		docker: "~{container_registry}/squidpy:1.6.2_1"
		cpu: 2
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 30
		zones: zones
	}
}
