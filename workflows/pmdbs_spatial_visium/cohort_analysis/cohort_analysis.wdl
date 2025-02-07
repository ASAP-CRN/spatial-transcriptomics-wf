version 1.0

# Merge and process adata object with QC, filtering, normalization, and clustering

import "../../../wf-common/wdl/tasks/write_cohort_sample_list.wdl" as WriteCohortSampleList
import "pmdbs-sc-rnaseq-wf/workflows/cohort_analysis/cluster_data/cluster_data.wdl" as ClusterData
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

		# Feature selection and clustering inputs and parameters
		String batch_key
		String scvi_latent_key
		Int n_top_genes
		File cell_type_markers_list

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
			container_registry = container_registry,
			zones = zones
	}

	call feature_selection {
		input:
			cohort_id = cohort_id,
			batch_key = batch_key,
			n_top_genes = n_top_genes,
			filtered_normalized_adata_object = filter_and_normalize.filtered_normalized_adata_object, #!FileCoercion
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call ClusterData.cluster_data {
		input:
			cohort_id = cohort_id,
			normalized_adata_object = feature_selection.feature_selection_adata_object, #!FileCoercion
			scvi_latent_key = scvi_latent_key,
			batch_key = batch_key,
			cell_type_markers_list = cell_type_markers_list,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call image_features {
		input:
			cohort_id = cohort_id,
			cell_annotated_adata_object = cluster_data.cell_annotated_adata_object, #!FileCoercion
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call SpatialStatistics.spatial_statistics {
		input:
			cohort_id = cohort_id,
			clustered_adata_object = image_features.image_features_adata_object, #!FileCoercion
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
			feature_selection.feature_dispersion_plot_png
		],
		[
			cluster_data.scvi_model_tar_gz,
			cluster_data.cell_types_csv
		],
		[
			image_features.image_features_spatial_scatter_plot_png
		],
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
		Array[File] qc_plots_png = merge_and_plot_qc_metrics.qc_plots_png #!FileCoercion

		# Filtered and normalized adata object
		File filtered_normalized_adata_object = filter_and_normalize.filtered_normalized_adata_object #!FileCoercion

		# Feature selection
		File feature_selection_adata_object = feature_selection.feature_selection_adata_object #!FileCoercion
		File feature_dispersion_plot_png = feature_selection.feature_dispersion_plot_png #!FileCoercion

		# Clustering output
		File integrated_adata_object = cluster_data.integrated_adata_object
		File scvi_model_tar_gz = cluster_data.scvi_model_tar_gz
		File umap_cluster_adata_object = cluster_data.umap_cluster_adata_object
		File cell_annotated_adata_object = cluster_data.cell_annotated_adata_object
		File cell_types_csv = cluster_data.cell_types_csv

		# Image features output
		File image_features_adata_object = image_features.image_features_adata_object #!FileCoercion
		File image_features_spatial_scatter_plot_png = image_features.image_features_spatial_scatter_plot_png #!FileCoercion

		# Spatial statistics output
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

		python3 /opt/scripts/merge_and_plot_visium_qc.py \
			--adata-paths-input ~{sep=' ' preprocessed_adata_objects} \
			--qc-plots-prefix ~{cohort_id} \
			--merged-adata-output ~{cohort_id}.merged_adata_object.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.merged_adata_object.h5ad" \
			-o "~{cohort_id}.qc_violin.png" \
			-o "~{cohort_id}.qc_scatter.png"
	>>>

	output {
		String merged_adata_object = "~{raw_data_path}/~{cohort_id}.merged_adata_object.h5ad"
		Array[String] qc_plots_png = [
			"~{raw_data_path}/~{cohort_id}.qc_violin.png",
			"~{raw_data_path}/~{cohort_id}.qc_scatter.png"
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

task filter_and_normalize {
	input {
		String cohort_id
		File merged_adata_object

		Int filter_cells_min_counts
		Int filter_genes_min_cells

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
	>>>

	output {
		File filtered_normalized_adata_object = "~{cohort_id}.filtered_normalized.h5ad"
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

task feature_selection {
	input {
		String cohort_id
		File filtered_normalized_adata_object

		String batch_key
		Int n_top_genes

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

		python3 /opt/scripts/feature_selection.py \
			--adata-input ~{filtered_normalized_adata_object} \
			--batch-key ~{batch_key} \
			--n-top-genes ~{n_top_genes} \
			--plots-prefix ~{cohort_id} \
			--adata-output ~{cohort_id}.hvg_pca_neighbors_umap.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.feature_dispersion.png"
	>>>

	output {
		File feature_selection_adata_object = "~{cohort_id}.hvg_pca_neighbors_umap.h5ad"
		String feature_dispersion_plot_png = "~{raw_data_path}/~{cohort_id}.umap.png"
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

task image_features {
	input {
		String cohort_id
		File cell_annotated_adata_object

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 4
	Int mem_gb = ceil(size(cell_annotated_adata_object, "GB") * 2 + 20)
	Int disk_size = ceil(size(cell_annotated_adata_object, "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		python3 /opt/scripts/image_features.py \
			--adata-input ~{cell_annotated_adata_object} \
			--n-jobs ~{threads} \
			--plots-prefix ~{cohort_id} \
			--adata-output ~{cohort_id}.image_features.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.image_features_spatial_scatter.png"
	>>>

	output {
		File image_features_adata_object = "~{cohort_id}.image_features.h5ad"
		String image_features_spatial_scatter_plot_png = "~{raw_data_path}/~{cohort_id}.image_features_spatial_scatter.png"
	}

	runtime {
		docker: "~{container_registry}/squidpy:1.6.2_1"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 30
		zones: zones
	}
}
