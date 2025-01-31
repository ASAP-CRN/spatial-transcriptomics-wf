version 1.0

# Merge and process adata object with QC, filtering, normalization, and clustering

import "../../../wf-common/wdl/tasks/write_cohort_sample_list.wdl" as WriteCohortSampleList

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

		# Feature selection parameters
		String batch_key
		Int n_top_genes

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

	output {
		File cohort_sample_list = write_cohort_sample_list.cohort_sample_list #!FileCoercion

		# Merged adata objects and QC plots
		File merged_adata_object = merge_and_plot_qc_metrics.merged_adata_object #!FileCoercion
		Array[File] qc_plots_png = merge_and_plot_qc_metrics.qc_plots_png #!FileCoercion

		# Filtered and normalized adata object
		File filtered_normalized_adata_object = filter_and_normalize.filtered_normalized_adata_object #!FileCoercion

		# Feature selection

		# Leiden clustered adata object and UMAP and spatial coordinates plots
		File umap_cluster_adata_object = cluster.umap_cluster_adata_object #!FileCoercion
		Array[File] umap_and_spatial_coord_plots_png = cluster.umap_and_spatial_coord_plots_png #!FileCoercion
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
			--adata-output ~{cohort_id}..h5ad

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
		zones: zones
	}
}
