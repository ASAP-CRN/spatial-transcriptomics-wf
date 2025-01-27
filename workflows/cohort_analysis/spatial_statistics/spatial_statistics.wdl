version 1.0

# Run spatial statistics

workflow spatial_statistics {
	input {
		String cohort_id
		File umap_cluster_adata_object

		String workflow_name
		String workflow_version
		String workflow_release
		String run_timestamp
		String raw_data_path_prefix
		String billing_project
		String container_registry
		String zones
	}

	String sub_workflow_name = "spatial_statistics"
	String sub_workflow_version = "1.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version, workflow_release]]

	String raw_data_path = "~{raw_data_path_prefix}/~{sub_workflow_name}/~{sub_workflow_version}/~{run_timestamp}"

	call neighbors_enrichment_analysis {
		input:
			cohort_id = cohort_id,
			umap_cluster_adata_object = umap_cluster_adata_object,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	output {
		# Neighborhood enrichment adata object and plot
		File nhood_enrichment_adata_object = neighbors_enrichment_analysis.nhood_enrichment_adata_object #!FileCoercion
		File nhood_enrichment_plot_png = neighbors_enrichment_analysis.nhood_enrichment_plot_png #!FileCoercion
	}
}

task neighbors_enrichment_analysis {
	input {
		String cohort_id
		File umap_cluster_adata_object

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(umap_cluster_adata_object, "GB") * 2 + 20)
	Int disk_size = ceil(size(umap_cluster_adata_object, "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		python3 /opt/scripts/neighbors_enrichment_analysis.py \
			--cohort-id ~{cohort_id} \
			--adata-input ~{umap_cluster_adata_object} \
			--adata-output ~{cohort_id}.nhood_enrichment_adata_object.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.nhood_enrichment_adata_object.h5ad" \
			-o "~{cohort_id}.nhood_enrichment.png"
	>>>

	output {
		String nhood_enrichment_adata_object = "~{raw_data_path}/~{cohort_id}.nhood_enrichment_adata_object.h5ad"
		String nhood_enrichment_plot_png = "~{raw_data_path}/~{cohort_id}.nhood_enrichment.png"
	}

	runtime {
		docker: "~{container_registry}/squidpy:1.6.2"
		cpu: 2
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		zones: zones
	}
}
