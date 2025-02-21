version 1.0

# Run spatial statistics

workflow spatial_statistics {
	input {
		String cohort_id
		File clustered_adata_object

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	call spatially_variable_gene_analysis {
		input:
			cohort_id = cohort_id,
			clustered_adata_object = clustered_adata_object,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	output {
		# Moran’s I global spatial auto-correlation statistics
		File final_adata_object = spatially_variable_gene_analysis.final_adata_object #!FileCoercion
		File moran_top_10_variable_genes_csv = spatially_variable_gene_analysis.moran_top_10_variable_genes_csv #!FileCoercion
		File moran_top_3_variable_genes_spatial_scatter_plot_png = spatially_variable_gene_analysis.moran_top_3_variable_genes_spatial_scatter_plot_png #!FileCoercion
	}
}

task spatially_variable_gene_analysis {
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

		python3 /opt/scripts/identify_spatially_variable_genes.py \
			--cohort-id ~{cohort_id} \
			--adata-input ~{clustered_adata_object} \
			--adata-output ~{cohort_id}.final_adata_object.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.final_adata_object.h5ad" \
			-o "~{cohort_id}.moran_top_10_variable_genes.csv" \
			-o "~{cohort_id}.moran_top_3_variable_genes_spatial_scatter.png"
	>>>

	output {
		String final_adata_object = "~{cohort_id}.final_adata_object.h5ad"
		String moran_top_10_variable_genes_csv = "~{raw_data_path}/~{cohort_id}.moran_top_10_variable_genes.csv"
		String moran_top_3_variable_genes_spatial_scatter_plot_png = "~{raw_data_path}/~{cohort_id}.moran_top_3_variable_genes_spatial_scatter.png"
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
