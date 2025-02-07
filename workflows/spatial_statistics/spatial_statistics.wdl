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

	call neighbors_enrichment_analysis {
		input:
			cohort_id = cohort_id,
			clustered_adata_object = clustered_adata_object,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call co_occurrence_probability {
		input:
			cohort_id = cohort_id,
			nhood_enrichment_adata_object = neighbors_enrichment_analysis.nhood_enrichment_adata_object,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	call spatially_variable_gene_analysis {
		input:
			cohort_id = cohort_id,
			co_occurrence_adata_object = co_occurrence_probability.co_occurrence_adata_object,
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

		# Co-occurrence probability adata object and plot
		File co_occurrence_adata_object = co_occurrence_probability.co_occurrence_adata_object #!FileCoercion
		File co_occurrence_plot_png = co_occurrence_probability.co_occurrence_plot_png #!FileCoercion

		# Moran’s I global spatial auto-correlation statistics
		File final_adata_object = spatially_variable_gene_analysis.final_adata_object #!FileCoercion
		File moran_top_10_variable_genes_csv = spatially_variable_gene_analysis.moran_top_10_variable_genes_csv #!FileCoercion
	}
}

task neighbors_enrichment_analysis {
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

		python3 /opt/scripts/neighbors_enrichment_analysis.py \
			--cohort-id ~{cohort_id} \
			--adata-input ~{clustered_adata_object} \
			--adata-output ~{cohort_id}.nhood_enrichment_adata_object.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.nhood_enrichment.png"
	>>>

	output {
		File nhood_enrichment_adata_object = "~{cohort_id}.nhood_enrichment_adata_object.h5ad"
		String nhood_enrichment_plot_png = "~{raw_data_path}/~{cohort_id}.nhood_enrichment.png"
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

task co_occurrence_probability {
	input {
		String cohort_id
		File nhood_enrichment_adata_object

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(nhood_enrichment_adata_object, "GB") * 2 + 20)
	Int disk_size = ceil(size(nhood_enrichment_adata_object, "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		python3 /opt/scripts/co_occurrence_probability.py \
			--cohort-id ~{cohort_id} \
			--adata-input ~{nhood_enrichment_adata_object} \
			--adata-output ~{cohort_id}.co_occurrence.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.co_occurrence.png"
	>>>

	output {
		File co_occurrence_adata_object = "~{cohort_id}.co_occurrence.h5ad"
		String co_occurrence_plot_png = "~{raw_data_path}/~{cohort_id}.co_occurrence.png"
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

task spatially_variable_gene_analysis {
	input {
		String cohort_id
		File co_occurrence_adata_object

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(co_occurrence_adata_object, "GB") * 2 + 20)
	Int disk_size = ceil(size(co_occurrence_adata_object, "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		python3 /opt/scripts/moran_i_score.py \
			--cohort-id ~{cohort_id} \
			--adata-input ~{co_occurrence_adata_object} \
			--adata-output ~{cohort_id}.final_adata_object.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.final_adata_object.h5ad" \
			-o "~{cohort_id}.moran_top_10_variable_genes.csv"
	>>>

	output {
		String final_adata_object = "~{raw_data_path}/~{cohort_id}.final_adata_object.h5ad"
		String moran_top_10_variable_genes_csv = "~{raw_data_path}/~{cohort_id}.moran_top_10_variable_genes.csv"
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
