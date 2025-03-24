version 1.0

# Perform dataset integration, UMAP, and clustering steps

workflow integrate_data {
	input {
		String cohort_id
		File processed_adata_object

		Int n_comps
		String batch_key
		Float leiden_resolution

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	call integrate_sample_data {
		input:
			cohort_id = cohort_id,
			processed_adata_object = processed_adata_object,
			batch_key = batch_key,
			container_registry = container_registry,
			zones = zones
	}

	call cluster {
		input:
			cohort_id = cohort_id,
			integrated_adata_object = integrate_sample_data.integrated_adata_object,
			n_comps = n_comps,
			leiden_resolution = leiden_resolution,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	output {
		File integrated_adata_object = integrate_sample_data.integrated_adata_object
		File clustered_adata_object = cluster.clustered_adata_object #!FileCoercion
		File umap_cluster_plots_png = cluster.umap_cluster_plots_png #!FileCoercion
	}
}

task integrate_sample_data {
	input {
		String cohort_id
		File processed_adata_object

		String batch_key

		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(processed_adata_object, "GB") * 5 + 20)
	Int disk_size = ceil(size(processed_adata_object, "GB") * 3 + 50)

	command <<<
		set -euo pipefail

		python3 /opt/scripts/integrate_harmony.py \
			--adata-input ~{processed_adata_object} \
			--batch-key ~{batch_key} \
			--adata-output ~{cohort_id}.harmony_integrated.h5ad
	>>>

	output {
		File integrated_adata_object = "~{cohort_id}.harmony_integrated.h5ad"
	}

	runtime {
		docker: "~{container_registry}/spatial_py:1.0.0"
		cpu: 4
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
		File integrated_adata_object

		Int n_comps
		Float leiden_resolution

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(integrated_adata_object, "GB") * 2 + 20)
	Int disk_size = ceil(size(integrated_adata_object, "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		python3 /opt/scripts/cluster.py \
			--adata-input ~{integrated_adata_object} \
			--n-comps ~{n_comps} \
			--resolution ~{leiden_resolution} \
			--plots-prefix ~{cohort_id} \
			--adata-output ~{cohort_id}.clustered.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.clustered.h5ad" \
			-o "~{cohort_id}.umap_cluster.png"
	>>>

	output {
		String clustered_adata_object = "~{raw_data_path}/~{cohort_id}.clustered.h5ad"
		String umap_cluster_plots_png = "~{raw_data_path}/~{cohort_id}.umap_cluster.png"
	}

	runtime {
		docker: "~{container_registry}/spatial_py:1.0.0"
		cpu: 4
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 30
		zones: zones
	}
}
