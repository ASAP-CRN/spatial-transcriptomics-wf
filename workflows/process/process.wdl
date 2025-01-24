version 1.0

# Process adata object with QC, filtering, and normalization

workflow process {
	input {
		String sample_id
		File preprocessed_adata_object

		String workflow_name
		String workflow_version
		String workflow_release
		String run_timestamp
		String raw_data_path_prefix
		String billing_project
		String container_registry
		String zones
	}

	String sub_workflow_name = "process"
	String sub_workflow_version = "1.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version, workflow_release]]

	String raw_data_path = "~{raw_data_path_prefix}/~{sub_workflow_name}/~{sub_workflow_version}/~{run_timestamp}"

	call qc_metrics {
		input:
			sample_id = sample_id,
			preprocessed_adata_object = preprocessed_adata_object,
			raw_data_path = raw_data_path,
			workflow_info = workflow_info,
			billing_project = billing_project,
			container_registry = container_registry,
			zones = zones
	}

	output {
		# QC adata object and plots
		File qc_adata_object = qc_metrics.qc_adata_object #!FileCoercion
		File qc_plots_png = qc_metrics.qc_plots_png #!FileCoercion
		Float qc_unassigned_ctrl_probes_percentage = qc_metrics.qc_unassigned_ctrl_probes_percentage
	}
}

task qc_metrics {
	input {
		String sample_id
		File preprocessed_adata_object

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(preprocessed_adata_object, "GB") * 2 + 20)
	Int disk_size = ceil(size(preprocessed_adata_object, "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		python3 /opt/scripts/qc.py \
			--sample-id ~{sample_id} \
			--adata-input ~{preprocessed_adata_object} \
			--qc-plots-output ~{sample_id}.qc.png \
			--qc-adata-output ~{sample_id}.qc.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{sample_id}.qc.png" \
			-o "~{sample_id}.qc.h5ad"
	>>>

	output {
		String qc_adata_object = "~{raw_data_path}/~{sample_id}.qc.h5ad"
		String qc_plots_png = "~{raw_data_path}/~{sample_id}.qc.png"
		Float qc_unassigned_ctrl_probes_percentage = read_float("unassigned_ctrl_probes_percentage.txt")
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
