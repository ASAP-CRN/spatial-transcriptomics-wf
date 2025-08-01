version 1.0

# Process RDS object with QC, filtering, and normalization, and convert to adata object

workflow process_to_adata {
	input {
		Array[File] preprocessed_rds_objects

		# Filtering parameters
		File cell_type_markers_list
		Float min_genes_detected_in_percent_segment

		String workflow_name
		String workflow_version
		String workflow_release
		String run_timestamp
		String raw_data_path_prefix
		String billing_project
		String container_registry
		String zones
	}

	String sub_workflow_name = "process_to_adata"
	String sub_workflow_version = "1.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version, workflow_release]]

	String raw_data_path = "~{raw_data_path_prefix}/~{sub_workflow_name}/~{sub_workflow_version}/~{run_timestamp}"

	scatter (preprocessed_rds_object in preprocessed_rds_objects) {
		call filter_and_normalize {
			input:
				preprocessed_rds_object = preprocessed_rds_object,
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
				processed_rds_object = filter_and_normalize.processed_rds_object,
				container_registry = container_registry,
				zones = zones
		}
	}

	output {
		# Processed RDS object and plots
		Array[File] processed_rds_objects = filter_and_normalize.processed_rds_object
		Array[File] segment_gene_detection_plot_png = filter_and_normalize.segment_gene_detection_plot_png #!FileCoercion
		Array[File] gene_detection_rate_csv = filter_and_normalize.gene_detection_rate_csv #!FileCoercion
		Array[File] q3_negprobe_plot_png = filter_and_normalize.q3_negprobe_plot_png #!FileCoercion
		Array[File] normalization_plot_png = filter_and_normalize.normalization_plot_png #!FileCoercion

		# Converted AnnData object
		Array[File] processed_adata_objects = rds_to_adata.processed_adata_object
	}

	meta {
		description: "Process the Nanostring GeoMx data by filtering, normalization, and convert to AnnData object."
	}

	parameter_meta {
		preprocessed_rds_objects: {help: "An array of preprocessed RDS objects to run cohort analysis on."}
		cell_type_markers_list: {help: "CSV file containing a list of major cell type markers; used for detecting genes of interest."}
		min_genes_detected_in_percent_segment: {help: "Minimum % of segments that detect the genes. [0.01]"}
		workflow_name: {help: "Workflow name; stored in the file-level manifest and final manifest with all saved files."}
		workflow_version: {help: "Workflow version; stored in the file-level manifest and final manifest with all saved files."}
		workflow_release: {help: "GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		run_timestamp: {help: "UTC timestamp; stored in the file-level manifest and final manifest with all saved files."}
		raw_data_path_prefix: {help: "Raw data bucket path prefix; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/cohort_analysis`)."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones where compute will take place. ['us-central1-c us-central1-f']"}
	}
}

task filter_and_normalize {
	input {
		File preprocessed_rds_object

		File cell_type_markers_list

		Float min_genes_detected_in_percent_segment

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	String slide_id = basename(preprocessed_rds_object, ".qc.rds")

	Int mem_gb = ceil(size([preprocessed_rds_object, cell_type_markers_list], "GB") * 2 + 20)
	Int disk_size = ceil(size([preprocessed_rds_object, cell_type_markers_list], "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		# Select ROI/AOI segments and genes based on LOQ and normalization
		geomx_process \
			--slide-id ~{slide_id} \
			--input ~{preprocessed_rds_object} \
			--celltype-markers ~{cell_type_markers_list} \
			--min-segment ~{min_genes_detected_in_percent_segment} \
			--output ~{slide_id}.processed.rds

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{slide_id}.segment_gene_detection_plot.png" \
			-o "~{slide_id}.gene_detection_rate.csv" \
			-o "~{slide_id}.q3_negprobe_plot.png" \
			-o "~{slide_id}.normalization_plot.png"
	>>>

	output {
		File processed_rds_object = "~{slide_id}.processed.rds"
		String segment_gene_detection_plot_png = "~{raw_data_path}/~{slide_id}.segment_gene_detection_plot.png"
		String gene_detection_rate_csv = "~{raw_data_path}/~{slide_id}.gene_detection_rate.csv"
		String q3_negprobe_plot_png = "~{raw_data_path}/~{slide_id}.q3_negprobe_plot.png"
		String normalization_plot_png = "~{raw_data_path}/~{slide_id}.normalization_plot.png"
	}

	runtime {
		docker: "~{container_registry}/spatial_r:1.0.0"
		cpu: 2
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 15
		zones: zones
	}

	meta {
		description: "Perform Limit of Quantification (LOQ), filter segments and/or genes with low signal, and normalize (Q3 and background)."
	}

	parameter_meta {
		preprocessed_rds_object: {help: "Preprocessed RDS object to run cohort analysis on."}
		cell_type_markers_list: {help: "CSV file containing a list of major cell type markers; used for detecting genes of interest."}
		min_genes_detected_in_percent_segment: {help: "Minimum % of segments that detect the genes. [0.01]"}
		raw_data_path: {help: "Raw data bucket path for processed RDS and gene detection and normalization plots outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/cohort_analysis/<cohort_analysis_version>/<run_timestamp>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones where compute will take place. ['us-central1-c us-central1-f']"}
	}
}

task rds_to_adata {
	input {
		File processed_rds_object

		String container_registry
		String zones
	}

	String slide_id = basename(processed_rds_object, ".processed.rds")

	Int mem_gb = ceil(size(processed_rds_object, "GB") * 2 + 20)
	Int disk_size = ceil(size(processed_rds_object, "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		geomx_rds_to_adata \
			--input ~{processed_rds_object} \
			--output-prefix ~{slide_id}.processed
	>>>

	output {
		File processed_adata_object = "~{slide_id}.processed.h5ad"
	}

	runtime {
		docker: "~{container_registry}/spatial_r:1.0.0"
		cpu: 2
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 15
		zones: zones
	}

	meta {
		description: "Convert processed NanoStringGeoMxSet (.RDS) object to AnnData object with Seurat."
	}

	parameter_meta {
		processed_rds_object: {help: "Processed RDS object to convert into an AnnData object."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones where compute will take place. ['us-central1-c us-central1-f']"}
	}
}
