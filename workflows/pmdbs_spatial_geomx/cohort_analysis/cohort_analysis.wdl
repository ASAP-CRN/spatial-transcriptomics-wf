version 1.0

# Merge and process RDS object with QC, filtering, and normalization, and convert to adata object, integrate and cluster

import "../../../wf-common/wdl/tasks/write_cohort_sample_list.wdl" as WriteCohortSampleList
import "../../integrate_data/integrate_data.wdl" as IntegrateData
import "../../../wf-common/wdl/tasks/upload_final_outputs.wdl" as UploadFinalOutputs

workflow cohort_analysis {
	input {
		String cohort_id
		Array[Array[String]] project_sample_ids
		Array[File] preprocessed_rds_objects

		# If provided, these files will be uploaded to the staging bucket alongside other intermediate files made by this workflow
		Array[String] preprocessing_output_file_paths = []

		# Filtering parameters
		File cell_type_markers_list
		Float min_genes_detected_in_percent_segment

		# Integrate and cluster parameters
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

	call merge_and_prep {
		input:
			cohort_id = cohort_id,
			processed_adata_objects = rds_to_adata.processed_adata_object,
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
			processed_adata_object = merge_and_prep.merged_adata_object, #!FileCoercion
			n_comps = n_comps,
			batch_key = batch_key,
			leiden_resolution = leiden_resolution,
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
		filter_and_normalize.segment_gene_detection_plot_png,
		filter_and_normalize.gene_detection_rate_csv,
		filter_and_normalize.q3_negprobe_plot_png,
		filter_and_normalize.normalization_plot_png,
		[
			merge_and_prep.merged_adata_object,
			merge_and_prep.hvg_plot_png,
			merge_and_prep.merged_adata_metadata_csv
		],
		[
			integrate_data.clustered_adata_object,
			integrate_data.umap_cluster_plots_png
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

		# Processed RDS object and plots
		Array[File] processed_rds_object = filter_and_normalize.processed_rds_object
		Array[File] segment_gene_detection_plot_png = filter_and_normalize.segment_gene_detection_plot_png #!FileCoercion
		Array[File] gene_detection_rate_csv = filter_and_normalize.gene_detection_rate_csv #!FileCoercion
		Array[File] q3_negprobe_plot_png = filter_and_normalize.q3_negprobe_plot_png #!FileCoercion
		Array[File] normalization_plot_png = filter_and_normalize.normalization_plot_png #!FileCoercion

		# Converted AnnData object
		Array[File] processed_adata_object = rds_to_adata.processed_adata_object

		# Merged and prepped AnnData object
		File merged_adata_object = merge_and_prep.merged_adata_object #!FileCoercion
		File hvg_plot_png = merge_and_prep.hvg_plot_png #!FileCoercion
		File merged_adata_metadata_csv = merge_and_prep.merged_adata_metadata_csv #!FileCoercion

		# Integrate data outputs
		File integrated_adata_object = integrate_data.integrated_adata_object
		File clustered_adata_object = integrate_data.clustered_adata_object
		File umap_cluster_plots_png = integrate_data.umap_cluster_plots_png

		Array[File] preprocess_manifest_tsvs = upload_preprocess_files.manifests #!FileCoercion
		Array[File] cohort_analysis_manifest_tsvs = upload_cohort_analysis_files.manifests #!FileCoercion
	}

	meta {
		description: "Run team-level and/or cross-team cohort analysis on the Nanostring GeoMx data by filtering, normalization, dimensionality reduction, sample integration, and clustering."
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files."}
		project_sample_ids: {help: "Associated team ID and sample ID; used to generate a sample list."}
		preprocessed_rds_objects: {help: "An array of preprocessed RDS objects to run cohort analysis on."}
		preprocessing_output_file_paths: {help: "Selected preprocessed output files to upload to the staging bucket alongside selected cohort analysis output files."}
		cell_type_markers_list: {help: "CSV file containing a list of major cell type markers; used for detecting genes of interest."}
		min_genes_detected_in_percent_segment: {help: "Minimum % of segments that detect the genes. [0.01]"}
		n_top_genes: {help: "Number of highly-variable genes to keep. [3000]"}
		n_comps: {help: "Number of principal components to compute. [30]"}
		batch_key: {help: "Key in AnnData object for batch information. ['batch_id']"}
		leiden_resolution: {help: "Value controlling the coarseness of the Leiden clustering. [0.4]"}
		workflow_name: {help: "Workflow name; stored in the file-level manifest and final manifest with all saved files."}
		workflow_version: {help: "Workflow version; stored in the file-level manifest and final manifest with all saved files."}
		workflow_release: {help: "GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		run_timestamp: {help: "UTC timestamp; stored in the file-level manifest and final manifest with all saved files."}
		raw_data_path_prefix: {help: "Raw data bucket path prefix; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/cohort_analysis`)."}
		staging_data_buckets: {help: "Array of staging data buckets to upload intermediate files to (i.e., DEV or UAT buckets depending on internal QC status)."}
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
		maxRetries: 2
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
		maxRetries: 2
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

task merge_and_prep {
	input {
		String cohort_id
		Array[File] processed_adata_objects

		Int n_top_genes
		Int n_comps

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int mem_gb = ceil(size(processed_adata_objects, "GB") * 2 + 20)
	Int disk_size = ceil(size(processed_adata_objects, "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		geomx_merge_and_prep \
			--adata-paths-input ~{sep=' ' processed_adata_objects} \
			--n-top-genes ~{n_top_genes} \
			--n-comps ~{n_comps} \
			--output-prefix ~{cohort_id} \
			--adata-output ~{cohort_id}.merged.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{cohort_id}.merged.h5ad" \
			-o "~{cohort_id}.hvg_dispersion.png" \
			-o "~{cohort_id}.merged_adata_metadata.csv"

	>>>

	output {
		String merged_adata_object = "~{raw_data_path}/~{cohort_id}.merged.h5ad"
		String hvg_plot_png = "~{raw_data_path}/~{cohort_id}.hvg_dispersion.png"
		String merged_adata_metadata_csv = "~{raw_data_path}/~{cohort_id}.merged_adata_metadata.csv"
	}

	runtime {
		docker: "~{container_registry}/spatial_py:1.0.0"
		cpu: 2
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		maxRetries: 2
		bootDiskSizeGb: 15
		zones: zones
	}

	meta {
		description: "Merge sample-level AnnData objects to a single cohort-level AnnData object, and prepare for downstream analysis by annotating highly-variable genes (HVG) and performing PCA."
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files."}
		processed_adata_objects: {help: "An array of processed AnnData object to merge."}
		n_top_genes: {help: "Number of highly-variable genes to keep. [3000]"}
		n_comps: {help: "Number of principal components to compute. [30]"}
		raw_data_path: {help: "Raw data bucket path for merged adata and HVG plot outputs; location of raw bucket to upload task outputs to (`<raw_data_bucket>/workflow_execution/cohort_analysis/<cohort_analysis_version>/<run_timestamp>`)."}
		workflow_info: {help: "UTC timestamp, workflow name, workflow version, and GitHub release; stored in the file-level manifest and final manifest with all saved files."}
		billing_project: {help: "Billing project to charge GCP costs."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones where compute will take place. ['us-central1-c us-central1-f']"}
	}
}
