version 1.0

# Harmonized human and non-human PMDBS spatial transcriptomics workflow entrypoint for Nanostring GeoMx data

import "../../wf-common/wdl/structs.wdl"
import "../../wf-common/wdl/tasks/get_workflow_metadata.wdl" as GetWorkflowMetadata
import "preprocess/preprocess.wdl" as Preprocess
import "cohort_analysis/cohort_analysis.wdl" as CohortAnalysis

workflow pmdbs_spatial_geomx_analysis {
	input {
		String cohort_id
		Array[Project] projects

		File geomxngs_config_pkc

		# QC parameters
		Int min_segment_reads = 1000
		Int min_percent_reads_trimmed = 80
		Int min_percent_reads_stitched = 80
		Int min_percent_reads_aligned = 80
		Int min_saturation = 50
		Int min_neg_ctrl_count = 10
		Int max_ntc_count = 1000
		Int min_nuclei = 100
		Int min_segment_area = 5000

		# Filter parameters
		Int filter_cells_min_counts = 5000
		Int filter_genes_min_cells = 10

		# Cohort analysis
		Boolean run_cross_team_cohort_analysis = false
		String cohort_raw_data_bucket
		Array[String] cohort_staging_data_buckets

		String container_registry
		String zones = "us-central1-c us-central1-f"
	}

	String workflow_execution_path = "workflow_execution"
	String workflow_name = "pmdbs_spatial_geomx"
	String workflow_version = "v1.0.0"
	String workflow_release = "https://github.com/ASAP-CRN/pmdbs-spatial-transcriptomics-wf/releases/tag/pmdbs_spatial_geomx_analysis-~{workflow_version}"

	call GetWorkflowMetadata.get_workflow_metadata {
		input:
			zones = zones
	}

	scatter (project in projects) {
		String project_raw_data_path_prefix = "~{project.raw_data_bucket}/~{workflow_execution_path}/~{workflow_name}"

		call Preprocess.preprocess {
			input:
				team_id = project.team_id,
				dataset_id = project.dataset_id,
				samples = project.samples,
				project_sample_metadata_csv = project.project_sample_metadata_csv,
				geomx_config_ini = project.geomx_config_ini,
				geomx_lab_annotation_xlsx = project.geomx_lab_annotation_xlsx,
				geomxngs_config_pkc = geomxngs_config_pkc,
				min_segment_reads = min_segment_reads,
				min_percent_reads_trimmed = min_percent_reads_trimmed,
				min_percent_reads_stitched = min_percent_reads_stitched,
				min_percent_reads_aligned = min_percent_reads_aligned,
				min_saturation = min_saturation,
				min_neg_ctrl_count = min_neg_ctrl_count,
				max_ntc_count = max_ntc_count,
				min_nuclei = min_nuclei,
				min_segment_area = min_segment_area,
				workflow_name = workflow_name,
				workflow_version = workflow_version,
				workflow_release = workflow_release,
				run_timestamp = get_workflow_metadata.timestamp,
				raw_data_path_prefix = project_raw_data_path_prefix,
				billing_project = get_workflow_metadata.billing_project,
				container_registry = container_registry,
				zones = zones
		}

		Array[String] preprocessing_output_file_paths = flatten([
			preprocess.geomxngs_dcc_zip,
			preprocess.geomxngs_output_tar_gz,
			preprocess.initial_rds_object,
			preprocess.sample_overview_sankey_html,
			preprocess.qc_rds_object,
			preprocess.segment_qc_summary_csv,
			preprocess.probe_qc_summary_csv,
			preprocess.gene_count_csv
		]) #!StringCoercion

		if (project.run_project_cohort_analysis) {
			call CohortAnalysis.cohort_analysis as project_cohort_analysis {
				input:
					cohort_id = project.team_id,
					project_sample_ids = preprocess.project_sample_ids,
					preprocessed_rds_objects = preprocess.qc_rds_object,
					preprocessing_output_file_paths = preprocessing_output_file_paths,
					filter_cells_min_counts = filter_cells_min_counts,
					filter_genes_min_cells = filter_genes_min_cells,
					workflow_name = workflow_name,
					workflow_version = workflow_version,
					workflow_release = workflow_release,
					run_timestamp = get_workflow_metadata.timestamp,
					raw_data_path_prefix = project_raw_data_path_prefix,
					staging_data_buckets = project.staging_data_buckets,
					billing_project = get_workflow_metadata.billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}
	}

	if (run_cross_team_cohort_analysis) {
		String cohort_raw_data_path_prefix = "~{cohort_raw_data_bucket}/~{workflow_execution_path}/~{workflow_name}"

		call CohortAnalysis.cohort_analysis as cross_team_cohort_analysis {
			input:
				cohort_id = cohort_id,
				project_sample_ids = flatten(preprocess.project_sample_ids),
				preprocessed_rds_objects = flatten(preprocess.qc_rds_object),
				preprocessing_output_file_paths = flatten(preprocessing_output_file_paths),
				filter_cells_min_counts = filter_cells_min_counts,
				filter_genes_min_cells = filter_genes_min_cells,
				workflow_name = workflow_name,
				workflow_version = workflow_version,
				workflow_release = workflow_release,
				run_timestamp = get_workflow_metadata.timestamp,
				raw_data_path_prefix = cohort_raw_data_path_prefix,
				staging_data_buckets = cohort_staging_data_buckets,
				billing_project = get_workflow_metadata.billing_project,
				container_registry = container_registry,
				zones = zones
		}
	}

	output {
		# Sample-level outputs
		## Sample list
		Array[Array[Array[String]]] project_sample_ids = preprocess.project_sample_ids

		## Preprocess - sample-level
		Array[Array[File]] geomxngs_dcc_zip = preprocess.geomxngs_dcc_zip
		Array[Array[File]] geomxngs_output_tar_gz = preprocess.geomxngs_output_tar_gz
		## Preprocess - team-level
		Array[File] initial_rds_object = preprocess.initial_rds_object
		Array[File] sample_overview_sankey_html = preprocess.sample_overview_sankey_html
		Array[File] qc_rds_object = preprocess.qc_rds_object
		Array[File] segment_qc_summary_csv = preprocess.segment_qc_summary_csv
		Array[File] probe_qc_summary_csv = preprocess.probe_qc_summary_csv
		Array[File] gene_count_csv = preprocess.gene_count_csv

		# Project cohort analysis outputs
		## List of samples included in the cohort
		Array[File?] project_cohort_sample_list = project_cohort_analysis.cohort_sample_list

		# Merged RDS objects, filtered and normalized RDS objects, clustered adata objects, and plots
		Array[File?] project_merged_rds_object = project_cohort_analysis.merged_rds_object
		Array[File?] project_filtered_normalized_adata_object = project_cohort_analysis.filtered_normalized_adata_object
		Array[File?] project_umap_cluster_adata_object = project_cohort_analysis.umap_cluster_adata_object
		Array[File?] project_umap_cluster_plot_png = project_cohort_analysis.umap_cluster_plot_png

		# Spatial statistics outputs
		Array[File?] project_moran_adata_object = project_cohort_analysis.moran_adata_object
		Array[File?] project_moran_top_10_variable_genes_csv = project_cohort_analysis.moran_top_10_variable_genes_csv
		Array[File?] project_nhood_enrichment_adata_object = project_cohort_analysis.nhood_enrichment_adata_object
		Array[File?] project_nhood_enrichment_plot_png = project_cohort_analysis.nhood_enrichment_plot_png
		Array[File?] project_final_adata_object = project_cohort_analysis.final_adata_object
		Array[File?] project_co_occurrence_plot_png = project_cohort_analysis.co_occurrence_plot_png

		Array[Array[File]?] preprocess_manifests = project_cohort_analysis.preprocess_manifest_tsvs
		Array[Array[File]?] project_manifests = project_cohort_analysis.cohort_analysis_manifest_tsvs

		# Cross-team cohort analysis outputs
		## List of samples included in the cohort
		File? cohort_cohort_sample_list = cross_team_cohort_analysis.cohort_sample_list

		# Merged RDS objects, filtered and normalized RDS objects, clustered adata objects, and plots
		File? cohort_merged_rds_object = cross_team_cohort_analysis.merged_rds_object
		File? cohort_filtered_normalized_adata_object = cross_team_cohort_analysis.filtered_normalized_adata_object
		File? cohort_umap_cluster_adata_object = cross_team_cohort_analysis.umap_cluster_adata_object
		File? cohort_umap_cluster_plot_png = cross_team_cohort_analysis.umap_cluster_plot_png

		# Spatial statistics outputs
		File? cohort_moran_adata_object = cross_team_cohort_analysis.moran_adata_object
		File? cohort_moran_top_10_variable_genes_csv = cross_team_cohort_analysis.moran_top_10_variable_genes_csv
		File? cohort_nhood_enrichment_adata_object = cross_team_cohort_analysis.nhood_enrichment_adata_object
		File? cohort_nhood_enrichment_plot_png = cross_team_cohort_analysis.nhood_enrichment_plot_png
		File? cohort_final_adata_object = cross_team_cohort_analysis.final_adata_object
		File? cohort_co_occurrence_plot_png = cross_team_cohort_analysis.co_occurrence_plot_png

		Array[File]? cohort_manifests = cross_team_cohort_analysis.cohort_analysis_manifest_tsvs
	}

	meta {
		description: "Harmonized human and non-human postmortem-derived brain sequencing (PMDBS) spatial transcriptomics workflow for Nanostring GeoMx data"
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files during cross-team downstream analysis."}
		projects: {help: "The project ID, set of samples and their associated reads and metadata, output bucket locations, and whether or not to run project-level downstream analysis."}
		geomxngs_config_pkc: {help: "The GeoMx DSP configuration file to associate assay targets with GeoMx HybCode barcodes and Seq Code primers."}
		min_segment_reads: {help: "Minimum number of segment reads. [1000]"}
		min_percent_reads_trimmed: {help: "Minimum % of reads trimmed. [80]"}
		min_percent_reads_stitched: {help: "Minimum % of reads stitched. [80]"}
		min_percent_reads_aligned: {help: "Minimum % of reads aligned. [80]"}
		min_saturation: {help: "Minimum sequencing saturation. [50]"}
		min_neg_ctrl_count: {help: "Minimum negative control counts. [10]"}
		max_ntc_count: {help: "Maximum counts observed in NTC well. [1000]"}
		min_nuclei: {help: "Minimum # of nuclei estimated. [100]"}
		min_segment_area: {help: "Minimum segment area. [5000]"}
		
		filter_cells_min_counts: {help: "Minimum number of counts required for a cell to pass filtering. [5000]"}
		filter_genes_min_cells: {help: "Minimum number of cells expressed required for a gene to pass filtering. [10]"}
		run_cross_team_cohort_analysis: {help: "Whether to run downstream harmonization steps on all samples across projects. If set to false, only preprocessing steps (GeoMxNGSPipeline and generating the initial adata object(s)) will run for samples. [false]"}
		cohort_raw_data_bucket: {help: "Bucket to upload cross-team downstream intermediate files to."}
		cohort_staging_data_buckets: {help: "Set of buckets to stage cross-team downstream analysis outputs in."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones where compute will take place."}
	}
}
