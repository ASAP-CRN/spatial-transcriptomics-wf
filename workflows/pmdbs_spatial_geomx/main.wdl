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

		File config_ini
		File geomxngs_config_pkc

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
				config_ini = config_ini,
				geomxngs_config_pkc = geomxngs_config_pkc,
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
			preprocess.initial_adata_object,
			preprocess.qc_adata_object
		]) #!StringCoercion

		if (project.run_project_cohort_analysis) {
			call CohortAnalysis.cohort_analysis as project_cohort_analysis {
				input:
					cohort_id = project.team_id,
					project_sample_ids = preprocess.project_sample_ids,
					preprocessed_adata_objects = preprocess.qc_adata_object,
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
				preprocessed_adata_objects = flatten(preprocess.qc_adata_object),
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

		## Preprocess
		Array[Array[File]] geomxngs_dcc_zip = preprocess.geomxngs_dcc_zip
		Array[Array[File]] geomxngs_output_tar_gz = preprocess.geomxngs_output_tar_gz
		Array[Array[File]] initial_adata_object = preprocess.initial_adata_object
		Array[Array[File]] qc_adata_object = preprocess.qc_adata_object
		Array[Array[Float?]] qc_unassigned_ctrl_probes_percentage = preprocess.qc_unassigned_ctrl_probes_percentage

		# Project cohort analysis outputs
		## List of samples included in the cohort
		Array[File?] project_cohort_sample_list = project_cohort_analysis.cohort_sample_list

		# Merged adata objects, filtered and normalized adata objects, clustered adata objects, and plots
		Array[File?] project_merged_adata_object = project_cohort_analysis.merged_adata_object
		Array[File?] project_qc_plots_png = project_cohort_analysis.qc_plots_png
		Array[File?] project_filtered_normalized_adata_object = project_cohort_analysis.filtered_normalized_adata_object
		Array[File?] project_umap_cluster_adata_object = project_cohort_analysis.umap_cluster_adata_object
		Array[Array[File]?] project_umap_and_spatial_coord_plots_png = project_cohort_analysis.umap_and_spatial_coord_plots_png

		# Spatial statistics outputs
		Array[File?] project_nhood_enrichment_adata_object = project_cohort_analysis.nhood_enrichment_adata_object
		Array[File?] project_nhood_enrichment_plot_png = project_cohort_analysis.nhood_enrichment_plot_png
		Array[File?] project_co_occurrence_adata_object = project_cohort_analysis.co_occurrence_adata_object
		Array[File?] project_co_occurrence_plot_png = project_cohort_analysis.co_occurrence_plot_png
		Array[File?] project_final_adata_object = project_cohort_analysis.final_adata_object
		Array[File?] project_moran_top_10_variable_genes_csv = project_cohort_analysis.moran_top_10_variable_genes_csv

		Array[Array[File]?] preprocess_manifests = project_cohort_analysis.preprocess_manifest_tsvs
		Array[Array[File]?] project_manifests = project_cohort_analysis.cohort_analysis_manifest_tsvs

		# Cross-team cohort analysis outputs
		## List of samples included in the cohort
		File? cohort_cohort_sample_list = cross_team_cohort_analysis.cohort_sample_list

		# Merged adata objects, filtered and normalized adata objects, clustered adata objects, and plots
		File? cohort_merged_adata_object = cross_team_cohort_analysis.merged_adata_object
		File? cohort_qc_plots_png = cross_team_cohort_analysis.qc_plots_png
		File? cohort_filtered_normalized_adata_object = cross_team_cohort_analysis.filtered_normalized_adata_object
		File? cohort_umap_cluster_adata_object = cross_team_cohort_analysis.umap_cluster_adata_object
		Array[File]? cohort_umap_and_spatial_coord_plots_png = cross_team_cohort_analysis.umap_and_spatial_coord_plots_png

		# Spatial statistics outputs
		File? cohort_nhood_enrichment_adata_object = cross_team_cohort_analysis.nhood_enrichment_adata_object
		File? cohort_nhood_enrichment_plot_png = cross_team_cohort_analysis.nhood_enrichment_plot_png
		File? cohort_co_occurrence_adata_object = cross_team_cohort_analysis.co_occurrence_adata_object
		File? cohort_co_occurrence_plot_png = cross_team_cohort_analysis.co_occurrence_plot_png
		File? cohort_final_adata_object = cross_team_cohort_analysis.final_adata_object
		File? cohort_moran_top_10_variable_genes_csv = cross_team_cohort_analysis.moran_top_10_variable_genes_csv

		Array[File]? cohort_manifests = cross_team_cohort_analysis.cohort_analysis_manifest_tsvs
	}

	meta {
		description: "Harmonized human and non-human postmortem-derived brain sequencing (PMDBS) spatial transcriptomics workflow for Nanostring GeoMx data"
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files during cross-team downstream analysis."}
		projects: {help: "The project ID, set of samples and their associated reads and metadata, output bucket locations, and whether or not to run project-level downstream analysis."}
		config_ini: {help: "The configuration (.ini) file, containing pipeline processing parameters."}
		geomxngs_config_pkc: {help: "The GeoMx DSP configuration file to associate assay targets with GeoMx HybCode barcodes and Seq Code primers."}
		filter_cells_min_counts: {help: "Minimum number of counts required for a cell to pass filtering. [5000]"}
		filter_genes_min_cells: {help: "Minimum number of cells expressed required for a gene to pass filtering. [10]"}
		run_cross_team_cohort_analysis: {help: "Whether to run downstream harmonization steps on all samples across projects. If set to false, only preprocessing steps (GeoMxNGSPipeline and generating the initial adata object(s)) will run for samples. [false]"}
		cohort_raw_data_bucket: {help: "Bucket to upload cross-team downstream intermediate files to."}
		cohort_staging_data_buckets: {help: "Set of buckets to stage cross-team downstream analysis outputs in."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones where compute will take place."}
	}
}
