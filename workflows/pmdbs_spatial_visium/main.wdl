version 1.0

# Harmonized human and non-human PMDBS spatial transcriptomics workflow entrypoint for 10x Visium data

import "../../wf-common/wdl/structs.wdl"
import "../../wf-common/wdl/tasks/get_workflow_metadata.wdl" as GetWorkflowMetadata
import "preprocess/preprocess.wdl" as Preprocess
import "image_analysis/image_analysis.wdl" as ImageAnalysis
import "cohort_analysis/cohort_analysis.wdl" as CohortAnalysis

workflow pmdbs_spatial_visium_analysis {
	input {
		String cohort_id
		Array[Project] projects

		File spaceranger_reference_data
		File visium_probe_set_csv

		# Filter parameters
		Int filter_cells_min_counts = 5000
		Int filter_genes_min_cells = 10

		# Feature selection and clustering inputs and parameters
		String batch_key = "batch_id"
		String scvi_latent_key = "X_scvi"
		Int n_top_genes = 3000
		File cell_type_markers_list

		# Cohort analysis
		Boolean run_cross_team_cohort_analysis = false
		String cohort_raw_data_bucket
		Array[String] cohort_staging_data_buckets

		String container_registry
		String zones = "us-central1-c us-central1-f"
	}

	String workflow_execution_path = "workflow_execution"
	String workflow_name = "pmdbs_spatial_visium"
	String workflow_version = "v1.0.0"
	String workflow_release = "https://github.com/ASAP-CRN/pmdbs-spatial-transcriptomics-wf/releases/tag/pmdbs_spatial_visium_analysis-~{workflow_version}"

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
				spaceranger_reference_data = spaceranger_reference_data,
				visium_probe_set_csv = visium_probe_set_csv,
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
			preprocess.raw_counts,
			preprocess.filtered_counts,
			preprocess.molecule_info,
			preprocess.metrics_summary_csv,
			preprocess.spatial_outputs_tar_gz,
			flatten(preprocess.spatial_images),
			preprocess.scalefactors_json,
			preprocess.tissue_positions_csv,
			preprocess.spatial_enrichment_csv,
			preprocess.initial_adata_object,
			preprocess.qc_adata_object
		]) #!StringCoercion

		call ImageAnalysis.image_analysis {
			input:
				preprocessed_adata_objects = preprocess.qc_adata_object,
				raw_data_path_prefix = project_raw_data_path_prefix,
				billing_project = get_workflow_metadata.billing_project,
				container_registry = container_registry,
				zones = zones
		}

		if (project.run_project_cohort_analysis) {
			call CohortAnalysis.cohort_analysis as project_cohort_analysis {
				input:
					cohort_id = project.team_id,
					project_sample_ids = preprocess.project_sample_ids,
					preprocessed_adata_objects = image_analysis.image_features_adata_object,
					preprocessing_output_file_paths = preprocessing_output_file_paths,
					filter_cells_min_counts = filter_cells_min_counts,
					filter_genes_min_cells = filter_genes_min_cells,
					batch_key = batch_key,
					scvi_latent_key = scvi_latent_key,
					n_top_genes = n_top_genes,
					cell_type_markers_list = cell_type_markers_list,
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
				preprocessed_adata_objects = flatten(image_analysis.image_features_adata_object),
				preprocessing_output_file_paths = flatten(preprocessing_output_file_paths),
				filter_cells_min_counts = filter_cells_min_counts,
				filter_genes_min_cells = filter_genes_min_cells,
				batch_key = batch_key,
				scvi_latent_key = scvi_latent_key,
				n_top_genes = n_top_genes,
				cell_type_markers_list = cell_type_markers_list,
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
		Array[Array[File]] raw_counts = preprocess.raw_counts
		Array[Array[File]] filtered_counts = preprocess.filtered_counts
		Array[Array[File]] molecule_info = preprocess.molecule_info
		Array[Array[File]] metrics_summary_csv = preprocess.metrics_summary_csv
		Array[Array[File]] spatial_outputs_tar_gz = preprocess.spatial_outputs_tar_gz
		Array[Array[Array[File]]] spatial_images = preprocess.spatial_images
		Array[Array[File]] scalefactors_json = preprocess.scalefactors_json
		Array[Array[File]] tissue_positions_csv = preprocess.tissue_positions_csv
		Array[Array[File]] spatial_enrichment_csv = preprocess.spatial_enrichment_csv
		Array[Array[File]] initial_adata_object = preprocess.initial_adata_object
		Array[Array[File]] qc_adata_object = preprocess.qc_adata_object

		## Image analysis
		Array[Array[File]] image_features_adata_object = image_analysis.image_features_adata_object

		# Project cohort analysis outputs
		## List of samples included in the cohort
		Array[File?] project_cohort_sample_list = project_cohort_analysis.cohort_sample_list

		# Merged adata objects, filtered and normalized adata objects, clustered adata objects, and plots
		Array[File?] project_merged_adata_object = project_cohort_analysis.merged_adata_object
		Array[Array[File]?] project_qc_plots_png = project_cohort_analysis.qc_plots_png
		Array[File?] project_filtered_normalized_adata_object = project_cohort_analysis.filtered_normalized_adata_object
		Array[File?] project_feature_selection_adata_object = project_cohort_analysis.feature_selection_adata_object
		Array[File?] project_feature_dispersion_plot_png = project_cohort_analysis.feature_dispersion_plot_png

		# Clustering outputs
		Array[File?] project_integrated_adata_object = project_cohort_analysis.integrated_adata_object
		Array[File?] project_scvi_model_tar_gz = project_cohort_analysis.scvi_model_tar_gz
		Array[File?] project_umap_cluster_adata_object = project_cohort_analysis.umap_cluster_adata_object
		Array[File?] project_cell_annotated_adata_object = project_cohort_analysis.cell_annotated_adata_object
		Array[File?] project_cell_types_csv = project_cohort_analysis.cell_types_csv

		# Spatial plots outputs
		Array[File?] project_features_umap_plot_png = project_cohort_analysis.features_umap_plot_png
		Array[File?] project_groups_umap_plot_png = project_cohort_analysis.groups_umap_plot_png
		Array[File?] project_image_features_spatial_scatter_plot_png = project_cohort_analysis.image_features_spatial_scatter_plot_png

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

		# Merged adata objects, filtered and normalized adata objects, clustered adata objects, and plots
		File? cohort_merged_adata_object = cross_team_cohort_analysis.merged_adata_object
		Array[File]? cohort_qc_plots_png = cross_team_cohort_analysis.qc_plots_png
		File? cohort_filtered_normalized_adata_object = cross_team_cohort_analysis.filtered_normalized_adata_object
		File? cohort_feature_selection_adata_object = cross_team_cohort_analysis.feature_selection_adata_object
		File? cohort_feature_dispersion_plot_png = cross_team_cohort_analysis.feature_dispersion_plot_png

		# Clustering outputs
		File? cohort_integrated_adata_object = cross_team_cohort_analysis.integrated_adata_object
		File? cohort_scvi_model_tar_gz = cross_team_cohort_analysis.scvi_model_tar_gz
		File? cohort_umap_cluster_adata_object = cross_team_cohort_analysis.umap_cluster_adata_object
		File? cohort_cell_annotated_adata_object = cross_team_cohort_analysis.cell_annotated_adata_object
		File? cohort_cell_types_csv = cross_team_cohort_analysis.cell_types_csv

		# Spatial plots outputs
		File? cohort_features_umap_plot_png = cross_team_cohort_analysis.features_umap_plot_png
		File? cohort_groups_umap_plot_png = cross_team_cohort_analysis.groups_umap_plot_png
		File? cohort_image_features_spatial_scatter_plot_png = cross_team_cohort_analysis.image_features_spatial_scatter_plot_png

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
		description: "Harmonized human and non-human postmortem-derived brain sequencing (PMDBS) spatial transcriptomics workflow for 10x Visium data"
	}

	parameter_meta {
		cohort_id: {help: "Name of the cohort; used to name output files during cross-team downstream analysis."}
		projects: {help: "The project ID, set of samples and their associated reads and metadata, output bucket locations, and whether or not to run project-level downstream analysis."}
		spaceranger_reference_data: {help: "Space Ranger transcriptome reference data; see https://www.10xgenomics.com/support/software/space-ranger/downloads."}
		visium_probe_set_csv: {help: "Visium probe-based assays target genes in Space Ranger transcriptome; see https://www.10xgenomics.com/support/software/space-ranger/downloads."}
		filter_cells_min_counts: {help: "Minimum number of counts required for a cell to pass filtering. [5000]"}
		filter_genes_min_cells: {help: "Minimum number of cells expressed required for a gene to pass filtering. [10]"}
		batch_key: {help: "Key in AnnData object for batch information so that highly-variable genes are selected within each batch separately and merged. ['batch_id']"}
		scvi_latent_key: {help: "Latent key to save the scVI latent to. ['X_scvi']"}
		n_top_genes: {help: "Number of highly-variable genes to keep. [3000]"}
		cell_type_markers_list: {help: "CSV file containing a list of major cell type markers; used to annotate clusters."}
		run_cross_team_cohort_analysis: {help: "Whether to run downstream harmonization steps on all samples across projects. If set to false, only preprocessing steps (GeoMxNGSPipeline and generating the initial adata object(s)) will run for samples. [false]"}
		cohort_raw_data_bucket: {help: "Bucket to upload cross-team downstream intermediate files to."}
		cohort_staging_data_buckets: {help: "Set of buckets to stage cross-team downstream analysis outputs in."}
		container_registry: {help: "Container registry where workflow Docker images are hosted."}
		zones: {help: "Space-delimited set of GCP zones where compute will take place."}
	}
}
