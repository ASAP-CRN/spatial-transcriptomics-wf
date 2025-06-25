version 1.0

struct Sample {
	String sample_id
	String? batch

	Array[File]+ fastq_R1s
	Array[File]+ fastq_R2s
	Array[File] fastq_I1s
	Array[File] fastq_I2s
}

struct Slide {
	String slide_id
	File geomx_lab_annotation_xlsx

	Array[Sample] samples
}

struct Project {
	String team_id
	String dataset_id
	String dataset_doi_url
	Array[Slide] slides

	File project_sample_metadata_csv
	File geomx_config_ini

	Boolean run_project_cohort_analysis

	String raw_data_bucket
	Array[String] staging_data_buckets
}
