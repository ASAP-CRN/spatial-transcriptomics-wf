version 1.0

# Generate a AnnData object by converting FASTQ files to counts to adata objects

import "../../../wf-common/wdl/structs.wdl"

workflow preprocess {
	input {
		String team_id
		String dataset_id
		Array[Sample] samples

		File spaceranger_reference_data

		String workflow_name
		String workflow_version
		String workflow_release
		String run_timestamp
		String raw_data_path_prefix
		String billing_project
		String container_registry
		String zones
	}

	# Task and subworkflow versions
	String sub_workflow_name = "preprocess"
	String spaceranger_count_task_version = "1.0.0"
	String counts_to_adata_task_version = "1.0.0"
	String qc_task_version = "1.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version, workflow_release]]

	String workflow_raw_data_path_prefix = "~{raw_data_path_prefix}/~{sub_workflow_name}"
	String spaceranger_raw_data_path = "~{workflow_raw_data_path_prefix}/spaceranger_count/~{spaceranger_count_task_version}"
	String adata_raw_data_path = "~{workflow_raw_data_path_prefix}/counts_to_adata/~{counts_to_adata_task_version}"
	String qc_raw_data_path = "~{workflow_raw_data_path_prefix}/qc/~{qc_task_version}"

	scatter (sample_object in samples) {
		String spaceranger_count_output = "~{spaceranger_raw_data_path}/~{sample_object.sample_id}.raw_feature_bc_matrix.h5"
		String counts_to_adata_output = "~{adata_raw_data_path}/~{sample_object.sample_id}.initial_adata_object.h5ad"
		String qc_output = "~{qc_raw_data_path}/~{sample_object.sample_id}.qc.h5ad"
	}

	# For each sample, outputs an array of true/false: [spaceranger_count_complete, counts_to_adata_complete, qc_complete]
	call check_output_files_exist {
		input:
			spaceranger_count_output_files = spaceranger_count_output,
			counts_to_adata_output_files = counts_to_adata_output,
			qc_output_files = qc_output,
			billing_project = billing_project,
			zones = zones
	}

	scatter (sample_index in range(length(samples))) {
		Sample sample = samples[sample_index]

		Array[String] project_sample_id = [team_id, sample.sample_id]

		String spaceranger_count_complete = check_output_files_exist.sample_preprocessing_complete[sample_index][0]
		String counts_to_adata_complete = check_output_files_exist.sample_preprocessing_complete[sample_index][1]
		String qc_complete = check_output_files_exist.sample_preprocessing_complete[sample_index][2]

		String spaceranger_raw_counts = "~{spaceranger_raw_data_path}/~{sample.sample_id}.raw_feature_bc_matrix.h5"
		String spaceranger_filtered_counts = "~{spaceranger_raw_data_path}/~{sample.sample_id}.filtered_feature_bc_matrix.h5"
		String spaceranger_molecule_info = "~{spaceranger_raw_data_path}/~{sample.sample_id}.molecule_info.h5"
		String spaceranger_metrics_summary_csv = "~{spaceranger_raw_data_path}/~{sample.sample_id}.metrics_summary.csv"
		String spaceranger_spatial_outputs_tar_gz = "~{spaceranger_raw_data_path}/~{sample.sample_id}.spaceranger_spatial_outputs.tar.gz"
		Array[String] spaceranger_spatial_images = [
			"~{spaceranger_raw_data_path}/~{sample.sample_id}.aligned_fiducials.jpg",
			"~{spaceranger_raw_data_path}/~{sample.sample_id}.detected_tissue_image.jpg",
			"~{spaceranger_raw_data_path}/~{sample.sample_id}.tissue_hires_image.png",
			"~{spaceranger_raw_data_path}/~{sample.sample_id}.tissue_lowres_image.png"
		]
		String spaceranger_scalefactors_json = "~{spaceranger_raw_data_path}/~{sample.sample_id}.scalefactors_json.json"
		String spaceranger_tissue_positions_csv = "~{spaceranger_raw_data_path}/~{sample.sample_id}.tissue_positions.csv"
		String spaceranger_spatial_enrichment_csv = "~{spaceranger_raw_data_path}/~{sample.sample_id}.spatial_enrichment.csv"

		if (spaceranger_count_complete == "false") {
			call spaceranger_count {
				input:
					sample_id = sample.sample_id,
					fastq_R1s = sample.fastq_R1s,
					fastq_R2s = sample.fastq_R2s,
					fastq_I1s = sample.fastq_I1s,
					fastq_I2s = sample.fastq_I2s,
					visium_brightfield_image = sample.visium_brightfield_image,
					spaceranger_reference_data = spaceranger_reference_data,
					raw_data_path = spaceranger_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		File raw_counts_output = select_first([spaceranger_count.raw_counts, spaceranger_raw_counts]) #!FileCoercion
		File filtered_counts_output = select_first([spaceranger_count.filtered_counts, spaceranger_filtered_counts]) #!FileCoercion
		File molecule_info_output = select_first([spaceranger_count.molecule_info, spaceranger_molecule_info]) #!FileCoercion
		File metrics_summary_csv_output = select_first([spaceranger_count.metrics_summary_csv, spaceranger_metrics_summary_csv]) #!FileCoercion
		File spatial_outputs_tar_gz_output = select_first([spaceranger_count.spatial_outputs_tar_gz, spaceranger_spatial_outputs_tar_gz]) #!FileCoercion
		Array[File] spatial_images_output = select_first([spaceranger_count.spatial_images, spaceranger_spatial_images]) #!FileCoercion
		File scalefactors_json_output = select_first([spaceranger_count.scalefactors_json, spaceranger_scalefactors_json]) #!FileCoercion
		File tissue_positions_csv_output = select_first([spaceranger_count.tissue_positions_csv, spaceranger_tissue_positions_csv]) #!FileCoercion
		File spatial_enrichment_csv_output = select_first([spaceranger_count.spatial_enrichment_csv, spaceranger_spatial_enrichment_csv]) #!FileCoercion

		String counts_to_adata_object = "~{adata_raw_data_path}/~{sample.sample_id}.initial_adata_object.h5ad"

		if (counts_to_adata_complete == "false") {
			call counts_to_adata {
				input:
					team_id = team_id,
					dataset_id = dataset_id,
					sample_id = sample.sample_id,
					batch = select_first([sample.batch]),
					spaceranger_counts = filtered_counts_output,
					spaceranger_spatial_tar_gz = spatial_outputs_tar_gz_output,
					raw_data_path = adata_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		File initial_adata_object_output = select_first([counts_to_adata.initial_adata_object, counts_to_adata_object]) #!FileCoercion

		String qc_metrics_adata_object = "~{qc_raw_data_path}/~{sample.sample_id}.qc.h5ad"

		if (qc_complete == "false") {
			call qc {
				input:
					sample_id = sample.sample_id,
					initial_adata_object = initial_adata_object_output,
					raw_data_path = qc_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		File qc_adata_object_output = select_first([qc.qc_adata_object, qc_metrics_adata_object]) #!FileCoercion
	}

	output {
		# Sample list
		Array[Array[String]] project_sample_ids = project_sample_id

		# Spaceranger count outputs
		Array[File] raw_counts = raw_counts_output #!FileCoercion
		Array[File] filtered_counts = filtered_counts_output #!FileCoercion
		Array[File] molecule_info = molecule_info_output #!FileCoercion
		Array[File] metrics_summary_csv = metrics_summary_csv_output #!FileCoercion
		Array[File] spatial_outputs_tar_gz = spatial_outputs_tar_gz_output #!FileCoercion
		Array[Array[File]] spatial_images = spatial_images_output #!FileCoercion
		Array[File] scalefactors_json = scalefactors_json_output #!FileCoercion
		Array[File] tissue_positions_csv = tissue_positions_csv_output #!FileCoercion
		Array[File] spatial_enrichment_csv = spatial_enrichment_csv_output #!FileCoercion

		# Initial adata object
		Array[File] initial_adata_object = initial_adata_object_output #!FileCoercion

		# QC adata object
		Array[File] qc_adata_object = qc_adata_object_output #!FileCoercion
	}
}

task check_output_files_exist {
	input {
		Array[String] spaceranger_count_output_files
		Array[String] counts_to_adata_output_files
		Array[String] qc_output_files

		String billing_project
		String zones
	}

	command <<<
		set -euo pipefail

		while read -r output_files || [[ -n "${output_files}" ]]; do
			spaceranger_counts_file=$(echo "${output_files}" | cut -f 1)
			adata_file=$(echo "${output_files}" | cut -f 2)
			qc_adata_file=$(echo "${output_files}" | cut -f 3)

			if gsutil -u ~{billing_project} ls "${spaceranger_counts_file}"; then
				if gsutil -u ~{billing_project} ls "${adata_file}"; then
					if gsutil -u ~{billing_project} ls "${qc_adata_file}"; then
						# If we find all outputs, don't rerun anything
						echo -e "true\ttrue\ttrue" >> sample_preprocessing_complete.tsv
					else
						# If we find spaceranger_counts and counts_to_adata outputs, then run (or rerun) qc
						echo -e "true\ttrue\tfalse" >> sample_preprocessing_complete.tsv
					fi
				else
					# If we only find spaceranger_counts, then run (or rerun) counts_to_adata and qc
					echo -e "true\tfalse\tfalse" >> sample_preprocessing_complete.tsv
				fi
			else
				# If we don't find spaceranger_counts output, we must need to run (or rerun) preprocessing
				echo -e "false\tfalse\tfalse" >> sample_preprocessing_complete.tsv
			fi
		done < <(paste ~{write_lines(spaceranger_count_output_files)} ~{write_lines(counts_to_adata_output_files)} ~{write_lines(qc_output_files)})
	>>>

	output {
		Array[Array[String]] sample_preprocessing_complete = read_tsv("sample_preprocessing_complete.tsv")
	}

	runtime {
		docker: "gcr.io/google.com/cloudsdktool/google-cloud-cli:444.0.0-slim"
		cpu: 2
		memory: "4 GB"
		disks: "local-disk 20 HDD"
		preemptible: 3
		zones: zones
	}
}

task spaceranger_count {
	input {
		String sample_id

		Array[File] fastq_R1s
		Array[File] fastq_R2s
		Array[File] fastq_I1s
		Array[File] fastq_I2s
		File visium_brightfield_image

		File spaceranger_reference_data

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 16
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil(size(flatten([fastq_R1s, fastq_R2s]), "GB") + size([visium_brightfield_image, spaceranger_reference_data], "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		# Unpack refdata
		mkdir spaceranger_refdata
		tar \
			-zxvf ~{spaceranger_reference_data} \
			-C spaceranger_refdata \
			--strip-components 1

		# Ensure fastqs are in the same directory
		mkdir fastqs
		while read -r fastq || [[ -n "${fastq}" ]]; do
			if [[ -n "${fastq}" ]]; then
				validated_fastq_name=$(fix_fastq_names --fastq "${fastq}" --sample-id "~{sample_id}")
				if [[ -e "fastqs/${validated_fastq_name}" ]]; then
					echo "[ERROR] Something's gone wrong with fastq renaming; trying to create fastq [${validated_fastq_name}] but it already exists. Exiting."
					exit 1
				else
					ln -s "${fastq}" "fastqs/${validated_fastq_name}"
				fi
			fi
		done < <(cat \
			~{write_lines(fastq_R1s)} \
			~{write_lines(fastq_R2s)} \
			~{write_lines(fastq_I1s)} \
			~{write_lines(fastq_I2s)})

		spaceranger --version

		/usr/bin/time \
		spaceranger count \
			--id=~{sample_id} \
			--transcriptome="$(pwd)/spaceranger_refdata" \
			--fastqs="$(pwd)/fastqs" \
			--image=~{visium_brightfield_image} \ # TODO - need to rename image file as well? One image per sample?
			--slide=V19J01-123 \ # TODO - Slide ID should be in metadata
			--area=A1 \ # TODO - Capture area should be in metadata
			--localcores ~{threads} \
			--localmem ~{mem_gb - 4}

		# Rename outputs to include sample ID
		mv ~{sample_id}/outs/raw_feature_bc_matrix.h5 ~{sample_id}.raw_feature_bc_matrix.h5
		mv ~{sample_id}/outs/filtered_feature_bc_matrix.h5 ~{sample_id}.filtered_feature_bc_matrix.h5
		mv ~{sample_id}/outs/molecule_info.h5 ~{sample_id}.molecule_info.h5
		mv ~{sample_id}/outs/metrics_summary.csv ~{sample_id}.metrics_summary.csv

		cp -r ~{sample_id}/outs/spatial spatial
		tar -czvf "~{sample_id}.spaceranger_spatial_outputs.tar.gz" spatial

		mv ~{sample_id}/outs/spatial/aligned_fiducials.jpg ~{sample_id}.aligned_fiducials.jpg
		mv ~{sample_id}/outs/spatial/detected_tissue_image.jpg ~{sample_id}.detected_tissue_image.jpg
		mv ~{sample_id}/outs/spatial/scalefactors_json.json ~{sample_id}.scalefactors_json.json
		mv ~{sample_id}/outs/spatial/tissue_hires_image.png ~{sample_id}.tissue_hires_image.png
		mv ~{sample_id}/outs/spatial/tissue_lowres_image.png ~{sample_id}.tissue_lowres_image.png
		mv ~{sample_id}/outs/spatial/tissue_positions.csv ~{sample_id}.tissue_positions.csv
		mv ~{sample_id}/outs/spatial/spatial_enrichment.csv ~{sample_id}.spatial_enrichment.csv

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{sample_id}.raw_feature_bc_matrix.h5" \
			-o "~{sample_id}.filtered_feature_bc_matrix.h5" \
			-o "~{sample_id}.molecule_info.h5" \
			-o "~{sample_id}.metrics_summary.csv" \
			-o "~{sample_id}.spaceranger_spatial_outputs.tar.gz" \
			-o "~{sample_id}.aligned_fiducials.jpg" \
			-o "~{sample_id}.detected_tissue_image.jpg" \
			-o "~{sample_id}.scalefactors_json.json" \
			-o "~{sample_id}.tissue_hires_image.png" \
			-o "~{sample_id}.tissue_lowres_image.png" \
			-o "~{sample_id}.tissue_positions.csv" \
			-o "~{sample_id}.spatial_enrichment.csv"
	>>>

	output {
		String raw_counts = "~{raw_data_path}/~{sample_id}.raw_feature_bc_matrix.h5"
		String filtered_counts = "~{raw_data_path}/~{sample_id}.filtered_feature_bc_matrix.h5"
		String molecule_info = "~{raw_data_path}/~{sample_id}.molecule_info.h5"
		String metrics_summary_csv = "~{raw_data_path}/~{sample_id}.metrics_summary.csv"
		String spatial_outputs_tar_gz = "~{raw_data_path}/~{sample_id}.spaceranger_spatial_outputs.tar.gz"
		Array[String] spatial_images = [
			"~{raw_data_path}/~{sample_id}.aligned_fiducials.jpg",
			"~{raw_data_path}/~{sample_id}.detected_tissue_image.jpg",
			"~{raw_data_path}/~{sample_id}.tissue_hires_image.png",
			"~{raw_data_path}/~{sample_id}.tissue_lowres_image.png"
		]
		String scalefactors_json = "~{raw_data_path}/~{sample_id}.scalefactors_json.json"
		String tissue_positions_csv = "~{raw_data_path}/~{sample_id}.tissue_positions.csv"
		String spatial_enrichment_csv = "~{raw_data_path}/~{sample_id}.spatial_enrichment.csv"
	}

	runtime {
		docker: "~{container_registry}/spaceranger:3.1.2"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		zones: zones
	}
}

task counts_to_adata {
	input {
		String team_id
		String dataset_id
		String sample_id
		String batch

		File spaceranger_counts
		File spaceranger_spatial_tar_gz

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int disk_size = ceil(size([spaceranger_counts, spaceranger_spatial_tar_gz], "GB") * 2 + 20)

	command <<<
		set -euo pipefail

		tar -xzvf ~{spaceranger_spatial_tar_gz}

		python3 /opt/scripts/counts_to_adata.py \
			--team ~{team_id} \
			--dataset ~{dataset_id} \
			--sample-id ~{sample_id} \
			--batch ~{batch} \
			--spaceranger-counts-input ~{spaceranger_counts} \
			--spaceranger-spatial-dir spatial \
			--adata-output ~{sample_id}.initial_adata_object.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{sample_id}.initial_adata_object.h5ad"
	>>>

	output {
		String initial_adata_object = "~{raw_data_path}/~{sample_id}.initial_adata_object.h5ad"
	}

	runtime {
		docker: "~{container_registry}/squidpy:1.6.2_1"
		cpu: 2
		memory: "16 GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 40
		zones: zones
	}
}

task qc {
	input {
		String sample_id

		File initial_adata_object

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 4
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil(size(initial_adata_object, "GB") * 2 + 30)

	command <<<
		set -euo pipefail

		python3 /opt/scripts/visium_qc.py \
			--adata-input ~{initial_adata_object} \
			--qc-adata-output ~{sample_id}.qc.h5ad

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{sample_id}.qc.h5ad"
	>>>

	output {
		String qc_adata_object = "~{raw_data_path}/~{sample_id}.qc.h5ad"
	}

	runtime {
		docker: "~{container_registry}/squidpy:1.6.2_1"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		zones: zones
	}
}
