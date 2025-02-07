version 1.0

# Generate a QC'ed AnnData object by converting FASTQ files to DCC (digital count conversion) files to count matrices

import "../../../wf-common/wdl/structs.wdl"

workflow preprocess {
	input {
		String team_id
		String dataset_id
		Array[Sample] samples

		File config_ini
		File geomxngs_config_pkc

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
	String fastq_to_dcc_task_version = "1.0.0"
	String dcc_to_adata_task_version = "1.0.0"
	String qc_task_version = "1.0.0"

	Array[Array[String]] workflow_info = [[run_timestamp, workflow_name, workflow_version, workflow_release]]

	String workflow_raw_data_path_prefix = "~{raw_data_path_prefix}/~{sub_workflow_name}"
	String dcc_raw_data_path = "~{workflow_raw_data_path_prefix}/fastq_to_dcc/~{fastq_to_dcc_task_version}"
	String adata_raw_data_path = "~{workflow_raw_data_path_prefix}/dcc_to_adata/~{dcc_to_adata_task_version}"
	String qc_raw_data_path = "~{workflow_raw_data_path_prefix}/qc/~{qc_task_version}"

	scatter (sample_object in samples) {
		String fastq_to_dcc_output = "~{dcc_raw_data_path}/~{sample_object.sample_id}.geomxngs_out_dir.tar.gz"
		String dcc_to_adata_output = "~{adata_raw_data_path}/~{sample_object.sample_id}.initial_adata_object.h5ad"
		String qc_output = "~{qc_raw_data_path}/~{sample_object.sample_id}.qc.h5ad"
	}

	# For each sample, outputs an array of true/false: [fastq_to_dcc_complete, dcc_to_adata_complete, qc_complete]
	call check_output_files_exist {
		input:
			fastq_to_dcc_output_files = fastq_to_dcc_output,
			dcc_to_adata_output_files = dcc_to_adata_output,
			qc_output_files = qc_output,
			billing_project = billing_project,
			zones = zones
	}

	scatter (sample_index in range(length(samples))) {
		Sample sample = samples[sample_index]

		Array[String] project_sample_id = [team_id, sample.sample_id]

		String fastq_to_dcc_complete = check_output_files_exist.sample_preprocessing_complete[sample_index][0]
		String dcc_to_adata_complete = check_output_files_exist.sample_preprocessing_complete[sample_index][1]
		String qc_complete = check_output_files_exist.sample_preprocessing_complete[sample_index][2]

		String fastq_to_dcc_geomxngs_dcc_zip = "~{dcc_raw_data_path}/~{sample.sample_id}.DCC.zip"
		String fastq_to_dcc_geomxngs_output_tar_gz = "~{dcc_raw_data_path}/~{sample.sample_id}.geomxngs_out_dir.tar.gz"

		if (fastq_to_dcc_complete == "false") {
			call fastq_to_dcc {
				input:
					sample_id = sample.sample_id,
					fastq_R1s = sample.fastq_R1s,
					fastq_R2s = sample.fastq_R2s,
					config_ini = config_ini,
					raw_data_path = dcc_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		File geomxngs_dcc_zip_output = select_first([fastq_to_dcc.geomxngs_dcc_zip, fastq_to_dcc_geomxngs_dcc_zip]) #!FileCoercion
		File geomxngs_output_tar_gz_output = select_first([fastq_to_dcc.geomxngs_output_tar_gz, fastq_to_dcc_geomxngs_output_tar_gz]) #!FileCoercion

		String dcc_to_adata_object = "~{adata_raw_data_path}/~{sample.sample_id}.initial_adata_object.h5ad"

		if (dcc_to_adata_complete == "false") {
			call dcc_to_adata {
				input:
					team_id = team_id,
					dataset_id = dataset_id,
					sample_id = sample.sample_id,
					batch = select_first([sample.batch]),
					geomxngs_dcc_zip = geomxngs_dcc_zip_output,
					config_ini = config_ini,
					geomxngs_config_pkc = geomxngs_config_pkc,
					raw_data_path = adata_raw_data_path,
					workflow_info = workflow_info,
					billing_project = billing_project,
					container_registry = container_registry,
					zones = zones
			}
		}

		File initial_adata_object_output = select_first([dcc_to_adata.initial_adata_object, dcc_to_adata_object]) #!FileCoercion

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

		# GeoMxNGSPipeline outputs including converted DCC files
		Array[File] geomxngs_dcc_zip = geomxngs_dcc_zip_output #!FileCoercion
		Array[File] geomxngs_output_tar_gz = geomxngs_output_tar_gz_output #!FileCoercion

		# Initial adata object
		Array[File] initial_adata_object = initial_adata_object_output #!FileCoercion

		# QC adata object
		Array[File] qc_adata_object = qc_adata_object_output #!FileCoercion
	}
}

task check_output_files_exist {
	input {
		Array[String] fastq_to_dcc_output_files
		Array[String] dcc_to_adata_output_files
		Array[String] qc_output_files

		String billing_project
		String zones
	}

	command <<<
		set -euo pipefail

		while read -r output_files || [[ -n "${output_files}" ]]; do
			dcc_file=$(echo "${output_files}" | cut -f 1)
			adata_file=$(echo "${output_files}" | cut -f 2)
			qc_adata_file=$(echo "${output_files}" | cut -f 3)

			if gsutil -u ~{billing_project} ls "${dcc_file}"; then
				if gsutil -u ~{billing_project} ls "${adata_file}"; then
					if gsutil -u ~{billing_project} ls "${qc_adata_file}"; then
						# If we find all outputs, don't rerun anything
						echo -e "true\ttrue\ttrue" >> sample_preprocessing_complete.tsv
					else
						# If we find fastq_to_dcc and dcc_to_adata outputs, then run (or rerun) qc
						echo -e "true\ttrue\tfalse" >> sample_preprocessing_complete.tsv
					fi
				else
					# If we only find fastq_to_dcc, then run (or rerun) dcc_to_adata and qc
					echo -e "true\tfalse\tfalse" >> sample_preprocessing_complete.tsv
				fi
			else
				# If we don't find fast_to_dcc output, we must need to run (or rerun) preprocessing
				echo -e "false\tfalse\tfalse" >> sample_preprocessing_complete.tsv
			fi
		done < <(paste ~{write_lines(fastq_to_dcc_output_files)} ~{write_lines(dcc_to_adata_output_files)} ~{write_lines(qc_output_files)})
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

task fastq_to_dcc {
	input {
		String sample_id

		Array[File] fastq_R1s
		Array[File] fastq_R2s

		File config_ini

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 8
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil(size(flatten([fastq_R1s, fastq_R2s]), "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		# Ensure fastqs are in the same directory
		mkdir fastqs
		while read -r fastq || [[ -n "${fastq}" ]]; do
			ln -s "${fastq}" "fastqs/$(basename "$fastq")"
		done < <(cat \
			~{write_lines(fastq_R1s)} \
			~{write_lines(fastq_R2s)})

		geomxngspipeline \
			--ini=~{config_ini} \
			--in="$(pwd)/fastqs" \
			--out="~{sample_id}_geomxngs_out_dir" \
			--check-illumina-naming=false

		# DCC zip file is automatically named following format: DCC-<YYYYMMDD>.zip
		dcc_file_to_rename=$(find ./~{sample_id}_geomxngs_out_dir -type f -name 'DCC-[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9].zip')
		cp "${dcc_file_to_rename}" ./~{sample_id}.DCC.zip

		tar -czvf "~{sample_id}.geomxngs_out_dir.tar.gz" "~{sample_id}_geomxngs_out_dir"

		upload_outputs \
			-b ~{billing_project} \
			-d ~{raw_data_path} \
			-i ~{write_tsv(workflow_info)} \
			-o "~{sample_id}.DCC.zip" \
			-o "~{sample_id}.geomxngs_out_dir.tar.gz"
	>>>

	output {
		String geomxngs_dcc_zip = "~{raw_data_path}/~{sample_id}.DCC.zip"
		String geomxngs_output_tar_gz = "~{raw_data_path}/~{sample_id}.geomxngs_out_dir.tar.gz"
	}

	runtime {
		docker: "~{container_registry}/geomx_ngs:3.1.1.6"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 30
		zones: zones
	}
}

task dcc_to_adata {
	input {
		String team_id
		String dataset_id
		String sample_id
		String batch

		File geomxngs_dcc_zip

		File config_ini
		File geomxngs_config_pkc

		String raw_data_path
		Array[Array[String]] workflow_info
		String billing_project
		String container_registry
		String zones
	}

	Int threads = 4
	Int mem_gb = ceil(threads * 2)
	Int disk_size = ceil(size([geomxngs_dcc_zip, config_ini, geomxngs_config_pkc], "GB") * 2 + 30)

	command <<<
		set -euo pipefail

		unzip -d ./dcc_files_dir ~{geomxngs_dcc_zip}

		python3 /opt/scripts/dcc_to_adata.py \
			--team ~{team_id} \
			--dataset ~{dataset_id} \
			--sample-id ~{sample_id} \
			--batch ~{batch} \
			--dcc-files-dir-input ./dcc_files_dir \
			--ini ~{config_ini} \
			--pkc ~{geomxngs_config_pkc} \
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
		docker: "~{container_registry}/geomx_ngs:3.1.1.6"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
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

		python3 /opt/scripts/geomx_qc.py \
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
		bootDiskSizeGb: 30
		zones: zones
	}
}
