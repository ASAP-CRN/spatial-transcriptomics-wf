version 1.0

# Calculate image features for all observations in a AnnData object

import "../../../wf-common/wdl/structs.wdl"

workflow image_analysis {
	input {
		Array[File] preprocessed_adata_objects

		String raw_data_path_prefix
		String billing_project
		String container_registry
		String zones
	}

	# Task and subworkflow versions
	String sub_workflow_name = "image_analysis"
	String image_features_task_version = "1.0.0"

	String workflow_raw_data_path_prefix = "~{raw_data_path_prefix}/~{sub_workflow_name}"
	String image_features_raw_data_path = "~{workflow_raw_data_path_prefix}/image_features/~{image_features_task_version}"

	scatter (preprocessed_adata_object in preprocessed_adata_objects) {
        String sample_id = basename(preprocessed_adata_object, ".qc.h5ad")
        String image_features_output = "~{image_features_raw_data_path}/~{sample_id}.image_features.h5ad"
    }

	Array[String] sample_ids = sample_id
	Array[Int] sample_indices = range(length(preprocessed_adata_objects))

	# For each sample, outputs an array of true/false: [image_features_complete]
	call check_output_files_exist {
		input:
			image_features_output_files = image_features_output,
			billing_project = billing_project,
			zones = zones
	}

	scatter (sample_index in sample_indices) {
		String image_features_complete = check_output_files_exist.sample_preprocessing_complete[sample_index][0]

		String squidpy_image_features_adata_object = "~{image_features_raw_data_path}/~{sample_ids[sample_index]}.image_features.h5ad"

		if (image_features_complete == "false") {
			call image_features {
				input:
					sample_id = sample_ids[sample_index],
					preprocessed_adata_object = preprocessed_adata_objects[sample_index],
					container_registry = container_registry,
					zones = zones
			}
		}

		File image_features_adata_object_output = select_first([image_features.image_features_adata_object, squidpy_image_features_adata_object]) #!FileCoercion
	}

	output {
		# Image features adata object
		Array[File] image_features_adata_object = image_features_adata_object_output #!FileCoercion
	}
}

task check_output_files_exist {
	input {
		Array[String] image_features_output_files

		String billing_project
		String zones
	}

	command <<<
		set -euo pipefail

		while read -r output_files || [[ -n "${output_files}" ]]; do
			image_features_file=$(echo "${output_files}" | cut -f 1)

			if gsutil -u ~{billing_project} ls "${image_features_file}"; then
				# If we find outputs, don't rerun anything
				echo -e "true" >> sample_preprocessing_complete.tsv
			else
				# If we don't find image_features_file output, we must need to run (or rerun) image_analysis
				echo -e "false" >> sample_preprocessing_complete.tsv
			fi
		done < <(paste ~{write_lines(image_features_output_files)})
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

task image_features {
	input {
		String sample_id
		File preprocessed_adata_object

		String container_registry
		String zones
	}

	Int threads = 4
	Int mem_gb = ceil(size(preprocessed_adata_object, "GB") * 2 + 20)
	Int disk_size = ceil(size(preprocessed_adata_object, "GB") * 2 + 50)

	command <<<
		set -euo pipefail

		python3 /opt/scripts/image_features.py \
			--sample-id ~{sample_id} \
			--adata-input ~{preprocessed_adata_object} \
			--n-jobs ~{threads} \
			--adata-output ~{sample_id}.image_features.h5ad
	>>>

	output {
		File image_features_adata_object = "~{sample_id}.image_features.h5ad"
	}

	runtime {
		docker: "~{container_registry}/spatial_py:1.0.0"
		cpu: threads
		memory: "~{mem_gb} GB"
		disks: "local-disk ~{disk_size} HDD"
		preemptible: 3
		bootDiskSizeGb: 30
		zones: zones
	}
}
