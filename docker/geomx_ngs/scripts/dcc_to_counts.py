import os
import re
import configparser
import json
import argparse
import pandas as pd
import numpy as np
from anndata import AnnData


# Modified from https://github.com/lilab-bcb/cumulus/blob/master/docker/geomxngs_dcc_to_count_matrix/1.0.0/src/create-counts.py

def main(args):
    ##########################
    ## PREP CONFIG INI FILE ##
    ##########################
    config = configparser.ConfigParser()
    config.optionxform = str # prevent conversion of keys to lowercase
    config.read(args.ini)

    targets = list(config["Targets"].keys())
    if len(targets) == 0:
        raise ValueError("No targets found")


    ###########################
    ## CONVERT DCC TO COUNTS ##
    ###########################
    target2index = {}
    for i in range(len(targets)):
        target2index[targets[i]] = i

    dcc_file_paths = [
        os.path.join(args.dcc_files_dir_input, file)
        for file in os.listdir(args.dcc_files_dir_input)
        if os.path.isfile(os.path.join(args.dcc_files_dir_input, file))
    ]

    # Maps targets to matrix indices and stores counts
    matrix = np.zeros(shape=(len(dcc_file_paths), len(targets)), dtype="int")
    for i in range(len(dcc_file_paths)):
        dcc_file_path = dcc_file_paths[i]
        with open(dcc_file_path, "rt") as dcc_in:
            for line in dcc_in:
                line = line.strip()
                if line == "<Code_Summary>":
                    break
            for line in dcc_in:
                line = line.strip()
                if line == "</Code_Summary>":
                    break
                tokens = line.split(",")
                assert len(tokens) == 2
                target = tokens[0]
                target_index = target2index[target]
                count = int(tokens[1])
                matrix[i, target_index] = count
                # e.g. RTS0052259,3


    ######################
    ## INCLUDE PKC FILE ##
    ######################
    pkcs = []
    with open(args.pkc, "rb") as pkc_file:
        pkcs.append(json.load(pkc_file))
    for pkc in pkcs:
        rts_ids = []
        genes = []
        probes = []
        genome = pkc["Name"]
        for target in pkc["Targets"]:
            display_name = target["DisplayName"]
            for probe in target["Probes"]:
                probes.append(probe["DisplayName"])
                rts_ids.append(probe["RTS_ID"])
                genes.append(display_name)

        pkc_df = pd.DataFrame(index=rts_ids, data={"gene": genes, "probe": probes})
        pkc_df["genome"] = genome
    var = pd.DataFrame(index=targets)
    obs = pd.DataFrame(index=[os.path.splitext(os.path.basename(path))[0] for path in dcc_file_paths])
    obs["plate_well"] = obs.index.str.split(r"[-_]").str[-2:].str.join("-")  # e.g. B-A02
    var = var.join(pkc_df)


    ###########################
    ## GENERATE ADATA OBJECT ##
    ###########################
    matrix = matrix.astype("int")
    adata = AnnData(X=matrix, obs=obs, var=var)
    for col in adata.obs.columns:
        if pd.api.types.is_object_dtype(adata.obs[col].dtype):
            adata.obs[col] = adata.obs[col].astype(str)

    # Add metadata
    adata.obs["team"] = args.team_id
    adata.obs["dataset"] = args.dataset_id
    adata.obs["sample"] = args.sample_id
    adata.obs["batch"] = args.batch
    adata.obs["batch_id"] = f"{args.team_id}_{args.dataset_id}_{args.batch}"
    
    # Save adata object
    adata.write_h5ad(filename=args.counts_output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert GeoMxNGSPipeline DCC files to count matrices"
    )
    parser.add_argument(
        "-t",
        "--team-id",
        type=str,
        required=True,
        help="Team ID"
    )
    parser.add_argument(
        "-d",
        "--dataset-id",
        type=str,
        required=True,
        help="Dataset ID"
    )
    parser.add_argument(
        "-s",
        "--sample-id",
        type=str,
        required=True,
        help="Sample ID"
    )
    parser.add_argument(
        "-b",
        "--batch",
        type=str,
        required=True,
        help="Batch from which the sample/dataset originated"
    )
    parser.add_argument(
        "-i",
        "--dcc-files-dir-input",
        type=str,
        required=True,
        help="Directory containing the DCC files for a sample"
    )
    parser.add_argument(
        "-c",
        "--ini",
        type=str,
        required=True,
        help="The configuration (.ini) file, containing pipeline processing parameters"
    )
    parser.add_argument(
        "-p",
        "--pkc",
        type=str,
        required=True,
        help="The GeoMx DSP configuration file to associate assay targets with GeoMx HybCode barcodes and Seq Code primers"
    )
    parser.add_argument(
        "-o",
        "--counts-output",
        type=str,
        required=True,
        help="Output file name for the counts matrix AnnData object"
    )

    args = parser.parse_args()

    main(args)
