import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc


# https://scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html

def main(args):
    #############################################
    ## GENERATE ADATA OBJECT WITH SPATIAL DATA ##
    #############################################
    adata = sc.read_visium(
        path=args.spaceranger_spatial_dir,
        library_id=args.sample_id, # TODO
    )

    # Add metadata
    adata.obs["team"] = args.team_id
    adata.obs["dataset"] = args.dataset_id
    adata.obs["sample"] = args.sample_id
    adata.obs["batch"] = args.batch
    adata.obs["batch_id"] = f"{args.team_id}_{args.dataset_id}_{args.batch}"
    
    # Save adata object
    adata.write_h5ad(filename=args.adata_output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert Space Ranger counts to adata objects"
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
        "-f",
        "--spaceranger-spatial-dir",
        type=str,
        required=True,
        help="Spaceranger spatial outputs directory/folder"
    )
    parser.add_argument(
        "-o",
        "--adata-output",
        type=str,
        required=True,
        help="Output file name for the AnnData object"
    )

    args = parser.parse_args()

    main(args)
