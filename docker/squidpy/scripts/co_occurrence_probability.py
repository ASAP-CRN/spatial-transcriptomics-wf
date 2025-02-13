import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq


def main(args):
    ###############################
    ## CO-OCCURRENCE PROBABILITY ##
    ###############################
    adata = sc.read_h5ad(args.adata_input)

    sq.gr.co_occurrence(
        adata,
        cluster_key="leiden_res_0.40",
    )

    sq.pl.co_occurrence(
        adata,
        cluster_key="leiden_res_0.40",
    )
    plt.savefig(f"{args.cohort_id}.co_occurrence.png", dpi=300, bbox_inches="tight")

    # Save adata object
    adata.write_h5ad(filename=args.counts_output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute the co-occurrence probability: measures how often two cell types appear near each other in spatial transcriptomics data"
    )
    parser.add_argument(
        "-c",
        "--cohort-id",
        type=str,
        required=True,
        help="Cohort ID"
    )
    parser.add_argument(
        "-i",
        "--adata-input",
        type=str,
        required=True,
        help="Adata object to compute co-occurrence probability on"
    )
    parser.add_argument(
        "-o",
        "--adata-output",
        type=str,
        required=True,
        help="Output file name for the co-occurrence probability ratio AnnData object"
    )

    args = parser.parse_args()

    main(args)
