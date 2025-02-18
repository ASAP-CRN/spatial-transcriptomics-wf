import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq


def main(args):
    ###################################
    ## NEIGHBORS ENRICHMENT ANALYSIS ##
    ###################################
    adata = sc.read_h5ad(args.adata_input)

    obs_columns = list(adata.obs.columns)
    if "cell_type" in obs_columns:
        cluster_key="cell_type"
    else:
        raise ValueError("Cell type is missing in adata.obs")

    plt.rcParams["figure.figsize"] = (8, 8)
    sq.gr.nhood_enrichment(
        adata,
        library_key="sample",
        cluster_key=cluster_key,
    )
    sq.pl.nhood_enrichment(
        adata,
        cluster_key=cluster_key,
        title="Neighborhood enrichment adata",
    )
    plt.savefig(f"{args.cohort_id}.nhood_enrichment.png", dpi=300, bbox_inches="tight")

    # Save adata object
    adata.write_h5ad(filename=args.adata_output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Perform neighbors enrichment analysis by using cell type annotations"
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
        help="Adata object to run neighbors enrichment analysis on"
    )
    parser.add_argument(
        "-o",
        "--adata-output",
        type=str,
        required=True,
        help="Output file name for the neighborhood enrichment AnnData object"
    )

    args = parser.parse_args()

    main(args)
