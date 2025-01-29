import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq


def main(args):
    ######################
    ## CLUSTER AND PLOT ##
    ######################
    adata = sc.read_h5ad(args.adata_input)

    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata)

    plt.rcParams["figure.figsize"] = (4, 4)
    sc.pl.umap(
        adata,
        color=[
            "total_counts",
            "n_genes_by_counts",
            "leiden",
        ],
        wspace=0.4,
        save=f"{args.cohort_id}.umap.png",
    )

    plt.rcParams["figure.figsize"] = (8, 8)
    sc.pl.spatial(
        adata,
        img_key="hires", # img_key: key where the img is stored in the adata.uns element
        color=[
            "total_counts",
            "n_genes_by_counts",
        ],
        save=f"{args.cohort_id}.spatial_coord_by_counts.png",
    )

    plt.rcParams["figure.figsize"] = (8, 8)
    sc.pl.spatial(
        adata,
        img_key="hires",
        color="clusters",
        size=1.5,
        save=f"{args.cohort_id}.spatial_coord_by_clusters.png",
    )

    # Save adata object
    adata.write_h5ad(filename=args.counts_output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Visualize annotation on UMAP and spatial coordinates"
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
        help="Adata object to cluster and plot"
    )
    parser.add_argument(
        "-o",
        "--adata-output",
        type=str,
        required=True,
        help="Output file name for the Leiden clustered AnnData object"
    )

    args = parser.parse_args()

    main(args)
