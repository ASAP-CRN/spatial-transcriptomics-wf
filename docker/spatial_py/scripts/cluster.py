import argparse
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scanpy as sc


def main(args):
    ##########################
    ## CLUSTER AND RUN UMAP ##
    ##########################
    adata = sc.read_h5ad(args.adata_input)

    sc.pp.neighbors(
        adata,
        n_pcs=args.n_comps,
    )
    sc.tl.umap(adata)

    # Visualize two covariates (total counts per spot and number of genes by counts) in UMAP space
    sc.tl.leiden(
        adata,
        resolution=args.resolution,
    )

    sc.pl.umap(
        adata,
        color=["total_counts", "n_genes_by_counts", "batch", "leiden"],
    )
    plt.title(f"UMAP - {args.plots_prefix}")
    plt.savefig(f"{args.plots_prefix}.umap_cluster.png", dpi=300, bbox_inches="tight")

    # Save outputs
    adata.write_h5ad(filename=args.adata_output, compression="gzip")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Perform clustering on adata objects and plot"
    )
    parser.add_argument(
        "-i",
        "--adata-input",
        type=str,
        required=True,
        help="Adata object to cluster"
    )
    parser.add_argument(
        "-c",
        "--n-comps",
        type=int,
        required=True,
        help="Number of principal components to use using `scanpy.pp.neighbors` [30]"
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=float,
        required=True,
        help="Value controlling the coarseness of the Leiden clustering (higher values lead to more clusters) using `scanpy.tl.leiden` [0.4]"
    )
    parser.add_argument(
        "-p",
        "--plots-prefix",
        type=str,
        required=True,
        help="Output file name prefix for the UMAP plots"
    )
    parser.add_argument(
        "-o",
        "--adata-output",
        type=str,
        required=True,
        help="Output file name for the clustered AnnData object"
    )

    args = parser.parse_args()

    main(args)
