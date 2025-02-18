import argparse
import pandas as pd
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scanpy as sc


def main(args):
    ##################
    ## ANNOTATE HVG ##
    ##################
    adata = sc.read_h5ad(args.adata_input)

    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat",
        n_top_genes=args.n_top_genes,
        inplace=True,
    )
    sc.pl.highly_variable_genes(
        adata,
    )
    plt.savefig(f"{args.plots_prefix}.hvg_dispersion.png", dpi=300, bbox_inches="tight")


    ##############################
    ## DIMENSIONALITY REDUCTION ##
    ##############################
    sc.pp.pca(
        adata,
        n_comps=50,
        use_highly_variable=True,
        svd_solver="arpack",
    )
    sc.pp.neighbors(adata)


    ##################
    ## UMAP CLUSTER ##
    ##################
    # Visualize two covariates (total counts per spot and number of genes by counts) in UMAP space
    sc.tl.umap(adata)
    sc.tl.louvain(adata, key_added="clusters")

    fig, axs = plt.subplots(1, 3, figsize=(15,4))
    sc.pl.umap(
        adata,
        color="total_counts",
        ax=axs[0],
    )
    sc.pl.umap(
        adata,
        color="n_genes_by_counts",
        ax=axs[1],
    )
    sc.pl.umap(
        adata,
        color="clusters",
        palette=sc.pl.palettes.default_20,
        ax=axs[2],
    )
    plt.savefig(f"{args.plots_prefix}.umap_cluster.png", dpi=300, bbox_inches="tight")

    # Save outputs
    adata.write_h5ad(filename=args.adata_output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Annotate highly variable genes (HVGs) and dimensionality reduction on adata objects"
    )
    parser.add_argument(
        "-i",
        "--adata-input",
        type=str,
        required=True,
        help="Adata object to reduce the dimensionality of the dataset and only include the most informative genes"
    )
    parser.add_argument(
        "-n",
        "--n-top-genes",
        type=int,
        required=True,
        help="Number of HVGs to keep using `scanpy.pp.highly_variable_genes [3000]"
    )
    parser.add_argument(
        "-p",
        "--plots-prefix",
        type=str,
        required=True,
        help="Output file name prefix for the highly variable genes and UMAP plots"
    )
    parser.add_argument(
        "-o",
        "--adata-output",
        type=str,
        required=True,
        help="Output file name for the HVGs selected and dimensionality reduced AnnData object"
    )

    args = parser.parse_args()

    main(args)
