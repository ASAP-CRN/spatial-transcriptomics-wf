import argparse
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc


def main(args):
    #########################
    ## MERGE ADATA OBJECTS ##
    #########################
    adatas = {}
    for file_path in args.adata_paths_input:
        adata = sc.read_h5ad(file_path)
        sample_id = adata.obs["sample"].unique()
        adatas[sample_id[0]] = adata

    merged_adata = ad.concat(adatas, index_unique="_", merge="same", uns_merge="unique")


    ##################
    ## ANNOTATE HVG ##
    ##################
    sc.pp.log1p(merged_adata)
    
    sc.pp.highly_variable_genes(
        merged_adata,
        flavor="seurat",
        n_top_genes=args.n_top_genes,
        inplace=True,
    )
    sc.pl.highly_variable_genes(
        merged_adata,
    )
    fig = plt.gcf()
    fig.suptitle(f"Highly variable genes dispersion plot - {args.plots_prefix}", va="center", ha="center", fontsize=16)
    plt.savefig(f"{args.plots_prefix}.hvg_dispersion.png", dpi=300, bbox_inches="tight")


    ##############################
    ## DIMENSIONALITY REDUCTION ##
    ##############################
    sc.pp.pca(
        merged_adata,
        n_comps=args.n_comps,
        svd_solver="arpack",
    )

    # Save outputs
    merged_adata.write_h5ad(filename=args.adata_output, compression="gzip")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge adata objects and process to prepare for downstream analysis"
    )
    parser.add_argument(
        "-i",
        "--adata-paths-input",
        nargs="+",
        required=True,
        help="List of preprocessed adata objects to merge"
    )
    parser.add_argument(
        "-n",
        "--n-top-genes",
        type=int,
        required=True,
        help="Number of HVGs to keep using `scanpy.pp.highly_variable_genes` [3000]"
    )
    parser.add_argument(
        "-c",
        "--n-comps",
        type=int,
        required=True,
        help="Number of principal components to compute using `scanpy.pp.pca` [30]"
    )
    parser.add_argument(
        "-p",
        "--plots-prefix",
        type=str,
        required=True,
        help="Output file name prefix for the HVGs plot"
    )
    parser.add_argument(
        "-o",
        "--adata-output",
        type=str,
        required=True,
        help="Output file name for the merged AnnData object"
    )

    args = parser.parse_args()

    main(args)
