import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc


def main(args):
    ##########################
    ## FILTER ADATA OBJECTS ##
    ##########################
    adata = sc.read_h5ad(args.adata_input)

    sc.pp.filter_cells(adata, min_counts=args.min_counts)
    print(f"Number of cells after min count filter ({args.min_counts}): {adata.n_obs}")

    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    print(f"Number of genes after cell filter ({args.min_cells}): {adata.n_vars}")

    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    print(f'Number of cells after gene filter ({args.min_genes}): {adata.n_obs}')

    adata = adata[adata.obs["mt_frac"] < args.mt_max_percent]
    print(f"Number of cells after MT filter ({args.mt_max_percent}): {adata.n_obs}")


    ###############
    ## NORMALIZE ##
    ###############
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(
        adata,
        target_sum=args.target_sum,
        inplace=True,
    )
    sc.pp.log1p(adata)


    ##################
    ## ANNOTATE HVG ##
    ##################
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
        n_comps=args.n_comps,
        svd_solver="arpack",
    )

    # Save outputs
    adata.write_h5ad(filename=args.adata_output, compression="gzip")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process adata objects by filtering, normalizing, identifying highly variable genes (HVGs), and PCA"
    )
    parser.add_argument(
        "-i",
        "--adata-input",
        type=str,
        required=True,
        help="Adata object to process"
    )
    parser.add_argument(
        "-m",
        "--min-counts",
        type=int,
        required=True,
        help="Minimum number of counts for filtering using `scanpy.pp.filter_cells` [5000]"
    )
    parser.add_argument(
        "-g",
        "--min-genes",
        type=int,
        required=True,
        help="Minimum number of genes for filtering using `scanpy.pp.filter_cells` [3000]"
    )
    parser.add_argument(
        "-l",
        "--min-cells",
        type=int,
        required=True,
        help="Minimum number of cells for filtering using `scanpy.pp.filter_gene` [10]"
    )
    parser.add_argument(
        "-a",
        "--mt-max-percent",
        type=float,
        required=True,
        help="Maximum percentage of mitochondrial read counts for filtering [0.2]"
    )
    parser.add_argument(
        "-t",
        "--target-sum",
        type=float,
        required=True,
        help="The total count to which each cell's gene expression values will be normalized using `scanpy.pp.normalize_total` [10000]"
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
        help="Output file name for the processed AnnData object"
    )

    args = parser.parse_args()

    main(args)
