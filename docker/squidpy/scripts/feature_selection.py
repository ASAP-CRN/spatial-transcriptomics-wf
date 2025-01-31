import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc


def main(args):
    #######################
    ## FEATURE SELECTION ##
    #######################
    adata = sc.read_10x_h5(args.adata_input)

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=args.n_top_genes,
        batch_key=args.batch_key,
    )

    sc.pl.highly_variable_genes(
        adata,
        save=f"{args.plots_prefix}.feature_dispersion.png",
    )


    ##############################
    ## DIMENSIONALITY REDUCTION ##
    ##############################
    sc.tl.pca(adata)

    # Save outputs
    adata.write_h5ad(filename=args.adata_output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Feature selection and dimensionality reduction adata objects"
    )
    parser.add_argument(
        "-i",
        "--adata-input",
        type=str,
        required=True,
        help="Adata object to reduce the dimensionality of the dataset and only include the most informative genes"
    )
    parser.add_argument(
        "-b",
        "--batch-key",
        type=int,
        required=True,
        help="Key in AnnData object for batch information so that highly-variable genes are selected within each batch separately and merged ['batch_id']"
    )
    parser.add_argument(
        "-n",
        "--n-top-genes",
        type=int,
        required=True,
        help="Number of highly-variable genes to keep using `scanpy.pp.highly_variable_genes [3000]"
    )
    parser.add_argument(
        "-p",
        "--plots-prefix",
        type=str,
        required=True,
        help="Output file name prefix for the highly variable genes and PCA plots"
    )
    parser.add_argument(
        "-o",
        "--adata-output",
        type=str,
        required=True,
        help="Output file name for the feature selected and dimensionality reduced AnnData object"
    )

    args = parser.parse_args()

    main(args)
