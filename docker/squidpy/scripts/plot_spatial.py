import argparse
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq


def main(args):
    adata = sc.read_h5ad(args.adata_input)

    # UMAP plots
    sc.settings.set_figure_params(
        dpi=100, fontsize=10, dpi_save=300, format="png", figsize=("12", "8")
    )

    sc.pl.embedding(
        adata,
        basis="umap",
        color=["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_rb", "doublet_score"],
        frameon=False,
        show=False,
        ncols=1,
    )
    plt.savefig(f"{args.plots_prefix}.features_umap.png", dpi=300, bbox_inches="tight")

    sc.pl.embedding(
        adata,
        basis="umap",
        color=["cell_type", "leiden_res_0.05", "leiden_res_0.10", "leiden_res_0.20", "leiden_res_0.40"],
        frameon=False,
        show=False,
        ncols=1,
    )
    plt.savefig(f"{args.plots_prefix}.groups_umap.png", dpi=300, bbox_inches="tight")

    # Compare feature and gene clusters
    sq.pl.spatial_scatter(
        adata,
        library_key="sample",
        color=["features_cluster", "cell_type"],
    )
    plt.savefig(f"{args.plots_prefix}.image_features_spatial_scatter.png", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate image features (e.g. morphological, intensity-based, and texture-based features) to help describe spatial patterns that can be linked to gene expression or cell organization."
    )
    parser.add_argument(
        "-i",
        "--adata-input",
        type=str,
        required=True,
        help="Adata object to calculate image features on"
    )
    parser.add_argument(
        "-p",
        "--plots-prefix",
        type=str,
        required=True,
        help="Output file name prefix for the feature and gene clusters spatial scatter plots"
    )

    args = parser.parse_args()

    main(args)
