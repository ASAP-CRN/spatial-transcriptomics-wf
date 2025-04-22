import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq


def main(args):
    adata = sc.read_h5ad(args.adata_input)

    samples = adata.obs["sample"].unique().tolist()
    colors = ["total_counts", "n_genes_by_counts", "batch", "leiden"]
    titles = [f"{sample} - {color}" for sample in samples for color in colors]

    sq.pl.spatial_scatter(
        adata,
        library_key="sample",
        color=colors,
        title=titles,
    )
    plt.savefig(f"{args.plots_prefix}.spatial_scatter.png", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot spatial omics data with data overlayed on top"
    )
    parser.add_argument(
        "-i",
        "--adata-input",
        type=str,
        required=True,
        help="Adata object to plot"
    )
    parser.add_argument(
        "-p",
        "--plots-prefix",
        type=str,
        required=True,
        help="Output file name prefix for the spatial scatter plots"
    )

    args = parser.parse_args()

    main(args)
