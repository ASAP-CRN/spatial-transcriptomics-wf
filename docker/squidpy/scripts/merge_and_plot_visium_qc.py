import argparse
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq


def main(args):
    #########################
    ## MERGE ADATA OBJECTS ##
    #########################
    adatas = {}
    for file_path in args.adata_paths_input:
        adata = sc.read_10x_h5(file_path)
        sample_id = adata.obs["sample"].unique()
        adatas[sample_id[0]] = adata

    merged_adata = ad.concat(adatas, index_unique="_", merge="same", uns_merge="same")


    #####################
    ## PLOT QC METRICS ##
    #####################
    # QC violin plots
    sc.pl.violin(
        merged_adata,
        [
            "n_genes_by_counts",
            "total_counts",
            "pct_counts_mt"
        ],
        jitter=0.4,
        multi_panel=True,
        save=f"{args.qc_plots_prefix}.qc_violin.png",
    )

    # QC scatter plot colored by pct_counts_mt
    sc.pl.scatter(
        merged_adata,
        "total_counts",
        "n_genes_by_counts",
        color="pct_counts_mt",
        save=f"{args.qc_plots_prefix}.qc_scatter.png",
    )

    # Save outputs
    merged_adata.write_h5ad(filename=args.merged_adata_output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge adata objects and plot quality control metrics"
    )
    parser.add_argument(
        "-i",
        "--adata-paths-input",
        nargs="+",
        required=True,
        help="List of preprocessed adata objects to merge"
    )
    parser.add_argument(
        "-p",
        "--qc-plots-prefix",
        type=str,
        required=True,
        help="Output file name prefix for the QC violin and scatter plots"
    )
    parser.add_argument(
        "-o",
        "--merged-adata-output",
        type=str,
        required=True,
        help="Output file name for the merged AnnData object"
    )

    args = parser.parse_args()

    main(args)
