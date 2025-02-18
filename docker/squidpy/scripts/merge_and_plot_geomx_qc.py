import argparse
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
import scanpy as sc
import squidpy as sq


def main(args):
    #########################
    ## MERGE ADATA OBJECTS ##
    #########################
    adatas = {}
    for file_path in args.adata_paths_input:
        adata = sc.read_h5ad(file_path)
        sample_id = adata.obs["sample"].unique()
        adatas[sample_id[0]] = adata

    merged_adata = ad.concat(adatas, index_unique="_", merge="same", uns_merge="same")


    #####################
    ## PLOT QC METRICS ##
    #####################
    # QC distribution plots
    fig, axs = plt.subplots(1, 3, figsize=(15, 4))
    axs[0].set_title("Total transcripts per cell")
    sb.histplot(
        merged_adata.obs["total_counts"],
        kde=False,
        ax=axs[0],
    )
    axs[1].set_title("Unique transcripts per cell")
    sb.histplot(
        merged_adata.obs["n_genes_by_counts"],
        kde=False,
        ax=axs[1],
    )
    fov_values = merged_adata.obs.get("fov")
    if fov_values is not None:
        axs[2].set_title("Transcripts per FOV")
        sb.histplot(
            merged_adata.obs.groupby("fov").sum()["total_counts"],
            kde=False,
            ax=axs[2],
        )

    # Save outputs
    plt.tight_layout()
    fig.savefig(args.qc_plots_output, dpi=300)

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
        "-o",
        "--merged-adata-output",
        type=str,
        required=True,
        help="Output file name for the merged AnnData object"
    )
    parser.add_argument(
        "-p",
        "--qc-plots-output",
        type=str,
        required=True,
        help="Output file name for the QC distribution plots"
    )

    args = parser.parse_args()

    main(args)
