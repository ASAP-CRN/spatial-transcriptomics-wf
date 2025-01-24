import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq


def main(args):
    adata = sc.read_h5ad(args.adata_input)
    adata.var["NegPrb"] = adata.var_names.str.startswith("NegPrb")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["NegPrb"], inplace=True)
    pd.set_option("display.max_columns", None)

    # Calculate the percentage of unassigned "NegPrb" transcripts
    unassigned_ctrl_probes_percentage = adata.obs["total_counts_NegPrb"].sum() / adata.obs["total_counts"].sum() * 100

    # QC distribution plots
    fig, axs = plt.subplots(1, 3, figsize=(15, 4))
    axs[0].set_title("Total transcripts per cell")
    sns.histplot(
        adata.obs["total_counts"],
        kde=False,
        ax=axs[0],
    )
    axs[1].set_title("Unique transcripts per cell")
    sns.histplot(
        adata.obs["n_genes_by_counts"],
        kde=False,
        ax=axs[1],
    )
    axs[2].set_title("Transcripts per FOV")
    sns.histplot(
        adata.obs.groupby("fov").sum()["total_counts"],
        kde=False,
        ax=axs[2],
    )

    # Save outputs
    with open("unassigned_ctrl_probes_percentage.txt", "w") as file:
        file.write(f"{unassigned_ctrl_probes_percentage}")

    plt.tight_layout()
    fig.savefig(args.qc_plots_output, dpi=300)

    adata.write_h5ad(filename=args.qc_adata_output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate quality control metrics"
    )
    parser.add_argument(
        "-s",
        "--sample-id",
        type=str,
        required=True,
        help="Sample ID"
    )
    parser.add_argument(
        "-i",
        "--adata-input",
        type=str,
        required=True,
        help="Preprocessed adata object to run QC on"
    )
    parser.add_argument(
        "-p",
        "--qc-plots-output",
        type=str,
        required=True,
        help="Output file name for the QC distribution plots"
    )
    parser.add_argument(
        "-o",
        "--qc-adata-output",
        type=str,
        required=True,
        help="Output file name for the QC'ed AnnData object"
    )

    args = parser.parse_args()

    main(args)
