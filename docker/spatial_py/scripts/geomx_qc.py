import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
import scanpy as sc


def main(args):
    ########
    ## QC ##
    ########
    adata = sc.read_h5ad(args.adata_input)
    adata.var["NegPrb"] = adata.var_names.str.startswith("NegProbe")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["NegPrb"], inplace=True)
    pd.set_option("display.max_columns", None)

    # Calculate the percentage of unassigned "NegPrb" transcripts
    unassigned_ctrl_probes_percentage = adata.obs["total_counts_NegPrb"].sum() / adata.obs["total_counts"].sum() * 100

    # Save outputs
    with open("unassigned_ctrl_probes_percentage.txt", "w") as file:
        file.write(f"{unassigned_ctrl_probes_percentage}")

    adata.write_h5ad(filename=args.qc_adata_output, compression="lzf")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate quality control metrics on adata objects"
    )
    parser.add_argument(
        "-i",
        "--adata-input",
        type=str,
        required=True,
        help="Adata object to run QC on"
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
