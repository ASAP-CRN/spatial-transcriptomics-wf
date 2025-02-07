import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc


def main(args):
    ########
    ## QC ##
    ########
    adata = sc.read_h5ad(args.adata_input)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
    )

    # Add doublet_score and predicted_doublet to .obs
    sc.pp.scrublet(adata)

    # Adata object must contain strings
    for col in adata.obs.columns:
        if pd.api.types.is_object_dtype(adata.obs[col].dtype):
            adata.obs[col] = adata.obs[col].astype(str) 
    
    adata.write_h5ad(filename=args.qc_adata_output)


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
