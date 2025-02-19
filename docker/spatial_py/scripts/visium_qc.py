import argparse
import pandas as pd
import numpy as np
import scipy as sp
import scanpy as sc

# https://scanpy-tutorials.readthedocs.io/en/multiomics/analysis-visualization-spatial.html


def main(args):
    ########
    ## QC ##
    ########
    adata = sc.read_h5ad(args.adata_input)
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["rb"] = adata.var_names.str.startswith(("RPS", "RPL"))
    adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

    adata.obs["mt_frac"] = adata[:, adata.var["mt"]].X.sum(1).A.squeeze()/adata.obs["total_counts"]

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "rb", "hb"], inplace=True, log1p=True
    )

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
