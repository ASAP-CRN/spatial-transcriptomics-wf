import argparse
import pandas as pd
import numpy as np
import scanpy as sc
import harmonypy as hm

# https://stlearn.readthedocs.io/en/latest/tutorials/Integration_multiple_datasets.html


def main(args):
    ##################################
    ## RUN INTEGRATION WITH HARMONY ##
    ##################################
    adata = sc.read_h5ad(args.adata_input)

    sc.external.pp.harmony_integrate(
        adata,
        key=args.batch_key,
    )

    # Save outputs
    adata.write_h5ad(filename=args.adata_output, compression="gzip")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run integration with Harmony on adata objects"
    )
    parser.add_argument(
        "-i",
        "--adata-input",
        type=str,
        required=True,
        help="Adata object to run Harmony on"
    )
    parser.add_argument(
        "-b",
        "--batch-key",
        type=str,
        required=True,
        help="Key in AnnData object for batch information ['batch_id']"
    )
    parser.add_argument(
        "-o",
        "--adata-output",
        type=str,
        required=True,
        help="Output file name for the Harmony integrated AnnData object"
    )

    args = parser.parse_args()

    main(args)
