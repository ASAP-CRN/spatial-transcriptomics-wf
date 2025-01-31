import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc


def main(args):
    ##########################
    ## FILTER ADATA OBJECTS ##
    ##########################
    adata = sc.read_h5ad(args.adata_input)

    sc.pp.filter_cells(adata, min_counts=args.min_counts)
    sc.pp.filter_genes(adata, min_cells=args.min_cells)


    ###############
    ## NORMALIZE ##
    ###############
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)

    # Save outputs
    adata.write_h5ad(filename=args.adata_output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter and normalize adata objects"
    )
    parser.add_argument(
        "-i",
        "--adata-input",
        type=str,
        required=True,
        help="Adata object to filter and normalize"
    )
    parser.add_argument(
        "-m",
        "--min-counts",
        type=int,
        required=True,
        help="Minimum number of counts required using `scanpy.pp.filter_cells` [5000]"
    )
    parser.add_argument(
        "-n",
        "--min-cells",
        type=int,
        required=True,
        help="Minimum number of cells required using `scanpy.pp.filter_gene` [10]"
    )
    parser.add_argument(
        "-o",
        "--adata-output",
        type=str,
        required=True,
        help="Output file name for the filtered and normalized AnnData object"
    )

    args = parser.parse_args()

    main(args)
