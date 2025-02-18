import argparse
import pandas as pd
import numpy as np
import scanpy as sc


def main(args):
    ##########################
    ## FILTER ADATA OBJECTS ##
    ##########################
    adata = sc.read_h5ad(args.adata_input)

    sc.pp.filter_cells(adata, min_counts=args.min_counts)
    print(f"Number of cells after min count filter ({args.min_counts}): {adata.n_obs}")

    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    print(f"Number of genes after cell filter ({args.min_cells}): {adata.n_vars}")

    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    print(f'Number of cells after gene filter ({args.min_genes}): {adata.n_obs}')

    adata = adata[adata.obs["mt_frac"] < args.mt_max_percent]
    print(f"Number of cells after MT filter ({args.mt_max_percent}): {adata.n_obs}")


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
        "-g",
        "--min-genes",
        type=int,
        required=True,
        help="Minimum number of genes required using `scanpy.pp.filter_cells` [3000]"
    )
    parser.add_argument(
        "-n",
        "--min-cells",
        type=int,
        required=True,
        help="Minimum number of cells required using `scanpy.pp.filter_gene` [10]"
    )
    parser.add_argument(
        "-p",
        "--mt-max-percent",
        type=float,
        required=True,
        help="Maximum percentage of mitochondrial read counts [0.2]"
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
