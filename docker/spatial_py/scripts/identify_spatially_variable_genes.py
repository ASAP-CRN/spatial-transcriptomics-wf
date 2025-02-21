import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sb
import scanpy as sc
import squidpy as sq


def main(args):
    ##########################################################
    ## MORAN'S I GLOBAL SPATIAL AUTO-CORRELATION STATISTICS ##
    ##########################################################
    adata = sc.read_h5ad(args.adata_input)

    genes = adata[:, adata.var.highly_variable]
    sq.gr.spatial_neighbors(
        adata,
        library_key="sample",
        coord_type="generic",
        delaunay=True,
    )
    sq.gr.spatial_autocorr(
        adata,
        mode="moran",
        genes=genes,
        n_perms=100,
        n_jobs=1,
    )
    top_10_variable_genes = adata.uns["moranI"].head(10)

    top_3_variable_gene_list = adata.uns["moranI"].head(3).index.tolist()
    top_3_variable_gene_list.append("leiden")
    sq.pl.spatial_scatter(
        adata,
        library_key="sample",
        color=top_3_variable_gene_list,
    )
    plt.savefig(f"{args.cohort_id}.moran_top_3_variable_genes_spatial_scatter.png", dpi=300, bbox_inches="tight")

    # Save table and adata object
    top_10_variable_genes.to_csv(f"{args.cohort_id}.moran_top_10_variable_genes.csv")
    adata.write_h5ad(filename=args.adata_output, compression="gzip")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Identify spatially variable genes by computing Moran's I Score"
    )
    parser.add_argument(
        "-c",
        "--cohort-id",
        type=str,
        required=True,
        help="Cohort ID"
    )
    parser.add_argument(
        "-i",
        "--adata-input",
        type=str,
        required=True,
        help="Adata object to compute Moran's I Score on"
    )
    parser.add_argument(
        "-o",
        "--adata-output",
        type=str,
        required=True,
        help="Output file name for the spatial auto-correction AnnData object"
    )

    args = parser.parse_args()

    main(args)
