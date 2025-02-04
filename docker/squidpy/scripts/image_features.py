import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq


# https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_visium_hne.html#image-features

def main(args):
    ##############################
    ## CALCULATE IMAGE FEATURES ##
    ##############################
    adata = sc.read_10x_h5(args.adata_input)

    # Calculate features for different scales (higher value means more context)
    for scale in [1.0, 2.0]:
        feature_name = f"features_summary_scale{scale}"
        sq.im.calculate_image_features(
            adata,
            img_key="hires",
            features="summary",
            key_added=feature_name,
            n_jobs=args.n_jobs,
            scale=scale,
        )


    # Combine features in one dataframe
    adata.obsm["features"] = pd.concat(
        [adata.obsm[f] for f in adata.obsm.keys() if "features_summary" in f],
        axis="columns",
    )
    # Make sure that we have no duplicated feature names in the combined table
    adata.obsm["features"].columns = ad.utils.make_index_unique(
        adata.obsm["features"].columns
    )

    ############################
    ## NEW CLUSTER ANNOTATION ##
    ############################
    # We can use the extracted image features to compute a new cluster annotation
    # This could be useful to gain insights in similarities across spots based on image morphology
    features = adata.obsm["features"].filter(like="summary")
    # Create temporary adata to calculate the clustering
    adata_tmp = ad.AnnData(features)
    # IMPORTANT - feature values are not scaled, so need to scale them before PCA
    sc.pp.scale(adata_tmp)
    sc.pp.pca(adata_tmp, n_comps=min(10, features.shape[1] - 1))
    sc.pp.neighbors(adata_tmp)
    sc.tl.leiden(adata_tmp)
    adata.obs["features_cluster"] = adata_tmp.obs["leiden"]

    # Compare feature and gene clusters
    sq.pl.spatial_scatter(
        adata,
        color=["features_cluster", "cluster"],
        save=f"{args.plots_prefix}.image_features_spatial_scatter.png",
    )

    # Save adata object
    adata.write_h5ad(filename=args.adata_output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate image features (e.g. morphological, intensity-based, and texture-based features) to help describe spatial patterns that can be linked to gene expression or cell organization."
    )
    parser.add_argument(
        "-i",
        "--adata-input",
        type=str,
        required=True,
        help="Adata object to calculate image features on"
    )
    parser.add_argument(
        "-n",
        "--n-jobs",
        type=str,
        required=True,
        help="Number of parallel jobs"
    )
    parser.add_argument(
        "-p",
        "--plots-prefix",
        type=str,
        required=True,
        help="Output file name prefix for the feature and gene clusters spatial scatter plots"
    )
    parser.add_argument(
        "-o",
        "--adata-output",
        type=str,
        required=True,
        help="Output file name for the extracted image features AnnData object"
    )

    args = parser.parse_args()

    main(args)
