import pandas as pd
import scanpy as sc
import numpy as np
import harmonypy as hm

def run_harmony(adata,
                batch_key="batch",
                max_iter_harmony=20,
                max_iter_kmeans=20,
                sigma=0.1,
                theta=2,
                nclust=50,
                verbose=True):

    ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, batch_key, max_iter_harmony=max_iter_harmony, max_iter_kmeans=max_iter_kmeans, sigma=sigma, theta=theta, nclust=nclust, verbose=verbose)
    x_pca_harmony = np.asarray(ho.Z_corr)

    if x_pca_harmony.shape[0] != adata.n_obs and x_pca_harmony.shape[1] == adata.n_obs:
        x_pca_harmony = x_pca_harmony.T

    if x_pca_harmony.shape[0] != adata.n_obs:
        raise ValueError(
            f"Harmony output has unexpected shape {x_pca_harmony.shape}; expected first dimension {adata.n_obs}."
        )

    return x_pca_harmony