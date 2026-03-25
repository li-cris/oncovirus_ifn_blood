import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import os
import random
# import harmonypy as hm

os.chdir("/Users/cristinali/Documents/Programming/PhD/onco_ifn")

from src.sctype.sctype_py import gene_sets_prepare, sctype_score, process_cluster, get_gene_symbols
from src.visualization.umap import plot_umap_with_subset_percentages

# Objective: Do the necessary processing of the individual AnnData objects to prepare them for integration, and then perform the integration using Scanpy's integration functions.

# Get list of GSM IDs from the directory
adata_path = "data/GSE310935_adata.h5ad"
adata = ad.read_h5ad(adata_path)

# --------------------------------- #
# Identifying what they do in the paper
# 1. Normalize raw gene expression: expression valyes for each gene normalized to the total expression in each cell, multiplied by a scale factor of 10,000, and followed by log-transformation
    # This is the standard normalization procedure in Scanpy, which we can do with sc.pp.normalize_total and sc.pp.log1p

# 2. Scaling data to center and scale each gene to have a mean of 0 and a variance of 1

# 3. Integration of individual samples (STS baseline, STS follow-up, LTS baseline and LTS follow-up)

# 4. Cluster annotation across integrated objects by single cell type

# 5. Quantification of abundance of each cell type in each sample, and comparison of abundance between STS and LTS groups at baseline and follow-up time points

# 6. If possible, quantify activation scores based on a list of genes:
    # "CD69", "CD44", "ICOS", "CD40LG", "TNF", "IFNG", "GZMB", "PRF1", "IL2", "FOS", "JUN", "CXCR4", "BATF"
# --------------------------------- #



# Integration and all that
# --------------------------------- #
# 1. Some preprocessing and adding batch information
# --------------------------------- #
# Identify the batches
# adata.obs["batch"] = adata.obs["LTS"].astype(str) + "_" + adata.obs["timepoint"].astype(str) # Create a batch column for integration, combining LTS/STS and baseline/follow-up
adata.obs["batch"] = adata.obs["og_id"].astype(str)
# Save the original counts in a separate layer before normalization, so we can use it for downstream analysis if needed
adata.layers["counts"] = adata.X.copy()
# Turn adata.X into a dense matrix for processing
# adata.X = adata.X.toarray()
# Up till now:
    # AnnData object with n_obs × n_vars = 34482 × 36601

# If can't find approved genes, do this
if not os.path.exists("results/approved_genes.csv"):
    list_of_genes = adata.var["HUGO_symbol"].values
    approved_genes = get_gene_symbols(list_of_genes) # Gives a DF, with "Symbol" being the approved name in the 2nd column
    approved_genes.to_csv("results/approved_genes.csv", sep = "\t", index=False)
else:
    approved_genes = pd.read_csv("results/approved_genes.csv", sep = "\t")

# --------------------------------- #
# 1.5. Filter to baseline or follo-up only, since that's what they do in the paper
# --------------------------------- #
# b or f
random.seed(42)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata = adata[adata.obs["timepoint"] == "b"].copy()

# --------------------------------- #
# 2. Filter out cells or genes with low counts (these are raw counts)
# --------------------------------- #

# After filtering:
    # AnnData object with n_obs × n_vars = 33185 × 25386
    # More genes were filtered out than cells, which is expected since many genes are not expressed in many cells


# --------------------------------- #
# 3. Normalization and log transformation, filtering highly variable genes (HVGs)
# --------------------------------- #
# Usual normalization and log transformation
sc.pp.normalize_total(adata,target_sum=1e4)
sc.pp.log1p(adata)
# pip install --user scikit-misc
# Use Seurat's method, but in Python since it's easier to work with anndata
sc.pp.highly_variable_genes(adata, n_top_genes=4000, flavor="seurat_v3",layer="counts", batch_key="batch")
adata.raw = adata
adata = adata[:, adata.var["highly_variable"]].copy()

# --------------------------------- #
# 4. Scale and run PCA
# --------------------------------- #
sc.pp.scale(adata,max_value=10)

# This part is for the ScType annotation later
scaled_data = pd.DataFrame(adata.X) # It works with DF for some reason
# change column to gene names

# Turn the current list of gene names into the approved version
# If a gene in approved list is not in adata, ignore it, but if a gene in adata is not in approved list, keep it as is
gene_name_mapping = dict(zip(approved_genes["Gene"], approved_genes["Symbol"]))
adata.var["approved_gene"] = adata.var["HUGO_symbol"].map(gene_name_mapping).fillna(adata.var["HUGO_symbol"])
scaled_data.columns = adata.var["approved_gene"].values
# Change the row cell names
scaled_data.index = adata.obs_names
scaled_data=scaled_data.T


# PCA
sc.tl.pca(adata)
# sc.tl.pca(adata,zero_center=False)

# --------------------------------- #
# 5. Harmony integration, UMAP and clustering
# --------------------------------- #

# Run Harmony and ensure the embedding has one row per cell (n_obs x n_pcs).
import harmonypy as hm
ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, "batch", max_iter_harmony=20, max_iter_kmeans=20, sigma=0.1, theta=2, nclust=50, verbose=True)
x_pca_harmony = np.asarray(ho.Z_corr)

if x_pca_harmony.shape[0] != adata.n_obs and x_pca_harmony.shape[1] == adata.n_obs:
    x_pca_harmony = x_pca_harmony.T

if x_pca_harmony.shape[0] != adata.n_obs:
    raise ValueError(
        f"Harmony output has unexpected shape {x_pca_harmony.shape}; expected first dimension {adata.n_obs}."
    )

adata.obsm["X_pca_harmony"] = x_pca_harmony

# sc.pp.neighbors(adata, n_neighbors=20, n_pcs=10,use_rep="X_pca")
sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=20, n_pcs=10)
sc.tl.umap(adata,min_dist=0.3)
sc.tl.leiden(adata, resolution=0.8)
sc.pl.umap(
    adata,
    color=["leiden", "patient", "LTS", "timepoint", "batch"],
    wspace=0.4
)
#sc.pl.umap(adata, color=['leiden'])


# --------------------------------- #
# 6. Python implementation of ScType
# --------------------------------- #

# 5. Assign with ScType annotator
scRNAseqData=scaled_data
gs_list=gene_sets_prepare(path_to_db_file="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",cell_type="Immune system")
# Just in case, also update identifiers in scRNAseqData to match updates gene names
es_max = sctype_score(scRNAseqData = scRNAseqData, scaled = True, gs = gs_list['gs_positive'], gs2 = gs_list['gs_negative'])

unique_clusters = adata.obs['leiden'].unique() # are we sure it's leiden?
# Apply the function to each unique cluster and combine the results into a DataFrame
cL_results = pd.concat([process_cluster(cluster,adata,es_max,'leiden') for cluster in unique_clusters])

# Group by cluster and select the top row based on scores
sctype_scores = cL_results.groupby('cluster').apply(lambda x: x.nlargest(1, 'scores')).reset_index(drop=True)

# Set low-confidence clusters to "Unknown"
sctype_scores.loc[sctype_scores['scores'] < sctype_scores['ncells'] / 4, 'type'] = 'Unknown'

# Iterate over unique clusters
adata.obs['sctype_classification'] = ""
for cluster in sctype_scores['cluster'].unique():
    # Filter sctype_scores for the current cluster
    cl_type = sctype_scores[sctype_scores['cluster'] == cluster]
    # Get the type for the current cluster
    cl_type_value = cl_type['type'].iloc[0]
    # Update 'sctype_classification' in pbmc.obs for cells belonging to the current cluster
    adata.obs.loc[adata.obs['leiden'] == cluster, 'sctype_classification'] = cl_type_value

# Plot the UMAP with sctype_classification as labels
# Add percentage of each cell type in each cluster to legend

# celltype_counts = adata.obs['sctype_classification'].value_counts()
# celltype_percent = celltype_counts / celltype_counts.sum() * 100

# ctab = pd.crosstab(
#     adata.obs['LTS'],
#     adata.obs['sctype_classification'],
#     normalize='index'
# ) * 100
# ctab.plot(kind='bar', stacked=True, figsize=(8,5))

# sc.pl.umap(adata, color='sctype_classification', title='UMAP with sctype_classification')


# --------------------------------- #
# 7. Plot UMAP with scType classification, adding percentage of each cell type in the legend
# --------------------------------- #

# Count percentages of each annotated cell type
# perc = (
#     adata.obs["sctype_classification"]
#     .value_counts(normalize=True)
#     .mul(100)
#     .round(1)
# )

# # Build new labels like "CD4 T cell (35.2%)"
# new_labels = {
#     ct: f"{ct} ({perc[ct]}%)"
#     for ct in perc.index
# }

# # Create a new column just for plotting
# adata.obs["sctype_legend"] = (
#     adata.obs["sctype_classification"]
#     .map(new_labels)
#     .astype("category")
# )

# # Optional: keep legend order consistent with abundance
# ordered_cats = [new_labels[ct] for ct in perc.index]
# adata.obs["sctype_legend"] = adata.obs["sctype_legend"].cat.set_categories(
#     ordered_cats
# )

# # Plot 2 UMAPS, one filtering by adata.obs["LTS"] == 1 and the other by adata.obs["LTS"] == 0, to compare the distribution of cell types between LTS and STS

# sc.pl.umap(
#     adata,
#     color="sctype_legend",
#     title="UMAP with scType classification",
#     legend_loc="right margin"
# )


#%%
import seaborn as sns
celltypes = adata.obs["sctype_classification"].value_counts().index.tolist()
palette = sns.color_palette("tab20", n_colors=len(celltypes))
color_map = dict(zip(celltypes, palette))

# import random
# random.seed(42)
plot_umap_with_subset_percentages(adata, "LTS", 1, label_colormap=color_map, label_types=celltypes)
plot_umap_with_subset_percentages(adata, "LTS", 0, label_colormap=color_map, label_types=celltypes)


# %% NO
# 6. HVGs, PCA, scaling
sc.pp.highly_variable_genes(adata, batch_key="batch", n_top_genes=2000)
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata)
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

sc.pl.umap(adata, color="batch")

# 7a. Integration option: BBKNN
import scanpy.external as sce
sce.pp.bbknn(adata, batch_key="batch")


# 8. UMAP and clustering
sc.tl.umap(adata)
sc.tl.leiden(adata)

# 9. Plot
sc.pl.umap(adata, color=["batch", "leiden"])



# Trying to query for approved genes in database

# import requests

# def hgnc_lookup(term):
#     url = f"https://rest.genenames.org/search/{term}"
#     headers = {"Accept": "application/json"}
#     r = requests.get(url, headers=headers, timeout=30)
#     r.raise_for_status()
#     return r.json()

# term = "TNFRSF8"
# data = hgnc_lookup(term)
# name = data['response']['docs'][0]['symbol'] if data['response']['docs'] else None

# if name:
#     # Use the gene symbol
# else:
#     # Remove said term from the list

# if name:
#     print(f"Found gene symbol: {name}")
# else:
#     print(f"No gene symbol found for the given term: {term}")
