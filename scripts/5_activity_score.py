import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import os
import random
import sys
import seaborn as sns
# import harmonypy as hm

os.chdir("/Users/cristinali/Documents/Programming/PhD/onco_ifn")

# Allow importing local packages from src/ without pip installation (temporary solution for development)
src_path = os.path.abspath("src")
if src_path not in sys.path:
    sys.path.insert(0, src_path)

from src.sctype.sctype_py import gene_sets_prepare, sctype_score, process_cluster, get_gene_symbols
from src.visualization.umap import plot_umap_with_subset_percentages
from src.integration.tools import run_harmony
from src.utils import check_missing_genes
from signaturescoring import score_signature

# Objective: Score cells for a specific gene signature using local signaturescoring package
# Load AnnData
adata_path = "data/GSE310935_adata.h5ad"
adata = ad.read_h5ad(adata_path)

# ============================================================================
# Step 1: Prepare approved gene mapping
# ============================================================================
if not os.path.exists("results/approved_genes.csv"):
    list_of_genes = adata.var["HUGO_symbol"].values
    approved_genes = get_gene_symbols(list_of_genes)  # Returns DF with approved symbols
    approved_genes.to_csv("results/approved_genes.csv", sep="\t", index=False)
else:
    approved_genes = pd.read_csv("results/approved_genes.csv", sep="\t")

# Map gene names to approved symbols
gene_name_mapping = dict(zip(approved_genes["Gene"], approved_genes["Symbol"]))
adata.var["approved_gene"] = adata.var["HUGO_symbol"].map(gene_name_mapping).fillna(adata.var["HUGO_symbol"])

# ============================================================================
# Step 2: Filter and normalize (Should be the same as in 2_anndata_integration.py)
# ============================================================================
adata.obs["batch"] = adata.obs["og_id"].astype(str)
adata.layers["counts"] = adata.X.copy()

random.seed(42)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

# Filter to timepoint "b" only
adata = adata[adata.obs["timepoint"] == "b"].copy()

# Normalize and log-transform
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Highly variable genes using Seurat_v3 method
sc.pp.highly_variable_genes(adata, n_top_genes=4000, flavor="seurat_v3", layer="counts", batch_key="batch")
adata.raw = adata
adata = adata[:, adata.var["highly_variable"]].copy()

# Scale data
sc.pp.scale(adata, max_value=10)

# Use approved gene names as var_names for scoring
adata.var_names = adata.var["approved_gene"].values

# ============================================================================
# Step 3: Score gene expression signature
# ============================================================================
# Verify genes are available, score_signature requires all given genes to be present
gene_list = ["CD69", "CD44", "ICOS", "CD40LG", "TNF", "IFNG", "GZMB", "PRF1", "IL2", "FOS", "JUN", "CXCR4", "BATF"]
gene_list = check_missing_genes(gene_list, adata.var_names)

# Make sure adata.val_names are unique, but for now just remove the duplicates REVIEW
adata = adata[:, ~adata.var_names.duplicated()].copy()

score_signature(
    adata=adata,
    gene_list=gene_list,
    method='seurat_scoring',
    ctrl_size=100,
    score_name='signature_score',
    use_raw=False # This is is because adata.raw contains original gene names which are still ENSEMBL
)

# Print results
print("\nSignature scoring complete!")
print(adata.obs['signature_score'].describe())
print(f"\nScores stored in: adata.obs['signature_score']")

# Okay so adata.obs["signature_score"] has a signature score for each cell
# Now I should be able to plot the overall score for each cell type in violin plots for instance?

# %%
# ============================================================================
# Step 4: Annotation of clusters as well as visualization of scores across cell types
# ============================================================================
scaled_data = pd.DataFrame(adata.X) # It works with DF for some reason
scaled_data.columns = adata.var["approved_gene"].values
# Change the row cell names
scaled_data.index = adata.obs_names
scaled_data=scaled_data.T
sc.tl.pca(adata)

# Run Harmony and ensure the embedding has one row per cell (n_obs x n_pcs).
adata.obsm["X_pca_harmony"] = run_harmony(adata,
                                          batch_key="batch",
                                          max_iter_harmony=20,
                                          max_iter_kmeans=20,
                                          sigma=0.1,
                                          theta=2,
                                          nclust=50,
                                          verbose=True)

sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_neighbors=20, n_pcs=10)
sc.tl.umap(adata,min_dist=0.3)
sc.tl.leiden(adata, resolution=0.8)
# plot_which = ["leiden", "patient", "LTS", "timepoint", "batch"]
plot_which = ["LTS"] # patient and batch are the same in this case
sc.pl.umap(
    adata,
    color=plot_which,
    wspace=0.4
)

# %%
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

# %%
# ============================================================================
# Step 5: Plot UMAP with subset percentages for LTS vs non-LTS
# ============================================================================
celltypes = adata.obs["sctype_classification"].value_counts().index.tolist()
palette = sns.color_palette("tab20", n_colors=len(celltypes))
color_map = dict(zip(celltypes, palette))

# import random
# random.seed(42)
plot_umap_with_subset_percentages(adata, "LTS", 1, label_colormap=color_map, label_types=celltypes)
plot_umap_with_subset_percentages(adata, "LTS", 0, label_colormap=color_map, label_types=celltypes)


# ============================================================================
# Step 6: Plot violin plots of signature scores across cell types
# ============================================================================
# %%
