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

# Repeat until section 4 in 2_anndata_integration.py
# Get list of GSM IDs from the directory
adata_path = "data/GSE310935_adata.h5ad"
adata = ad.read_h5ad(adata_path)
# Gene list: "CD69", "CD44", "ICOS", "CD40LG", "TNF", "IFNG", "GZMB", "PRF1", "IL2", "FOS", "JUN", "CXCR4", "BATF"


gene_list = ["CD69", "CD44", "ICOS", "CD40LG", "TNF", "IFNG", "GZMB", "PRF1", "IL2", "FOS", "JUN", "CXCR4", "BATF"]
# Make sure all genes are in adata.var["approved_gene"]
missing_genes = [gene for gene in gene_list if gene not in adata.var["approved_gene"].values]


