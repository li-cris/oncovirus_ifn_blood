import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import os
import random
import sys
import seaborn as sns
import matplotlib.pyplot as plt

from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parents[1]
import os
os.chdir(PROJECT_ROOT)

# Allow importing local packages from src/ without pip installation (temporary solution for development)
src_path = os.path.abspath("src")
if src_path not in sys.path:
    sys.path.insert(0, src_path)

from src.visualization.differential_expression import volcano_plot, ma_plot

# A count matrix of shape ‘number of samples’ x ‘number of genes’, containing read counts (non-negative integers),
# Metadata (or “column” data) of shape ‘number of samples’ x ‘number of variables’, containing sample annotations that will be used to split the data in cohorts.

# ============================================================================
# Step 1: Load AnnData and do some basic preprocessing. Keep counts as integers for DESeq2.
# ============================================================================
adata_path = "data/GSE310935_adata.h5ad"
adata = ad.read_h5ad(adata_path)

adata.obs["batch"] = adata.obs["og_id"].astype(str)
adata.layers["counts"] = adata.X.copy()

random.seed(42)
sc.pp.filter_cells(adata, min_genes=1000) # Tutorial is 10; I also ised min_genes=200
sc.pp.filter_genes(adata, min_cells=100) # Tutorial is 10; I also used min_cells=3

counts_df = pd.DataFrame(adata.X.toarray(), index=adata.obs_names, columns=adata.var_names)
metadata_df = pd.DataFrame(adata.obs.copy())

# ============================================================================
# Step 1.5: Turn into pseudobulk based on "og_id"
# ============================================================================
# REVIEW. This is just a test but next time don't sum all pseudobulks without considering sample number
adata.obs["og_id"] = adata.obs["og_id"].astype(str)
pseudobulk_counts = counts_df.groupby(adata.obs["og_id"]).sum()
pseudobulk_metadata = metadata_df.drop_duplicates(subset=["og_id"]).set_index("og_id")

# Order them so indexes match
pseudobulk_counts = pseudobulk_counts.loc[pseudobulk_metadata.index]


# ============================================================================
# Step 2: Counts remodeling wuth DeSeq2
# ============================================================================

inference = DefaultInference(n_cpus=3)
dds = DeseqDataSet(
    counts=pseudobulk_counts,  # Ensure counts are in the same order as metadata
    metadata=pseudobulk_metadata,
    design="~LTS",
    refit_cooks=True,
    inference=inference,
)
dds.deseq2() # It's an AnnData object

print(dds.var["dispersions"])

# contrast, which is a list of three strings of the form ["variable", "tested_level", "control_level"]
ds = DeseqStats(dds, contrast=["LTS", 1, 0], inference=inference)

# Wald test
ds.summary()


# ============================================================================
# Step 3: Check results and do some basic filtering for significant genes
# ============================================================================

res = ds.results_df.copy()

res["significant"] = (res["padj"] < 0.05) & (res["log2FoldChange"].abs() > 1) # Significant genes with a change of at least 2-fold (log2FC > 1 or < -1)
res["direction"] = "NS"
res.loc[(res["padj"] < 0.05) & (res["log2FoldChange"] > 1), "direction"] = "Up"
res.loc[(res["padj"] < 0.05) & (res["log2FoldChange"] < -1), "direction"] = "Down"

volcano_plot(res)
ma_plot(res)

# %%

# ============================================================================
# Step 4: Turn ENSEMBL IDs to gene symbols for enrichment analysis
# ============================================================================

# Full dictionary in "data/ens_hg_dict.tsv", col 1 is ENSEMBL ID, col2 is gene symbol. This is a bit faster than using the HGNC API for each gene
# Get gene symbols
ens_hg_dict = pd.read_csv("data/ens_hg_dict.tsv", sep="\t", header=None, names=["ensembl_id", "gene_symbol"])

# Get res["gene_symbol"] by mapping res.index (which is ENSEMBL ID) to ens_hg_dict
res["gene_symbol"] = res.index.map(ens_hg_dict.set_index("ensembl_id")["gene_symbol"])


# Get list of overexpressed gene symbols
up_genes = res.loc[res["direction"] == "Up", "gene_symbol"].dropna().tolist()
# Save as txt
with open("results/up_genes.txt", "w") as f:
    for gene in up_genes:
        f.write(gene + "\n")
# Get list of underexpressed gene symbols
down_genes = res.loc[res["direction"] == "Down", "gene_symbol"].dropna().tolist()
# Save as txt
with open("results/down_genes.txt", "w") as f:
    for gene in down_genes:
        f.write(gene + "\n")