import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import os
import random
import sys
import seaborn as sns


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


# if not os.path.exists("results/approved_genes.csv"):
#     list_of_genes = adata.var["HUGO_symbol"].values
#     approved_genes = get_gene_symbols(list_of_genes)  # Returns DF with approved symbols
#     approved_genes.to_csv("results/approved_genes.csv", sep="\t", index=False)
# else:
#     approved_genes = pd.read_csv("results/approved_genes.csv", sep="\t")

# # Map gene names to approved symbols
# gene_name_mapping = dict(zip(approved_genes["Gene"], approved_genes["Symbol"]))
# adata.var["approved_gene"] = adata.var["HUGO_symbol"].map(gene_name_mapping).fillna(adata.var["HUGO_symbol"])


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

res["significant"] = (res["padj"] < 0.05) & (res["log2FoldChange"].abs() > 1)
res["direction"] = "NS"
res.loc[(res["padj"] < 0.05) & (res["log2FoldChange"] > 1), "direction"] = "Up"
res.loc[(res["padj"] < 0.05) & (res["log2FoldChange"] < -1), "direction"] = "Down"




import numpy as np
import matplotlib.pyplot as plt

res["minus_log10_padj"] = -np.log10(res["padj"])

plt.figure(figsize=(8,6))

for group, color in [("NS", "lightgray"), ("Up", "red"), ("Down", "blue")]:
    subset = res[res["direction"] == group]
    plt.scatter(
        subset["log2FoldChange"],
        subset["minus_log10_padj"],
        s=10,
        alpha=0.7,
        label=group,
        c=color
    )

plt.axvline(1, linestyle="--")
plt.axvline(-1, linestyle="--")
plt.axhline(-np.log10(0.05), linestyle="--")
plt.xlabel("log2 Fold Change")
plt.ylabel("-log10 adjusted p-value")
plt.legend()
plt.title("Volcano plot")
plt.show()
# 

# %%
plt.figure(figsize=(8,6))
plt.scatter(res["baseMean"], res["log2FoldChange"], s=10, alpha=0.5)
plt.xscale("log")
plt.axhline(0, linestyle="--")
plt.xlabel("Mean normalized expression")
plt.ylabel("log2 Fold Change")
plt.title("MA plot")
plt.show()

# %%
import requests

def ensembl_to_hgnc(ensembl_id, session=None, timeout=30, cache=None):
    """Map Ensembl gene ID to HGNC symbol using HGNC REST API."""
    if cache is None:
        cache = {}
    
    if ensembl_id in cache:
        return cache[ensembl_id]
    
    url = f"https://rest.genenames.org/search/ensembl_gene_id/{ensembl_id}"
    headers = {"Accept": "application/json"}
    http = session or requests
    
    response = http.get(url, headers=headers, timeout=timeout)
    response.raise_for_status()
    result = response.json()
    
    cache[ensembl_id] = result
    return result


dds.var["hgnc_symbol"] = dds.var_names.map(lambda x: ensembl_to_hgnc(x)["response"]["docs"][0]["symbol"] if ensembl_to_hgnc(x)["response"]["numFound"] > 0 else None)