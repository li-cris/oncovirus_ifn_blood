import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp
import os
import random
import sys
import seaborn as sns
import matplotlib.pyplot as plt

from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parents[1]
import os
os.chdir(PROJECT_ROOT)

# Allow importing local packages from src/ without pip installation (temporary solution for development)
src_path = os.path.abspath("src")
if src_path not in sys.path:
    sys.path.insert(0, src_path)


# %%

# ============================================================================
# Step 1: Harllmark gene setss
# ============================================================================

from gseapy import Msigdb

msig = Msigdb()
gmt = msig.get_gmt(category="h.all", dbver="2026.1.Hs")

# msig.list_dbver()
# msig.list_category(dbver="2026.1.Hs")

# ============================================================================
# Step 2
# ============================================================================

names = gp.get_library_name(organism="Human")
go_bp = gp.get_library(name='GO_Biological_Process_2025', organism='Human')
# go_cc = gp.get_library(name='GO_Cellular_Component_2025', organism='Human')
# go_mf = gp.get_library(name='GO_Molecular_Function_2025', organism='Human')



# ============================================================================
# Step 3.1: Over representation analysis. Gene list is enough. UP
# ============================================================================
# up genes and down genes from 6

up = up_genes.copy()
down = down_genes.copy()


# if you are only intrested in dataframe that enrichr returned, please set outdir=None
enr = gp.enrichr(gene_list=up, # or "./tests/data/gene_list.txt",
                 gene_sets=['GO_Biological_Process_2025','KEGG_2021_Human'],
                 organism='human', # don't forget to set organism to the one you desired! e.g. Yeast
                 outdir=None, # don't write to disk
                )

enrdf = pd.DataFrame(enr.results[enr.results["Adjusted P-value"] < 0.05])

from gseapy import barplot, dotplot
ax = dotplot(enrdf,
              column="Adjusted P-value",
              x='Gene_set', # set x axis, so you could do a multi-sample/library comparsion
              size=10,
              top_term=10,
              figsize=(3,5),
              title = "Terms enriched in upregulated genes",
              xticklabels_rot=45, # rotate xtick labels
              show_ring=True, # set to False to revmove outer ring
              marker='o',
             )


# ============================================================================
# Step 3.2: Over representation analysis. Gene list is enough. DOWN
# ============================================================================
down = down_genes.copy()

enr = gp.enrichr(gene_list=down, # or "./tests/data/gene_list.txt",
                 gene_sets=['GO_Biological_Process_2025','KEGG_2021_Human'],
                 organism='human', # don't forget to set organism to the one you desired! e.g. Yeast
                 outdir=None, # don't write to disk
                )

enrdf = pd.DataFrame(enr.results[enr.results["Adjusted P-value"] < 0.05])

from gseapy import barplot, dotplot
ax = dotplot(enrdf,
              column="Adjusted P-value",
              x='Gene_set', # set x axis, so you could do a multi-sample/library comparsion
              size=10,
              top_term=10,
              figsize=(3,5),
              title = "Terms enriched in downregulated genes",
              xticklabels_rot=45, # rotate xtick labels
              show_ring=True, # set to False to revmove outer ring
              marker='o',
             )

