import anndata as ad
import numpy as np
import pandas as pd
import os

os.chdir("/Users/cristinali/Documents/Programming/PhD/onco_ifn")
# add src to temporary path to import the functions for gene name conversion
import sys
sys.path.insert(0, ".") # temporary, should do something else when doing this
from src.sctype.sctype_py import get_gene_symbols

adata_path = "data/GSE310935_adata.h5ad"
adata = ad.read_h5ad(adata_path)

list_of_genes = adata.var["HUGO_symbol"].values
approved_genes = get_gene_symbols(list_of_genes) # Gives a DF, with "Symbol" being the approved name in the 2nd column
approved_genes.to_csv("results/approved_genes.csv", sep = "\t", index=False)