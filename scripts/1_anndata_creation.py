import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parents[1]
import os
os.chdir(PROJECT_ROOT)

from src.datatools.anndata_builder import AnndataBuilder

# Objective: Create an AnnData object for each sample, and then concatenate them into a single AnnData object for the whole dataset.

# ============================================================================
# Step 1: Set the directories and get the samples in the GEO war data folder
# ============================================================================
# Get list of GSM IDs from the directory
maindir = "data/raw/GSE310935_RAW/"
dict_path = "data/dict.xlsx"
dirlist = os.listdir(maindir)
# Getting a list of all GSM IDs in the directory
gsm_ids = set([filename.split("_")[0] for filename in dirlist if filename.endswith("_features.tsv")])
dict_df = pd.read_excel(dict_path)


# ============================================================================
# Step 2: Create an AnnData object for each sample, and concatenate them into a single AnnData object for the whole dataset
# ============================================================================
databuilder = AnndataBuilder(maindir, dict_path)
adatas = []
for id in gsm_ids:
    print(f"Processing sample {id}...")
    # Join the GSM IDs with each of the feature, barcode and matrix files to get the full paths
    # Probably just searching for the file name that matches id + something +"_features.tsv" etc. would be more robust than assuming the order of files in the directory
    feature_path = [filename for filename in dirlist if filename.startswith(id) and filename.endswith("_features.tsv")][0]
    barcode_path = [filename for filename in dirlist if filename.startswith(id) and filename.endswith("_barcodes.tsv")][0]
    matrix_path = [filename for filename in dirlist if filename.startswith(id) and filename.endswith("_matrix.mtx")][0]

    # feature_path = "data/raw/GSE310935_RAW/GSM9312617_SC299_9_features.tsv"
    # barcode_path = "data/raw/GSE310935_RAW/GSM9312617_SC299_9_barcodes.tsv"
    # matrix_path = "data/raw/GSE310935_RAW/GSM9312617_SC299_9_matrix.mtx"

    features, barcodes, matrix = databuilder.load_anndata_components(feature_path, barcode_path, matrix_path)
    # features: n_genes x number (depending on what columns we have of gene names)
    # barcodes: n_cells x 1 (cell barcodes)
    # matrix: n_genes x n_cells sparse matrix
    # adata: n_cells x n_genes AnnData object with X as the matrix, obs as barcodes, and var as features

    # Building a single anndata object for one sample, we can then concatenate them later
    barcodes.columns = ["barcode"]
    features.columns = ["ensembl_id", "HUGO_symbol", "gene_type"]
    adata = databuilder.build_anndata(features, barcodes, matrix)

    # This part is more specific to each dataset, since the obs you want depend on what's available
    adata.obs["patient"] = databuilder.get_sample_id(dict_df, id, column_name="patient_ID")
    adata.obs["LTS"] = databuilder.get_sample_id(dict_df, id, column_name="LTS")
    adata.obs["timepoint"] = databuilder.get_sample_id(dict_df, id, column_name="point")
    adata.obs["time"] = databuilder.get_sample_id(dict_df, id, column_name="time")
    adata.obs["og_id"] = databuilder.get_sample_id(dict_df, id, column_name="full_ID")

    # Join adata objects that we keep getting by .var, since what we want to increase is .obs (cells)
    # But keep .var in adata_all
    # Batch categories is patient, so no need to add a batch column, we can just use the patient column as the batch key when concatenating
    adatas.append(adata)

adata_all = ad.concat(adatas, join="outer", merge="first", index_unique="-")  # Concatenate along the obs axis (cells), keeping all genes (outer join) and merging var annotations by taking the first non-null value
adata_all.obs.index.name = "identifier" # Make obs names unique by adding a suffix if there are duplicates, which there will be since barcodes are not unique across samples

# adata_all 
# AnnData object with n_obs × n_vars = 34482 × 36601
#     obs: 'barcode', 'patient', 'LTS', 'timepoint', 'time', 'og_id', 'batch'
#     var: 'ensembl_id', 'HUGO_symbol', 'gene_type'

# ============================================================================
# Step 3: Save the AnnData object as a .h5ad file for later use in the analysis pipeline
# ============================================================================
# adata_all.write_h5ad("data/GSE310935_adata.h5ad")

