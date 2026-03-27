
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import os
from scipy.io import mmread

class AnndataBuilder:
    """
    A class to build and AnnData from 10X Genomics data, given the paths to the:
    - features file (genes)
    - barcodes file (cells)
    - matrix file (counts, or expression in general)

    Requirements:
    - maindir: the directory where the raw data is stored, which contains the feature, barcode, and matrix files for each sample
    - sample_info_path: the path to the excel file that contains the sample information, which is used to extract the sample IDs and other metadata for the obs of the AnnData object

    """
    def __init__(self, maindir, sample_info_path):
        self.maindir = maindir
        self.sample_info_path = sample_info_path

    def load_anndata_components(self, feature_path, barcode_path, matrix_path):
        """
        - features (genes)
        - barcodes (cells)
        - matrix (counts, or expression in general)
        """

        feature_path = os.path.join(self.maindir, feature_path)
        barcode_path = os.path.join(self.maindir, barcode_path)
        matrix_path = os.path.join(self.maindir, matrix_path)

        features = pd.read_csv(feature_path, sep="\t", header=None) # should have ENSEMBL, gene symbol, and gene type
        barcodes = pd.read_csv(barcode_path, sep="\t", header=None)
        # matrix is a mtx sparse matrix, but we can read it as a dense matrix for simplicity
        matrix = mmread(matrix_path)
        matrix = matrix.tocsr()  # Convert to Compressed Sparse Row format for efficient slicing

        return features, barcodes, matrix


    def build_anndata(self, features, barcodes, matrix):
        """
        Build an AnnData object from the features, barcodes, and matrix.
        - features: n_genes x n_columns (at least one column must have "ensembl_id")
        - barcodes: n_cells x 1 (column "barcode")
        """

        barcodes = barcodes.copy()
        features = features.copy()

        # Use stable identifiers so concatenation aligns genes/cells correctly. str change avoids warning message
        # Barcodes are not actually unique across samples, but we will make them unique later by adding a sample ID prefix, so we can just use the original barcodes as they are for now
        barcodes.index = barcodes["barcode"].astype(str)
        features.index = features["ensembl_id"].astype(str)

        adata = ad.AnnData(X=matrix.T, obs=barcodes, var=features)
        adata.obs_names_make_unique()
        adata.var_names_make_unique()
        
        return adata


    def get_sample_id(self, dict_df, sample_name, column_name="patient_ID"):
        """
        Extract the sample ID from the file path, which is an excel file
        sample name is the GSM portion
        """        
        sample_id = dict_df[dict_df["GSM_ID"] == sample_name][column_name].values[0]  # Get the sample ID corresponding to the sample name
        return sample_id


# def load_anndata_components(maindir, feature_path, barcode_path, matrix_path):
#     """
#     - features (genes)
#     - barcodes (cells)
#     - matrix (counts, or expression in general)
#     """
#     from scipy.io import mmread

#     feature_path = os.path.join(maindir, feature_path)
#     barcode_path = os.path.join(maindir, barcode_path)
#     matrix_path = os.path.join(maindir, matrix_path)

#     features = pd.read_csv(feature_path, sep="\t", header=None) # should have ENSEMBL, gene symbol, and gene type
#     barcodes = pd.read_csv(barcode_path, sep="\t", header=None)
#     # matrix is a mtx sparse matrix, but we can read it as a dense matrix for simplicity
#     matrix = mmread(matrix_path)
#     matrix = matrix.tocsr()  # Convert to Compressed Sparse Row format for efficient slicing

#     return features, barcodes, matrix


# def build_anndata(features, barcodes, matrix):
#     """
#     Build an AnnData object from the features, barcodes, and matrix.
#     - features: n_genes x n_columns (at least one column must have "ensembl_id")
#     - barcodes: n_cells x 1 (column "barcode")
#     """

#     barcodes = barcodes.copy()
#     features = features.copy()

#     # Use stable identifiers so concatenation aligns genes/cells correctly. str change avoids warning message
#     # Barcodes are not actually unique across samples, but we will make them unique later by adding a sample ID prefix, so we can just use the original barcodes as they are for now
#     barcodes.index = barcodes["barcode"].astype(str)
#     features.index = features["ensembl_id"].astype(str)

#     adata = ad.AnnData(X=matrix.T, obs=barcodes, var=features)
#     adata.obs_names_make_unique()
#     adata.var_names_make_unique()
    
#     return adata


# def get_sample_id(dict_df, sample_name, column_name="patient_ID"):
#     """
#     Extract the sample ID from the file path, which is an excel file
#     sample name is the GSM portion
#     """
#     sample_id = dict_df[dict_df["GSM_ID"] == sample_name][column_name].values[0]  # Get the sample ID corresponding to the sample name
#     return sample_id