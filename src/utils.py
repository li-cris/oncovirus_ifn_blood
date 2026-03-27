import pandas as pd
from .sctype_utils import gene_sets_prepare, sctype_score, process_cluster
from dataclasses import dataclass
@dataclass
class ScTypeOptions:
    """Class for keeping track of ScType options."""
    path_to_db_file: str = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
    cell_type: str = "Immune system" # Annotates cells from immune system
    scaled: bool =  True # Whether to use scaled data for scoring (True) or raw data (False)
    gene_names_to_uppercase: bool = True 
    cluster_key: str = 'leiden' # Type of clustering used and name in adata.obs[cluster_key]
    unknown_threshold: float = 0.25 # Threshold for setting low-confidence clusters to "Unknown" based on ncells


def check_missing_genes(gene_list, var_names, replace: bool = True):
    """
    Should work with anything as long as both inputs are lists of strings.
    """
    missing_genes = [gene for gene in gene_list if gene not in var_names]
    if missing_genes:
        print(f"Warning: Missing genes in dataset: {missing_genes}")
        if replace:
            print("Removing missing genes from gene list.")
            gene_list = [gene for gene in gene_list if gene in var_names]
    else:
        print(f"All {len(gene_list)} genes found in dataset")
    return gene_list



def run_sctype_scoring(scRNAseqData, adata, SctypeOptions: ScTypeOptions = None):
    if SctypeOptions is None:
        SctypeOptions = ScTypeOptions()
    gs_list=gene_sets_prepare(path_to_db_file=SctypeOptions.path_to_db_file,cell_type=SctypeOptions.cell_type)
    # Just in case, also update identifiers in scRNAseqData to match updates gene names
    es_max = sctype_score(scRNAseqData = scRNAseqData, scaled = SctypeOptions.scaled, gs = gs_list['gs_positive'], gs2 = gs_list['gs_negative'])

    unique_clusters = adata.obs[SctypeOptions.cluster_key].unique() # are we sure it's leiden?
    # Apply the function to each unique cluster and combine the results into a DataFrame
    cL_results = pd.concat([process_cluster(cluster,adata,es_max,SctypeOptions.cluster_key) for cluster in unique_clusters])

    # Group by cluster and select the top row based on scores
    sctype_scores = cL_results.groupby('cluster').apply(lambda x: x.nlargest(1, 'scores')).reset_index(drop=True)

    # Set low-confidence clusters to "Unknown"
    sctype_scores.loc[sctype_scores['scores'] < sctype_scores['ncells'] * SctypeOptions.unknown_threshold, 'type'] = 'Unknown'

    # Iterate over unique clusters
    adata.obs['sctype_classification'] = ""
    for cluster in sctype_scores['cluster'].unique():
        # Filter sctype_scores for the current cluster
        cl_type = sctype_scores[sctype_scores['cluster'] == cluster]
        # Get the type for the current cluster
        cl_type_value = cl_type['type'].iloc[0]
        # Update 'sctype_classification' in pbmc.obs for cells belonging to the current cluster
        adata.obs.loc[adata.obs[SctypeOptions.cluster_key] == cluster, 'sctype_classification'] = cl_type_value

    return adata