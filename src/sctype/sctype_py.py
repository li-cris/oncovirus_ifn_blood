import numpy as np
import pandas as pd
from sklearn.preprocessing import scale, MinMaxScaler
import requests
import scanpy as sc
from collections import defaultdict
import concurrent.futures
import multiprocessing
from functools import partial
from concurrent.futures import ThreadPoolExecutor





def gene_sets_prepare(path_to_db_file, cell_type):
    # Read data from Excel file
    cell_markers = pd.read_excel(path_to_db_file)
    # Filter by cell_type
    cell_markers = cell_markers[cell_markers['tissueType'] == cell_type]
    # Preprocess geneSymbolmore1 and geneSymbolmore2 columns
    #for col in ['geneSymbolmore1', 'geneSymbolmore2']:
    #    cell_markers[col] = cell_markers[col].str.replace(" ", "").str.upper()
    for col in ['geneSymbolmore1', 'geneSymbolmore2']:
        if col in cell_markers.columns:
            cell_markers[col] = cell_markers[col].fillna('').str.replace(" ", "").str.upper()
    # Stack and drop duplicates to get unique gene names
    gene_names = pd.concat([cell_markers['geneSymbolmore1'], cell_markers['geneSymbolmore2']]).str.split(',', expand=True).stack().drop_duplicates().reset_index(drop=True)
    gene_names = gene_names[gene_names.str.strip() != '']
    gene_names = gene_names[gene_names != 'None'].unique()
    # Get approved symbols for gene names
    res = get_gene_symbols(set(gene_names))
    res = dict(zip(res['Gene'], res['Symbol']))
    # Process gene symbols
    for col in ['geneSymbolmore1', 'geneSymbolmore2']:
        cell_markers[col] = cell_markers[col].apply(lambda row: process_gene_symbols(row, res)).str.replace("///", ",").str.replace(" ", "")
    # Group by cellName and create dictionaries of gene sets
    gs_positive = cell_markers.groupby('cellName')['geneSymbolmore1'].apply(lambda x: list(set(','.join(x).split(',')))).to_dict()
    gs_negative = cell_markers.groupby('cellName')['geneSymbolmore2'].apply(lambda x: list(set(','.join(x).split(',')))).to_dict()
    return {'gs_positive': gs_positive, 'gs_negative': gs_negative}


def process_gene_symbols(gene_symbols, res):
    if pd.isnull(gene_symbols):
        return ""
    markers_all = gene_symbols.upper().split(',')
    markers_all = [marker.strip().upper() for marker in markers_all if marker.strip().upper() not in ['NA', '']]
    markers_all = sorted(markers_all)
    if len(markers_all) > 0:
        markers_all = [res.get(marker) for marker in markers_all]
        markers_all = [symbol for symbol in markers_all if symbol is not None]
        markers_all = list(set(markers_all))
        return ','.join(markers_all)
    else:
        return ""


def hgnc_lookup(term, session=None, timeout=30, cache=None):
    """Look up gene symbol from HGNC REST API with optional caching."""
    if cache is None:
        cache = {}
    
    if term in cache:
        return cache[term]
    
    url = f"https://rest.genenames.org/search/{term}"
    headers = {"Accept": "application/json"}
    http = session or requests
    response = http.get(url, headers=headers, timeout=timeout)
    response.raise_for_status()
    result = response.json()
    cache[term] = result
    return result


def get_gene_symbols(genes):
    """Resolve gene symbols via HGNC, printing unresolved genes."""
    data = {"Gene": [], "Symbol": []}
    unresolved = []
    cache = {}
    
    with requests.Session() as session:
        for gene in genes:
            try:
                payload = hgnc_lookup(gene, session=session, cache=cache)
                docs = payload.get("response", {}).get("docs", [])
                symbol = docs[0].get("symbol") if docs else None
                if symbol:
                    data["Gene"].append(gene)
                    data["Symbol"].append(symbol)
                else:
                    unresolved.append(gene)
            except requests.RequestException as e:
                # Track API errors as unresolved.
                unresolved.append(gene)
    
    if unresolved:
        print(f"Genes not found in HGNC: {', '.join(sorted(unresolved))}")
    
    return pd.DataFrame(data)



def sctype_score(scRNAseqData, scaled=True, gs=None, gs2=None, gene_names_to_uppercase=True, *args, **kwargs):
    marker_stat = defaultdict(int, {gene: sum(gene in genes for genes in gs.values()) for gene in set(gene for genes in gs.values() for gene in genes)})
    marker_sensitivity = pd.DataFrame({'gene_': list(marker_stat.keys()), 'score_marker_sensitivity': list(marker_stat.values())})
    # Rescaling the score_marker_sensitivity column
    # grab minimum and maximum
    min_value=1
    max_value= len(gs)
    # Apply the formula to the column
    marker_sensitivity['score_marker_sensitivity'] = 1 - (marker_sensitivity['score_marker_sensitivity'] - min_value) / (max_value - min_value)
    # Convert gene names to Uppercase
    if gene_names_to_uppercase:
        scRNAseqData.index = scRNAseqData.index.str.upper()
    # Subselect genes only found in data
    names_gs_cp = list(gs.keys())
    names_gs_2_cp = list(gs2.keys())
    gs_ = {key: [gene for gene in scRNAseqData.index if gene in gs[key]] for key in gs}
    gs2_ = {key: [gene for gene in scRNAseqData.index if gene in gs2[key]] for key in gs2}
    gs__ = dict(zip(names_gs_cp, gs_.values()))
    gs2__ = dict(zip(names_gs_2_cp, gs2_.values()))
    cell_markers_genes_score = marker_sensitivity[marker_sensitivity['gene_'].isin(set.union(*map(set, gs__.values())))]
    # Z-scale if not
    if not scaled:
        Z = scale(scRNAseqData.T).T
    else:
        Z = scRNAseqData
    # Multiply by marker sensitivity
    for _, row in cell_markers_genes_score.iterrows():
        Z.loc[row['gene_']] *= row['score_marker_sensitivity']
    marker_genes=list(set().union(*gs__.values(), *gs2__.values()))
    Z = Z.loc[marker_genes]
    # Combine scores
    es = pd.DataFrame(index=gs__.keys(),columns=Z.columns,data=np.zeros((len(gs__), Z.shape[1])))
    for gss_, genes in gs__.items():
        for j in range(Z.shape[1]):
            gs_z = Z.loc[genes, Z.columns[j]]
            gz_2 = Z.loc[gs2__[gss_], Z.columns[j]] * -1 if gs2__ and gss_ in gs2__ else pd.Series(dtype=np.float64)
            sum_t1 = np.sum(gs_z) / np.sqrt(len(gs_z))
            sum_t2 = np.sum(gz_2) / np.sqrt(len(gz_2)) if not gz_2.empty else 0
            if pd.isna(sum_t2):
                sum_t2 = 0
            es.loc[gss_, Z.columns[j]] = sum_t1 + sum_t2
    es = es.dropna(how='all')
    return es

def process_cluster(cluster,adata,es_max,clustering):
    cluster_data = es_max.loc[:, adata.obs.index[adata.obs[clustering] == cluster]]
    es_max_cl = cluster_data.sum(axis=1).sort_values(ascending=False)
    top_scores = es_max_cl.head(10)
    ncells = sum(adata.obs[clustering] == cluster)
    return pd.DataFrame({
        'cluster': cluster,
        'type': top_scores.index,
        'scores': top_scores.values,
        'ncells': ncells
    })
