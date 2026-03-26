

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
