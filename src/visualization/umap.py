import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

def plot_umap_with_subset_percentages(
    adata,
    filter_col,
    filter_value,
    label_colormap,
    label_types,
    label_col="sctype_classification",
):
    adata_sub = adata[adata.obs[filter_col] == filter_value].copy()

    # Keep a fixed category order for the original label column
    adata_sub.obs[label_col] = pd.Categorical(
        adata_sub.obs[label_col],
        categories=label_types,
        ordered=True
    )

    # Percentages within the subset
    perc = (
        adata_sub.obs[label_col]
        .value_counts(normalize=True)
        .mul(100)
        .round(1)
    )

    # Build legend labels in the global order, even if some are absent
    new_labels = {
        ct: f"{ct} ({perc.get(ct, 0)}%)"
        for ct in label_types
    }

    plot_col = f"{label_col}_legend"

    # Create plotting column with fixed order
    adata_sub.obs[plot_col] = pd.Categorical(
        adata_sub.obs[label_col].map(new_labels),
        categories=[new_labels[ct] for ct in label_types],
        ordered=True
    )

    # Assign colors to the plotted column
    adata_sub.uns[f"{plot_col}_colors"] = [
        label_colormap[ct] for ct in label_types
    ]

    sc.pl.umap(
        adata_sub,
        color=plot_col,
        title=f"UMAP for {filter_col} = {filter_value}",
        legend_loc="right margin"
    )