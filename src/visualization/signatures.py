import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import seaborn as sns
import pandas as pd
import numpy as np

def p_to_stars(p):
    if p < 1e-4:
        return "****"
    elif p < 1e-3:
        return "***"
    elif p < 1e-2:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"


def reorder_stats_significance(df, celltype_order):
    """
    df: DataFrame with columns "cell_type", "p-value"
    """
    if "cell_type" not in df.columns or "p-value" not in df.columns:
        raise ValueError("DataFrame must contain 'cell_type' and 'p-value' columns")

    df["significance"] = df["p-value"].apply(p_to_stars)
    df["p_label"] = df["p-value"].apply(lambda p: f"p={p:.2e}")
    df = df.set_index("cell_type").loc[celltype_order].reset_index()

    return df


def plot_score_on_violin(adata, celltype_order,
                                hue_order: list = None,
                                my_palette = "Set2",
                                plot_sig: bool = True,
                                stats_df = None):

    fig, ax = plt.subplots(figsize=(12, 6))

    if hue_order is None:
        hue_order = sorted(adata.obs["OS_label"].unique())

    sns.violinplot(
        x="sctype_classification",
        y="signature_score",
        hue="OS_label",
        hue_order=hue_order,
        order=celltype_order,
        data=adata.obs,
        palette=my_palette,
        split=False,
        inner="quartile",
        ax=ax
    )

    # Put grid/line behind violins
    ax.set_axisbelow(True)

    # Faint horizontal grid every 0.5
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
    ax.grid(which="major", linestyle="--", linewidth=0.6, alpha=0.25, color="gray", zorder=0)

    # y=0 reference line behind violins
    n_types = adata.obs["sctype_classification"].nunique()
    ax.hlines(
        y=0,
        xmin=-0.5,
        xmax=n_types - 0.5,
        colors="lightgray",
        linestyles="dashed",
        linewidth=1.0,
        zorder=0
    )

    # Significance stars
    if plot_sig:
        for i, row in enumerate(stats_df.itertuples()):
            celltype = row._1  # or row[1] depending on your df
            label = row.significance  # or row.p_label

            # Get max y-value for this cell type
            subset = adata.obs[adata.obs["sctype_classification"] == celltype]
            y_max = subset["signature_score"].max()

            ax.text(
                i,
                y_max + 0.2,  # adjust spacing
                label,
                ha="center",
                va="bottom",
                fontsize=10,
                color="black"
            )

    # Remove extra corner spacing
    ax.set_xlim(-0.5, n_types - 0.5)
    ax.set_ylim(adata.obs["signature_score"].min() - 0.5, adata.obs["signature_score"].max() + 1)
    ax.margins(x=0)

    ax.set_title("Distribution of Signature Scores by Cell Type and LTS Status")
    ax.set_xlabel("Cell Type")
    ax.set_ylabel("Signature Score")
    plt.xticks(rotation=45)
    plt.legend(title="Survival group", loc="upper right")
    plt.tight_layout()
    plt.show()