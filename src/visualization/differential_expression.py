import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def volcano_plot(res, xlabel="log2 Fold Change", ylabel="-log10 adjusted p-value", title="Volcano plot"):

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
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    plt.title(title)
    plt.show()



def ma_plot(res, xlabel="Mean normalized expression", ylabel="log2 Fold Change", title="MA plot"):
    plt.figure(figsize=(8,6))
    plt.scatter(res["baseMean"], res["log2FoldChange"], s=10, alpha=0.5)
    plt.xscale("log")
    plt.axhline(0, linestyle="--")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.show()