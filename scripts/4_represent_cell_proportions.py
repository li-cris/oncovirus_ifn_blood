import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from pathlib import Path
PROJECT_ROOT = Path(__file__).resolve().parents[1]
import os
os.chdir(PROJECT_ROOT)

prop_path = "results/umap_percentages.xlsx"

df = pd.read_excel(prop_path, index_col=0)
df = df.iloc[:, :4]
# Fill missing cell types (important!)
df = df.fillna(0)

# --- 1. Heatmap ---
plt.figure(figsize=(10, 8))
sns.heatmap(df, cmap="BuPu", linewidths=0.5)
plt.title("Immune Cell Proportions")
plt.ylabel("Cell Type")
plt.xlabel("Condition")
plt.tight_layout()
plt.show()

# --- 2. Delta heatmap ---
new_df = df[(df["STS_b"] != 0) & (df["STS_f"] != 0) & (df["LTS_b"] != 0) & (df["LTS_f"] != 0)]

delta = pd.DataFrame({
    "STS_change": new_df["STS_f"] - new_df["STS_b"],
    "LTS_change": new_df["LTS_f"] - new_df["LTS_b"]
})

plt.figure(figsize=(6, 8))
sns.heatmap(delta, cmap="coolwarm", center=0, linewidths=0.5)
plt.title("Change in Cell Proportions (After - Before)")
plt.tight_layout()
plt.show()


# --- 3. Delta heatmap (within-time) ---
# filter df to keep only rows wihtout zeros
# new_df = df[(df["STS_b"] != 0) & (df["STS_f"] != 0) & (df["LTS_b"] != 0) & (df["LTS_f"] != 0)]
# delta = pd.DataFrame({
#     "baseline": new_df["LTS_b"] - new_df["STS_f"],
#     "follow-up": new_df["LTS_f"] - new_df["STS_f"]
# })

delta = pd.DataFrame({
    "baseline": df["LTS_b"] - df["STS_f"],
    "follow-up": df["LTS_f"] - df["STS_f"]
})

plt.figure(figsize=(6, 8))
sns.heatmap(delta, cmap="coolwarm", center=0, linewidths=0.5)
plt.title("Change in Cell Proportions (LTS - STS)")
plt.tight_layout()
plt.show()
