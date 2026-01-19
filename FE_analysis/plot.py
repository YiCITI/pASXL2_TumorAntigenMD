import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import numpy as np
import matplotlib.patches as mpatches

def load_total_as_value(path):
    df = pd.read_csv(path)
    values = df["DELTA TOTAL"]
    return pd.DataFrame({"value": values})

df4     = load_total_as_value("4.csv")
df_prot = load_total_as_value("prot.csv")
df_no   = load_total_as_value("no.csv")
df6     = load_total_as_value("6.csv")
df10    = load_total_as_value("10.csv")

vals_4    = df4["value"].to_numpy()
vals_prot = df_prot["value"].to_numpy()
vals_no   = df_no["value"].to_numpy()
vals_6    = df6["value"].to_numpy()
vals_10   = df10["value"].to_numpy()

# Add labels
df4["Condition"]     = "Site 4 phosphorylation"
df_prot["Condition"] = "Protonated site 4 phosphorylation"
df_no["Condition"]   = "No phosphorylation"
df6["Condition"]     = "Site 6 phosphorylation"
df10["Condition"]    = "Site 10 phosphorylation"

df_all = pd.concat([df4, df_prot, df_no, df6, df10], ignore_index=True)

order = [
    "Site 4 phosphorylation",
    "No phosphorylation",                    # <-- moved up
    "Protonated site 4 phosphorylation",     # <-- moved down
    "Site 6 phosphorylation",
    "Site 10 phosphorylation",
]

short_names = [
    "Site 4",
    "None",                 # <-- aligned with new order
    "Protonated site 4",
    "Site 6",
    "Site 10",
]

color_map = {
    "Site 4 phosphorylation": "#800080",
    "Protonated site 4 phosphorylation": "#800080",
    "No phosphorylation": "#808080",
    "Site 6 phosphorylation": "#C3B1E1",
    "Site 10 phosphorylation": "#C3B1E1",
}

alphas = {
    "Site 4 phosphorylation": 1.0,
    "Protonated site 4 phosphorylation": 1.0,
    "No phosphorylation": 1.0,
    "Site 6 phosphorylation": 0.6,
    "Site 10 phosphorylation": 0.6,
}

palette_list = [color_map[c] for c in order]

TICK_FONTSIZE  = 16
LABEL_FONTSIZE = 20
MEAN_FONTSIZE  = 16

BOX_WIDTH = 0.35  # <<< 越小越“瘦”，比如 0.45 / 0.40

BASE_RIGHT_EDGE = BOX_WIDTH / 2.0  
MEAN_X_SHIFT    = 0.10             
MEAN_Y_OFFSET   = 0.10             
sns.set(style="whitegrid")
plt.figure(figsize=(9, 5))

ax = sns.boxplot(
    data=df_all,
    x="Condition",
    y="value",
    order=order,
    palette=palette_list,
    width=BOX_WIDTH,  
    linewidth=1.5,
    flierprops=dict(marker='o', markersize=5, markerfacecolor='black', alpha=0.5)
)

box_patches = [p for p in ax.patches if isinstance(p, mpatches.PathPatch)]
box_patches = box_patches[:len(order)]

for patch, cond in zip(box_patches, order):
    r, g, b, _ = patch.get_facecolor()
    patch.set_facecolor((r, g, b, alphas[cond]))
    patch.set_hatch(None)

ax.tick_params(axis='both', labelsize=TICK_FONTSIZE)
ax.set_xticklabels(short_names, fontsize=TICK_FONTSIZE)

ax.set_xlabel("")
ax.set_ylabel("Binding free energy change (kcal/mol)", fontsize=LABEL_FONTSIZE)

means = df_all.groupby("Condition")["value"].mean()
q1s   = df_all.groupby("Condition")["value"].quantile(0.25)

y_min = df_all["value"].min()
y_max = df_all["value"].max()
y_range = (y_max - y_min) if (y_max > y_min) else 1.0

for i, cond in enumerate(order):
    mean_val = float(means.loc[cond])
    q1 = float(q1s.loc[cond]) 

    x_text = i + BASE_RIGHT_EDGE + MEAN_X_SHIFT
    y_text = q1 - MEAN_Y_OFFSET * y_range
    y_text = max(y_text, y_min + 0.02 * y_range)

    ax.text(
        x_text, y_text, f"{mean_val:.2f}",
        ha="right", va="top",
        fontsize=MEAN_FONTSIZE,
        color="black",
        zorder=10,
        bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.75)
    )

def p_to_stars(p):
    if p <= 0.001: return "***"
    elif p <= 0.01: return "**"
    elif p <= 0.05: return "*"
    else: return "ns"

comparisons = [
    ("Site 4 phosphorylation", vals_4,
     "Protonated site 4 phosphorylation", vals_prot, 0, 2),

    ("Site 4 phosphorylation", vals_4,
     "No phosphorylation", vals_no, 0, 1),

    ("Protonated site 4 phosphorylation", vals_prot,
     "No phosphorylation", vals_no, 2, 1),
]

height_step = 0.08 * y_range
start_height = y_max + 0.05 * y_range

for j, (n1, v1, n2, v2, x1c, x2c) in enumerate(comparisons):
    stat, p_val = mannwhitneyu(v1, v2, alternative="two-sided")
    stars = p_to_stars(p_val)

    h = start_height + j * height_step
    ax.plot([x1c, x1c, x2c, x2c],
            [h, h + 0.01*y_range, h + 0.01*y_range, h],
            linewidth=1.2, color="black")
    ax.text((x1c + x2c) / 2, h + 0.02*y_range, stars,
            ha="center", va="bottom", fontsize=12)

ax.set_ylim(y_min, start_height + len(comparisons) * height_step + 0.05 * y_range)

plt.tight_layout()
plt.savefig("gb.svg", dpi=300, bbox_inches="tight")
plt.show()
