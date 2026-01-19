import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import matplotlib.patches as mpatches

CSV_NAME = "rmsf.csv"
OUT_SVG = "rmsf_pos_neg_boxplot.svg"
OUT_STAT = "rmsf_statistics.csv"

TICK_FONTSIZE = 22
LABEL_FONTSIZE = 26
LEGEND_FONTSIZE = 22   

STAR_FONTSIZE = 16
LINE_WIDTH = 1.2

POS_COLOR = "#4C72B0"  
NEG_COLOR = "#DD8452"  


pos_data = {}
neg_data = {}

for d in os.listdir("."):
    if not os.path.isdir(d):
        continue

    csv_path = os.path.join(d, CSV_NAME)
    if not os.path.isfile(csv_path):
        continue

    df = pd.read_csv(csv_path)

    for _, row in df.iterrows():
        resid = int(row.iloc[0])
        rmsf = float(row.iloc[1])

        if d.startswith("pos"):
            pos_data.setdefault(resid, []).append(rmsf)
        elif d.startswith("neg"):
            neg_data.setdefault(resid, []).append(rmsf)

all_resids = sorted(set(pos_data.keys()) | set(neg_data.keys()))

stat_rows = []
pval_dict = {}

for resid in all_resids:
    pos_vals = pos_data.get(resid, [])
    neg_vals = neg_data.get(resid, [])

    if len(pos_vals) > 0 and len(neg_vals) > 0:
        _, pval = mannwhitneyu(pos_vals, neg_vals, alternative="two-sided")
    else:
        pval = None

    stat_rows.append({
        "Residue_ID": resid,
        "pos_n": len(pos_vals),
        "neg_n": len(neg_vals),
        "pos_mean": sum(pos_vals) / len(pos_vals) if pos_vals else None,
        "neg_mean": sum(neg_vals) / len(neg_vals) if neg_vals else None,
        "p_value": pval
    })

    if pval is not None:
        pval_dict[resid] = pval

pd.DataFrame(stat_rows).to_csv(OUT_STAT, index=False)

plt.figure(figsize=(1.3 * len(all_resids), 6))

plot_data = []
positions = []
group_labels = []  # "pos"/"neg"

xticks = []
xtick_labels = []

pair_xpos = {}  # resid -> {"pos": x, "neg": x}
x = 1.0

for resid in all_resids:
    pos_vals = pos_data.get(resid, [])
    neg_vals = neg_data.get(resid, [])

    this_positions = []

    if len(pos_vals) > 0:
        plot_data.append(pos_vals)
        positions.append(x)
        this_positions.append(x)
        group_labels.append("pos")
        pair_xpos.setdefault(resid, {})["pos"] = x
        x += 0.4

    if len(neg_vals) > 0:
        plot_data.append(neg_vals)
        positions.append(x)
        this_positions.append(x)
        group_labels.append("neg")
        pair_xpos.setdefault(resid, {})["neg"] = x
        x += 0.4

    if this_positions:
        xticks.append(sum(this_positions) / len(this_positions))  
        xtick_labels.append(str(resid))
        x += 0.6  

if len(plot_data) == 0:
    raise RuntimeError("No valid RMSF data found for plotting.")

box = plt.boxplot(
    plot_data,
    positions=positions,
    widths=0.3,
    patch_artist=True,
    showfliers=False
)

for patch, lab in zip(box["boxes"], group_labels):
    patch.set_facecolor(POS_COLOR if lab == "pos" else NEG_COLOR)

def p_to_star(p):
    if p is None:
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"

for resid in all_resids:
    if resid not in pair_xpos:
        continue
    if "pos" not in pair_xpos[resid] or "neg" not in pair_xpos[resid]:
        continue

    pval = pval_dict.get(resid, None)
    if pval is None:
        continue

    x1 = pair_xpos[resid]["pos"]
    x2 = pair_xpos[resid]["neg"]
    star = p_to_star(pval)

    pos_vals = pos_data.get(resid, [])
    neg_vals = neg_data.get(resid, [])
    local_max = max(max(pos_vals) if pos_vals else 0, max(neg_vals) if neg_vals else 0)

    y = local_max * 1.12 + 1e-9
    h = max(local_max * 0.03, 0.01)

    plt.plot([x1, x1, x2, x2],
             [y, y + h, y + h, y],
             lw=LINE_WIDTH, c="black")

    plt.text((x1 + x2) / 2, y + h,
             star, ha="center", va="bottom", fontsize=STAR_FONTSIZE)

plt.xticks(xticks, xtick_labels, fontsize=TICK_FONTSIZE)
plt.yticks(fontsize=TICK_FONTSIZE)
plt.tick_params(axis="both", which="major", labelsize=TICK_FONTSIZE)

plt.xlabel("Residue ID", fontsize=LABEL_FONTSIZE)
plt.ylabel("RMSF (ns)", fontsize=LABEL_FONTSIZE)

blue_patch = mpatches.Patch(color=POS_COLOR, label="positive samples")
red_patch  = mpatches.Patch(color=NEG_COLOR, label="negative samples")

plt.legend(
    handles=[blue_patch, red_patch],
    fontsize=LEGEND_FONTSIZE,
    loc="upper right",
    frameon=True,
    framealpha=0.85,
    borderpad=0.4,
    labelspacing=0.3,
    handlelength=1.2
)

plt.tight_layout()
plt.savefig(OUT_SVG, format="svg", bbox_inches="tight", pad_inches=0.25)
plt.close()

print("✅ Finished successfully.")
print(f"Figure saved to: {OUT_SVG}")
print(f"Statistics saved to: {OUT_STAT}")

