from pathlib import Path
from io import StringIO
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

try:
    from scipy import stats
    SCIPY_AVAILABLE = True
except Exception:
    SCIPY_AVAILABLE = False


# ============================================================
# 1. Basic settings
# ============================================================

# Only analyze these three systems.
FOLDERS = ["phos", "prot", "wt"]
N_REPLICATES = 5

# Colors consistent with previous scripts:
# phos = purple, prot = blue, wt = gray.
GROUP_COLORS = {
    "phos": "#800080",
    "prot": "#0072B2",
    "wt":   "#666666",
}

# Five parallel trajectories are distinguished by alpha.
REPLICA_ALPHAS = [0.30, 0.42, 0.54, 0.66, 0.80]

BOX_WIDTH = 0.18
WITHIN_GROUP_SPACING = 0.23
GROUP_SPACING = 1.45

LABEL_FONTSIZE = 32
TICK_FONTSIZE  = 24
LEGEND_FONTSIZE = 18

# Mean line style: mean of the five trajectory means for each system.
MEAN_LINE_STYLE = "--"
MEAN_LINE_WIDTH = 2.0
MEAN_LINE_ALPHA = 0.95
MEAN_TEXT_FONTSIZE = 14

# Significance settings.
# Statistical unit: one independent trajectory mean.
# Test: two-sided Welch's t-test, n = 5 per system.
# Multiple-comparison correction: Holm correction within each plotted energy term.
SIGNIFICANCE_TEST = "welch_ttest"
P_VALUE_CORRECTION = "holm"

# Plot Holm-adjusted p-values on the figure.
# The CSV keeps both raw and adjusted p-values.
P_VALUE_TO_PLOT = "p_adj"

# If True, only p_adj < 0.05 comparisons are labeled.
# If False, all three pairwise p-values are labeled.
PLOT_ONLY_SIGNIFICANT = False

P_VALUE_FONTSIZE = 13
SIG_LINE_WIDTH = 1.1
SIG_LINE_GAP_FRACTION = 0.045
SIG_TEXT_GAP_FRACTION = 0.010
Y_TOP_EXTRA_FRACTION = 0.28


# ============================================================
# 2. Energy terms
# ============================================================

# - DELTA uses the DELTA TOTAL column in DELTA Energy Terms
# - complex/receptor/ligand use TOTAL in their corresponding sections
ENERGY_TERMS = {
    "delta_total": {
        "section": "DELTA Energy Terms",
        "column": "DELTA TOTAL",
        "ylabel": r"$\Delta G_{\mathrm{total}}$  (kcal/mol)",
        "output": "delta_total_boxplot.svg",
        "print_name": "binding_delta_total",
    },
    "ligand": {
        "section": "Ligand Energy Terms",
        "column": "TOTAL",
        "ylabel": r"$G_{\mathrm{ligand}}$  (kcal/mol)",
        "output": "ligand_free_energy_boxplot.svg",
        "print_name": "ligand_total",
    },
    "receptor": {
        "section": "Receptor Energy Terms",
        "column": "TOTAL",
        "ylabel": r"$G_{\mathrm{receptor}}$  (kcal/mol)",
        "output": "receptor_free_energy_boxplot.svg",
        "print_name": "receptor_total",
    },
    "complex": {
        "section": "Complex Energy Terms",
        "column": "TOTAL",
        "ylabel": r"$G_{\mathrm{complex}}$  (kcal/mol)",
        "output": "complex_free_energy_boxplot.svg",
        "print_name": "complex_total",
    },
}


# ============================================================
# 3. IO helpers
# ============================================================

def read_energy_column_from_mmpbsa_csv(csv_file, section_name, value_column):
    """
    Read one energy column from an AMBER MMPBSA.py csv output.

    If the same section appears more than once, the last section is used,
    following the original script logic.
    """
    csv_file = Path(csv_file)

    with open(csv_file, "r", encoding="utf-8", errors="ignore") as f:
        lines = f.readlines()

    start_indices = [
        i for i, line in enumerate(lines)
        if line.strip() == section_name
    ]

    if not start_indices:
        raise ValueError(f"Cannot find '{section_name}' in {csv_file}")

    start = start_indices[-1]
    header_idx = start + 1
    data_start = start + 2

    if header_idx >= len(lines):
        raise ValueError(f"Missing header after '{section_name}' in {csv_file}")

    section_lines = [lines[header_idx]]

    for line in lines[data_start:]:
        if line.strip() == "":
            break
        if "," not in line:
            break
        section_lines.append(line)

    section_text = "".join(section_lines)
    df = pd.read_csv(StringIO(section_text), skipinitialspace=True)
    df.columns = [c.strip() for c in df.columns]

    if value_column not in df.columns:
        raise ValueError(
            f"'{value_column}' column not found in section '{section_name}' of {csv_file}. "
            f"Columns are: {df.columns.tolist()}"
        )

    return df[value_column].dropna().to_numpy(dtype=float)


def get_replica_offsets():
    """
    Symmetric offsets for five trajectory boxplots within one system.
    """
    center = (N_REPLICATES - 1) / 2.0
    return np.array(
        [(i - center) * WITHIN_GROUP_SPACING for i in range(N_REPLICATES)],
        dtype=float,
    )


def collect_energy_data(term_config):
    """
    Collect frame-level values for each trajectory boxplot and trajectory-level
    means for statistics.
    """
    all_data = []
    positions = []
    box_colors = []
    box_alphas = []
    group_centers = []
    group_labels = []
    summary_rows = []
    group_mean_rows = []
    trajectory_means_by_group = {}

    offsets = get_replica_offsets()

    for g_idx, folder in enumerate(FOLDERS):
        folder_path = Path(folder)

        if not folder_path.is_dir():
            raise FileNotFoundError(f"Folder not found: {folder}")

        csv_files = sorted(folder_path.glob("*.csv"))

        if len(csv_files) != N_REPLICATES:
            raise ValueError(
                f"{folder} should contain {N_REPLICATES} csv files, "
                f"but found {len(csv_files)}: {[f.name for f in csv_files]}"
            )

        group_center = g_idx * GROUP_SPACING
        group_centers.append(group_center)
        group_labels.append(folder)

        trajectory_means = []

        for r_idx, csv_file in enumerate(csv_files):
            values = read_energy_column_from_mmpbsa_csv(
                csv_file=csv_file,
                section_name=term_config["section"],
                value_column=term_config["column"],
            )

            all_data.append(values)
            positions.append(group_center + offsets[r_idx])
            box_colors.append(GROUP_COLORS[folder])
            box_alphas.append(REPLICA_ALPHAS[r_idx])

            mean_value = np.mean(values)
            std_value = np.std(values, ddof=1) if len(values) > 1 else np.nan
            trajectory_means.append(mean_value)

            summary_rows.append({
                "energy_type": term_config["print_name"],
                "folder": folder,
                "replica": r_idx + 1,
                "csv_file": csv_file.name,
                "n_frames": len(values),
                "trajectory_mean": mean_value,
                "trajectory_std": std_value,
            })

            print(
                f"{term_config['print_name']} | "
                f"{folder} | replica {r_idx + 1} | {csv_file.name} | "
                f"n = {len(values)} | "
                f"mean = {mean_value:.3f} | "
                f"std = {std_value:.3f}"
            )

        trajectory_means = np.asarray(trajectory_means, dtype=float)
        trajectory_means_by_group[folder] = trajectory_means

        group_mean_rows.append({
            "energy_type": term_config["print_name"],
            "folder": folder,
            "n_replicates": len(trajectory_means),
            "mean_of_trajectory_means": float(np.mean(trajectory_means)),
            "sd_of_trajectory_means": float(np.std(trajectory_means, ddof=1)),
            "sem_of_trajectory_means": float(
                np.std(trajectory_means, ddof=1) / np.sqrt(len(trajectory_means))
            ),
        })

    return {
        "all_data": all_data,
        "positions": positions,
        "box_colors": box_colors,
        "box_alphas": box_alphas,
        "group_centers": group_centers,
        "group_labels": group_labels,
        "summary_rows": summary_rows,
        "group_mean_rows": group_mean_rows,
        "trajectory_means_by_group": trajectory_means_by_group,
    }


# ============================================================
# 4. Statistics
# ============================================================

def welch_ttest_pvalue(x, y):
    """
    Two-sided Welch's t-test p-value.
    Falls back to a small exact permutation test if scipy is unavailable.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    if SCIPY_AVAILABLE:
        return float(stats.ttest_ind(x, y, equal_var=False, alternative="two-sided").pvalue)

    # Exact permutation fallback for n=5 vs n=5.
    pooled = np.concatenate([x, y])
    n_x = len(x)
    observed = abs(x.mean() - y.mean())

    more_extreme = 0
    total = 0

    from itertools import combinations

    for idx_x in combinations(range(len(pooled)), n_x):
        mask = np.zeros(len(pooled), dtype=bool)
        mask[list(idx_x)] = True

        x_perm = pooled[mask]
        y_perm = pooled[~mask]
        stat = abs(x_perm.mean() - y_perm.mean())

        if stat >= observed - 1e-15:
            more_extreme += 1

        total += 1

    return float(more_extreme / total)


def holm_adjust(p_values):
    """
    Holm step-down correction.
    """
    p_values = np.asarray(p_values, dtype=float)
    m = len(p_values)

    order = np.argsort(p_values)
    adjusted = np.empty(m, dtype=float)

    running_max = 0.0

    for rank, idx in enumerate(order):
        adj = (m - rank) * p_values[idx]
        running_max = max(running_max, adj)
        adjusted[idx] = min(running_max, 1.0)

    return adjusted


def compute_pairwise_pvalues(term_config, trajectory_means_by_group):
    """
    Pairwise system comparisons using trajectory means as independent samples.
    """
    rows = []

    for folder_a, folder_b in combinations(FOLDERS, 2):
        x = trajectory_means_by_group[folder_a]
        y = trajectory_means_by_group[folder_b]

        p_raw = welch_ttest_pvalue(x, y)

        rows.append({
            "energy_type": term_config["print_name"],
            "group_1": folder_a,
            "group_2": folder_b,
            "n_1": len(x),
            "n_2": len(y),
            "mean_1": float(np.mean(x)),
            "mean_2": float(np.mean(y)),
            "difference_mean_1_minus_mean_2": float(np.mean(x) - np.mean(y)),
            "test": "two-sided Welch's t-test on trajectory means",
            "p_raw": p_raw,
        })

    stats_df = pd.DataFrame(rows)

    if P_VALUE_CORRECTION.lower() == "holm":
        stats_df["p_adj"] = holm_adjust(stats_df["p_raw"].to_numpy())
        stats_df["p_adjust_method"] = "Holm"
    else:
        stats_df["p_adj"] = stats_df["p_raw"]
        stats_df["p_adjust_method"] = "none"

    return stats_df


def format_pvalue(p):
    if not np.isfinite(p):
        return "p=NA"
    if p < 1e-4:
        return "p<1e-4"
    if p < 1e-3:
        return f"p={p:.1e}"
    return f"p={p:.3f}"


# ============================================================
# 5. Plotting
# ============================================================

def add_pairwise_pvalue_bar(ax, x1, x2, y, text, line_height):
    ax.plot(
        [x1, x1, x2, x2],
        [y, y + line_height, y + line_height, y],
        color="black",
        linewidth=SIG_LINE_WIDTH,
        clip_on=False,
    )

    ax.text(
        (x1 + x2) / 2.0,
        y + line_height,
        text,
        ha="center",
        va="bottom",
        fontsize=P_VALUE_FONTSIZE,
    )


def plot_boxplot(
    all_data,
    positions,
    box_colors,
    box_alphas,
    group_centers,
    group_labels,
    trajectory_means_by_group,
    pvalue_df,
    ylabel,
    output_file,
):
    """
    Plot individual trajectory boxplots, system-level mean dashed lines,
    and pairwise p-values.
    """
    fig, ax = plt.subplots(figsize=(14, 6))

    bp = ax.boxplot(
        all_data,
        positions=positions,
        widths=BOX_WIDTH,
        patch_artist=True,
        showfliers=False,
        whis=1.5,
        manage_ticks=False,
    )

    for patch, color, alpha in zip(bp["boxes"], box_colors, box_alphas):
        patch.set_facecolor(color)
        patch.set_alpha(alpha)
        patch.set_edgecolor("black")
        patch.set_linewidth(1.4)

    for key in ["whiskers", "caps", "medians"]:
        for item in bp[key]:
            item.set_color("black")
            item.set_linewidth(1.3)

    ax.set_xticks(group_centers)
    ax.set_xticklabels(group_labels, fontsize=TICK_FONTSIZE)

    ax.tick_params(axis="y", labelsize=TICK_FONTSIZE)
    ax.tick_params(axis="x", labelsize=TICK_FONTSIZE)

    ax.set_ylabel(ylabel, fontsize=LABEL_FONTSIZE)

    # Dividers between systems.
    for i in range(len(group_centers) - 1):
        x_sep = (group_centers[i] + group_centers[i + 1]) / 2
        ax.axvline(x_sep, color="gray", linewidth=0.8, alpha=0.3)

    # Current y-range before adding annotations.
    y_data_min = min(np.min(values) for values in all_data)
    y_data_max = max(np.max(values) for values in all_data)
    y_range = max(y_data_max - y_data_min, 1e-6)

    # Dashed mean lines: mean of five trajectory means.
    offsets = get_replica_offsets()
    line_half_width = max(abs(offsets[0]), abs(offsets[-1])) + BOX_WIDTH * 0.65

    for center, folder in zip(group_centers, group_labels):
        group_mean = float(np.mean(trajectory_means_by_group[folder]))

        ax.hlines(
            group_mean,
            center - line_half_width,
            center + line_half_width,
            colors=GROUP_COLORS[folder],
            linestyles=MEAN_LINE_STYLE,
            linewidth=MEAN_LINE_WIDTH,
            alpha=MEAN_LINE_ALPHA,
        )

        text_y = group_mean + 0.018 * y_range
        ax.text(
            center,
            text_y,
            f"mean={group_mean:.2f}",
            ha="center",
            va="bottom",
            fontsize=MEAN_TEXT_FONTSIZE,
            color=GROUP_COLORS[folder],
        )

    # Pairwise p-value annotations.
    sig_rows = pvalue_df.copy()

    if PLOT_ONLY_SIGNIFICANT:
        sig_rows = sig_rows[sig_rows["p_adj"] < 0.05].copy()

    if not sig_rows.empty:
        center_map = {folder: center for folder, center in zip(group_labels, group_centers)}

        sig_rows["_plot_p"] = sig_rows[P_VALUE_TO_PLOT]
        sig_rows = sig_rows.sort_values("_plot_p", ascending=True)

        sig_gap = SIG_LINE_GAP_FRACTION * y_range
        text_gap = SIG_TEXT_GAP_FRACTION * y_range
        line_height = 0.018 * y_range

        start_y = y_data_max + sig_gap

        for i, (_, row) in enumerate(sig_rows.iterrows()):
            x1 = center_map[row["group_1"]]
            x2 = center_map[row["group_2"]]
            y = start_y + i * (sig_gap + text_gap + line_height)

            if P_VALUE_TO_PLOT == "p_adj":
                p_text = format_pvalue(row["p_adj"]).replace("p", "p_adj")
            else:
                p_text = format_pvalue(row["p_raw"])

            add_pairwise_pvalue_bar(
                ax=ax,
                x1=x1,
                x2=x2,
                y=y,
                text=p_text,
                line_height=line_height,
            )

    # Expand y-axis to leave room for p-value labels.
    y_extra = Y_TOP_EXTRA_FRACTION * y_range
    ax.set_ylim(y_data_min - 0.05 * y_range, y_data_max + y_extra)

    legend_handles = [
        Patch(facecolor="gray", edgecolor="black", alpha=REPLICA_ALPHAS[i], label=f"replica {i + 1}")
        for i in range(N_REPLICATES)
    ]

    mean_line_handle = plt.Line2D(
        [0],
        [0],
        color="black",
        linestyle=MEAN_LINE_STYLE,
        linewidth=MEAN_LINE_WIDTH,
        label="mean of 5 trajectory means",
    )

    ax.legend(
        handles=legend_handles + [mean_line_handle],
        frameon=False,
        fontsize=LEGEND_FONTSIZE,
        loc="best",
    )

    for side in ["top", "right"]:
        ax.spines[side].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_file, dpi=300)
    plt.close()


# ============================================================
# 6. Main
# ============================================================

def main():
    all_summary_rows = []
    all_group_mean_rows = []
    all_pvalue_rows = []

    print("\nStatistical method:")
    print("  Unit of analysis: mean value of each independent MD trajectory")
    print(f"  Number of trajectories per system: n = {N_REPLICATES}")
    print("  Pairwise test: two-sided Welch's t-test")
    print("  Multiple-comparison correction: Holm correction within each energy term")
    print("  Figure labels: Holm-adjusted p-values (p_adj)")
    if not SCIPY_AVAILABLE:
        print("  WARNING: scipy is unavailable; using exact permutation p-values as fallback.")

    for _, term_config in ENERGY_TERMS.items():
        print("\n" + "=" * 80)
        print(f"Plotting: {term_config['print_name']}")
        print(f"Section: {term_config['section']} | Column: {term_config['column']}")

        collected = collect_energy_data(term_config)
        pvalue_df = compute_pairwise_pvalues(
            term_config=term_config,
            trajectory_means_by_group=collected["trajectory_means_by_group"],
        )

        plot_boxplot(
            all_data=collected["all_data"],
            positions=collected["positions"],
            box_colors=collected["box_colors"],
            box_alphas=collected["box_alphas"],
            group_centers=collected["group_centers"],
            group_labels=collected["group_labels"],
            trajectory_means_by_group=collected["trajectory_means_by_group"],
            pvalue_df=pvalue_df,
            ylabel=term_config["ylabel"],
            output_file=term_config["output"],
        )

        all_summary_rows.extend(collected["summary_rows"])
        all_group_mean_rows.extend(collected["group_mean_rows"])
        all_pvalue_rows.extend(pvalue_df.to_dict("records"))

        print(f"Saved: {term_config['output']}")

    summary_df = pd.DataFrame(all_summary_rows)
    group_mean_df = pd.DataFrame(all_group_mean_rows)
    pvalue_all_df = pd.DataFrame(all_pvalue_rows)

    summary_df.to_csv("free_energy_boxplot_summary.csv", index=False)
    group_mean_df.to_csv("free_energy_group_mean_lines.csv", index=False)
    pvalue_all_df.to_csv("free_energy_pairwise_pvalues.csv", index=False)

    print("\nSaved: free_energy_boxplot_summary.csv")
    print("Saved: free_energy_group_mean_lines.csv")
    print("Saved: free_energy_pairwise_pvalues.csv")


if __name__ == "__main__":
    main()
