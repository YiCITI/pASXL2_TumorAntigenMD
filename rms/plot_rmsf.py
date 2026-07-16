#!/usr/bin/env python3
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


# ============================================================
# Basic settings
# ============================================================

SYSTEMS = ["phos", "prot", "wt"]

SYSTEM_COLORS = {
    "phos": "purple",
    "prot": "blue",
    "wt": "grey",
}

SYSTEM_LABELS = {
    "phos": "Phosphorylated",
    "prot": "Protonated",
    "wt": "WT",
}

REPLICAS = [1, 2, 3, 4, 5]

TOP_FILE = "r.pdb"
TRAJ_TEMPLATE = "t{}.dcd"

TOTAL_TRAJ_TIME_NS = 550.0
ANALYSIS_TIME_NS   = 500.0
RELAX_TIME_NS      = TOTAL_TRAJ_TIME_NS - ANALYSIS_TIME_NS

# 三个片段
# 第二条链去掉 N/C cap，所以用 13-287
CHAIN_RANGES = [
    (1,   11,  "Peptide"),
    (13,  287, "HLA chain 2"),
    (289, 388, "HLA chain 3"),
]

# 用于整体叠合的 CA 范围，同样排除 cap
ALIGN_RANGES = [
    (1,   11),
    (13,  287),
    (289, 388),
]

LINE_ALPHA = 0.85
LINE_WIDTH = 2.0
MARKER_SIZE = 3.5
ERROR_ALPHA = 0.45
ERROR_LINEWIDTH = 0.8
CAPSIZE = 1.5

LABEL_FONTSIZE = 32
TICK_FONTSIZE  = 22
LEGEND_FONTSIZE = 20

Y_LIM = (0.0, 1.0)

OUTPUT_FIG = "rmsf_phos_prot_wt_mean_sd.svg"


# ============================================================
# Helper functions
# ============================================================

def select_ca_by_resseq(top, res_start, res_end):
    sel = top.select(f"name CA and resSeq {res_start} to {res_end}")

    if len(sel) == 0:
        raise ValueError(f"No CA atoms found in resSeq {res_start} to {res_end}")

    res_ids = np.array(
        [top.atom(i).residue.resSeq for i in sel],
        dtype=int
    )

    return sel, res_ids


def select_ca_by_ranges(top, ranges):
    all_indices = []

    for res_start, res_end in ranges:
        sel, _ = select_ca_by_resseq(top, res_start, res_end)
        all_indices.extend(sel)

    return np.array(all_indices, dtype=int)


def make_ticks(res_start, res_end, tick_step):
    ticks = list(np.arange(res_start, res_end + 1, tick_step))

    if res_start not in ticks:
        ticks.insert(0, res_start)

    if res_end not in ticks:
        ticks.append(res_end)

    return ticks


def get_last_500ns(traj):
    """
    对总长 550 ns 的轨迹，只保留最后 500 ns。
    如果轨迹是 5500 帧，则保留最后 5000 帧。
    """
    n_frames_total = traj.n_frames

    analysis_frames = int(round(n_frames_total * ANALYSIS_TIME_NS / TOTAL_TRAJ_TIME_NS))
    analysis_frames = min(analysis_frames, n_frames_total)

    start_frame = n_frames_total - analysis_frames

    traj_prod = traj[start_frame:]

    return traj_prod, start_frame, n_frames_total


def compute_rmsf_one_traj(traj_path, top_path):
    """
    对一条 trajectory 计算三个片段的 RMSF。
    返回：
    {
        (res_start, res_end): (res_ids, rmsf)
    }
    """
    traj = md.load(str(traj_path), top=str(top_path))
    top = traj.topology

    traj_prod, start_frame, n_frames_total = get_last_500ns(traj)

    align_indices = select_ca_by_ranges(top, ALIGN_RANGES)

    traj_prod.superpose(
        traj_prod,
        frame=0,
        atom_indices=align_indices
    )

    print(
        f"{traj_path}: total frames = {n_frames_total}, "
        f"start frame = {start_frame + 1}, "
        f"analysis frames = {traj_prod.n_frames}"
    )

    result = {}

    for res_start, res_end, chain_name in CHAIN_RANGES:
        rmsf_indices, res_ids = select_ca_by_resseq(
            top,
            res_start,
            res_end
        )

        rmsf = md.rmsf(
            traj_prod,
            traj_prod,
            frame=0,
            atom_indices=rmsf_indices
        )

        result[(res_start, res_end)] = (res_ids, rmsf)

    return result


def compute_system_rmsf(system):
    """
    对某个体系的 5 条 replica 计算 RMSF，并返回 mean ± SD。
    """
    system_dir = Path(system)
    top_path = system_dir / TOP_FILE

    if not top_path.exists():
        raise FileNotFoundError(f"Cannot find topology file: {top_path}")

    replica_results = []

    for rep in REPLICAS:
        traj_path = system_dir / TRAJ_TEMPLATE.format(rep)

        if not traj_path.exists():
            raise FileNotFoundError(f"Cannot find trajectory file: {traj_path}")

        one_result = compute_rmsf_one_traj(traj_path, top_path)
        replica_results.append(one_result)

    summary = {}

    for res_start, res_end, chain_name in CHAIN_RANGES:
        key = (res_start, res_end)

        ref_res_ids = replica_results[0][key][0]

        rmsf_list = []

        for one_result in replica_results:
            res_ids, rmsf = one_result[key]

            if len(res_ids) != len(ref_res_ids) or not np.all(res_ids == ref_res_ids):
                raise ValueError(
                    f"Residue ID mismatch in {system}, range {res_start}-{res_end}"
                )

            rmsf_list.append(rmsf)

        rmsf_array = np.vstack(rmsf_list)

        rmsf_mean = rmsf_array.mean(axis=0)
        rmsf_sd   = rmsf_array.std(axis=0, ddof=1)

        summary[key] = {
            "res_ids": ref_res_ids,
            "mean": rmsf_mean,
            "sd": rmsf_sd,
        }

    return summary


# ============================================================
# Main plotting
# ============================================================

def main():
    all_summary = {}

    for system in SYSTEMS:
        print("\n" + "=" * 80)
        print(f"Processing system: {system}")
        print("=" * 80)

        all_summary[system] = compute_system_rmsf(system)

    fig, axes = plt.subplots(
        nrows=3,
        ncols=1,
        figsize=(18, 8),
        sharey=True
    )

    for ax, (res_start, res_end, chain_name) in zip(axes, CHAIN_RANGES):
        key = (res_start, res_end)

        for system in SYSTEMS:
            res_ids = all_summary[system][key]["res_ids"]
            rmsf_mean = all_summary[system][key]["mean"]
            rmsf_sd = all_summary[system][key]["sd"]

            color = SYSTEM_COLORS[system]
            label = SYSTEM_LABELS[system]

            ax.errorbar(
                res_ids,
                rmsf_mean,
                yerr=rmsf_sd,
                fmt="o-",
                color=color,
                ecolor=color,
                elinewidth=ERROR_LINEWIDTH,
                capsize=CAPSIZE,
                markersize=MARKER_SIZE,
                linewidth=LINE_WIDTH,
                alpha=LINE_ALPHA,
                label=label
            )

        ax.set_xlim(res_start, res_end)
        ax.set_ylim(*Y_LIM)

        ax.set_ylabel("RMSF (nm)", fontsize=LABEL_FONTSIZE)

        ax.tick_params(
            axis="both",
            labelsize=TICK_FONTSIZE,
            width=1.5,
            length=6
        )

        ax.text(
            0.02,
            0.82,
            chain_name,
            transform=ax.transAxes,
            fontsize=24,
            verticalalignment="top"
        )

        if res_start == 1 and res_end == 11:
            xticks = np.arange(res_start, res_end + 1, 1)
        elif res_start == 13 and res_end == 287:
            xticks = make_ticks(res_start, res_end, tick_step=50)
        else:
            xticks = make_ticks(res_start, res_end, tick_step=20)

        ax.set_xticks(xticks)

        for side in ["top", "right"]:
            ax.spines[side].set_visible(False)

    axes[-1].set_xlabel("Residue index", fontsize=LABEL_FONTSIZE)

    axes[0].legend(
        frameon=False,
        fontsize=LEGEND_FONTSIZE,
        loc="upper right",
        ncol=3
    )

    plt.tight_layout(h_pad=0.8)
    plt.savefig(OUTPUT_FIG, dpi=300)
    plt.close()

    print(f"Saved: {OUTPUT_FIG}")


if __name__ == "__main__":
    main()
