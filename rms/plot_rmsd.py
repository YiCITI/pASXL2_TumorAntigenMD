#!/usr/bin/env python3
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
from pathlib import Path


# ============================================================
# Basic settings
# ============================================================

SYSTEMS = ["phos", "prot", "wt"]

TRAJ_FILES = [f"t{i}.dcd" for i in range(1, 6)]
TOP_FILE = "r.pdb"

TOTAL_TRAJ_TIME_NS = 550.0
ANALYSIS_TIME_NS   = 500.0

COLORS = ["red", "orange", "green", "blue", "purple"]

MARKER_SIZE  = 3
MARKER_ALPHA = 0.6
HIST_ALPHA   = 0.35
HIST_BINS    = 30

LABEL_FONTSIZE = 32
TICK_FONTSIZE  = 26

RMSD_YLIM = (0.0, 0.5)
TIME_XLIM = (0.0, 500.0)
HIST_XLIM = (0.0, 1000.0)


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

    time_ns = np.linspace(
        0.0,
        ANALYSIS_TIME_NS,
        traj_prod.n_frames
    )

    return traj_prod, time_ns, start_frame, n_frames_total


def plot_one_system(system):
    system_dir = Path(system)
    top_path = system_dir / TOP_FILE

    if not top_path.exists():
        raise FileNotFoundError(f"Cannot find topology file: {top_path}")

    fig = plt.figure(figsize=(11, 6))
    gs = gridspec.GridSpec(
        1,
        2,
        width_ratios=[4, 1],
        wspace=0.0
    )

    ax_main = plt.subplot(gs[0, 0])
    ax_dist = plt.subplot(gs[0, 1])

    for traj_name, color in zip(TRAJ_FILES, COLORS):
        traj_path = system_dir / traj_name

        if not traj_path.exists():
            raise FileNotFoundError(f"Cannot find trajectory file: {traj_path}")

        traj = md.load(str(traj_path), top=str(top_path))

        atom_indices = traj.topology.select("name CA")
        if len(atom_indices) == 0:
            raise ValueError(f"Cannot find CA atoms in topology: {top_path}")

        traj_prod, time_ns, start_frame, n_frames_total = get_last_500ns(traj)

        rmsd = md.rmsd(
            traj_prod,
            traj_prod,
            frame=0,
            atom_indices=atom_indices
        )

        print(
            f"{system}/{traj_name}: "
            f"total frames = {n_frames_total}, "
            f"start frame = {start_frame + 1}, "
            f"analysis frames = {traj_prod.n_frames}, "
            f"RMSD reference = first frame of last 500 ns"
        )

        ax_main.plot(
            time_ns,
            rmsd,
            marker="o",
            linestyle="None",
            markersize=MARKER_SIZE,
            alpha=MARKER_ALPHA,
            color=color,
            label=traj_name.replace(".dcd", "")
        )

        ax_dist.hist(
            rmsd,
            bins=HIST_BINS,
            orientation="horizontal",
            density=False,
            color=color,
            alpha=HIST_ALPHA,
            edgecolor="black"
        )

    ax_main.set_xlabel("Time (ns)", fontsize=LABEL_FONTSIZE)
    ax_main.set_ylabel("RMSD (nm)", fontsize=LABEL_FONTSIZE)

    ax_main.set_xlim(*TIME_XLIM)
    ax_main.set_ylim(*RMSD_YLIM)

    ax_main.tick_params(
        axis="both",
        labelsize=TICK_FONTSIZE
    )

    ax_main.legend(
        frameon=False,
        fontsize=14,
        loc="upper right"
    )

    ax_dist.set_ylim(ax_main.get_ylim())

    ax_dist.set_xlabel("")
    ax_dist.set_ylabel("")
    ax_dist.set_xticks([])
    ax_dist.set_yticks([])

    ax_dist.tick_params(
        axis="both",
        which="both",
        bottom=False,
        top=False,
        left=False,
        right=False,
        labelbottom=False,
        labelleft=False
    )

    ax_dist.set_xlim(*HIST_XLIM)

    for side in ["top", "bottom", "left", "right"]:
        ax_dist.spines[side].set_visible(False)

    ax_dist.set_frame_on(False)

    ax_main.spines["right"].set_visible(True)

    plt.tight_layout()

    output_name = f"rmsd_{system}.svg"
    plt.savefig(output_name, dpi=300)
    plt.close()

    print(f"Saved: {output_name}")


def main():
    for system in SYSTEMS:
        print("\n" + "=" * 80)
        print(f"Processing system: {system}")
        print("=" * 80)

        plot_one_system(system)


if __name__ == "__main__":
    main()
