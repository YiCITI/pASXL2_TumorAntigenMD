#!/usr/bin/env python3
import re
from pathlib import Path

import mdtraj as md
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


# ============================================================
# 1. Basic settings
# ============================================================

# Exactly two valid folders are expected in the current directory.
# The folder names can only be phos, prot, or wt.
ALLOWED_SYSTEMS = ["phos", "prot", "wt"]

# Five independent MD replicas per system: t1/r1 ... t5/r5
REPLICAS = [1, 2, 3, 4, 5]

# Current trajectories: total 550 ns, only analyze the last 500 ns
TOTAL_TRAJ_TIME_NS = 550.0
ANALYSIS_TIME_NS   = 500.0

TRAJ_TEMPLATE = "t{}.dcd"
TOP_TEMPLATE  = "r{}.pdb"


def safe_filename_text(text):
    """
    Convert a folder/system name into a safe filename component.
    """
    text = str(text).strip()
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text)
    text = text.strip("._-")
    return text or "system"


def required_files_for_system(system_dir):
    """
    Return all required trajectory/topology files for one system folder.
    """
    system_dir = Path(system_dir)

    required = []
    for rep in REPLICAS:
        required.append(system_dir / TRAJ_TEMPLATE.format(rep))
        required.append(system_dir / TOP_TEMPLATE.format(rep))

    return required


def is_valid_system_folder(system_dir):
    """
    A valid system folder contains all required t*.dcd and r*.pdb files.
    """
    system_dir = Path(system_dir)
    return system_dir.is_dir() and all(p.is_file() for p in required_files_for_system(system_dir))


def discover_systems(base_dir="."):
    """
    Automatically detect exactly two valid systems from phos/prot/wt.

    Rules:
    1. Only immediate subfolders named phos, prot, or wt are considered.
    2. A valid system folder must contain t1.dcd-t5.dcd and r1.pdb-r5.pdb.
    3. Exactly two valid folders must be present.
    4. The system order follows ALLOWED_SYSTEMS, so phos/prot/wt ordering is stable.
    """
    base_dir = Path(base_dir)

    valid_systems = []
    print("\nDetected allowed system folders:")

    for system in ALLOWED_SYSTEMS:
        folder = base_dir / system

        if not folder.is_dir():
            print(f"  {system}: not found")
            continue

        missing = [
            str(p.relative_to(folder))
            for p in required_files_for_system(folder)
            if not p.is_file()
        ]

        if missing:
            print(f"  {system}: incomplete, missing {missing}")
        else:
            print(f"  {system}: valid")
            valid_systems.append(system)

    if len(valid_systems) != 2:
        raise RuntimeError(
            "\nExactly two valid system folders are required.\n"
            "Folder names must be selected from: phos, prot, wt.\n"
            "Each valid system folder must contain:\n"
            "  t1.dcd, t2.dcd, t3.dcd, t4.dcd, t5.dcd\n"
            "  r1.pdb, r2.pdb, r3.pdb, r4.pdb, r5.pdb\n"
            f"Detected valid systems: {valid_systems}"
        )

    return valid_systems


SYSTEMS = discover_systems(".")
SYSTEM_1, SYSTEM_2 = SYSTEMS
SYSTEM_TAG = "_vs_".join(safe_filename_text(s) for s in SYSTEMS)

print(f"\nSystems to compare: {SYSTEM_1} vs {SYSTEM_2}")
print(f"Replicas per system: {REPLICAS}")


# ============================================================
# 2. Interaction settings
# ============================================================

HBOND_DIST_CUTOFF = 0.25
HBOND_ANGLE_DEG   = 120.0

SALT_BRIDGE_CUTOFF = 0.40
HYDROPHOBIC_CUTOFF = 0.45

# Residue ranges
# If your real system is 1-386, change 388 -> 386 here.
RES_START, RES_END = 1, 388

# peptide: 1-11
# HLA: 12-388
RANGE_1_START, RANGE_1_END = 1, 11
RANGE_2_START, RANGE_2_END = 12, 388

HYDROPHOBIC_RESNAMES = {
    "ALA", "VAL", "LEU", "ILE", "MET",
    "PHE", "TYR", "TRP", "PRO", "CYS"
}

# Modified phosphorylated residues that mdtraj may not classify as standard protein.
# SEP/TPO/PTR are the commonly used phosphorylated residues.
# S1P/T1P/Y1P are protonated, monoanionic phospho-residue variants in some AMBER libraries.
PHOSPHO_RESNAMES = {"SEP", "TPO", "PTR", "S1P", "T1P", "Y1P"}
TARGET_NON_PROTEIN_RESNAMES = PHOSPHO_RESNAMES

# Phosphate oxygen atom names used by different PDB/AMBER/CHARMM conventions.
# These atoms are treated as possible negatively charged atoms for salt-bridge detection.
PHOSPHO_NEGATIVE_ATOM_NAMES = {
    "O1P", "O2P", "O3P",
    "OP1", "OP2", "OP3",
    "O1", "O2", "O3",
}

# Only keep interactions above this occupancy in each replica.
# Missing interactions are later filled as 0.0 in the summary table.
OCC_THRESHOLD = 0.01

TARGET_TYPES = ["hbond", "salt_bridge", "hydrophobic"]
MAX_INTERACTIONS_PER_TYPE = 5
MAX_PLOTTED_INTERACTIONS = MAX_INTERACTIONS_PER_TYPE * len(TARGET_TYPES)

PAIR_BATCH_SIZE = 20000


# ============================================================
# 3. Plot settings
# ============================================================

# Colors fixed by system name:
# phos = purple, prot = blue, wt = gray.
SYSTEM_COLORS = {
    "phos": "#800080",
    "prot": "#0072B2",
    "wt":   "#666666",
}

AXIS_LABEL_FONTSIZE = 32
TICK_FONTSIZE = 24
TYPE_TITLE_FONTSIZE = 28
PANEL_TITLE_FONTSIZE = 34
LEGEND_FONTSIZE = 22

BAR_ALPHA = 0.85
SHOW_ERRORBARS = True

# Difference heatmap:
# The heatmap is the mean occupancy difference SYSTEM_1 - SYSTEM_2.
# red = positive difference, blue = negative difference, white = 0.
HEATMAP_ABS_LIMIT = 1.0

PLOT_TARGETS = [
    ("peptide_internal", "A. Peptide–Peptide contacts", "A_peptide_peptide_contacts"),
    ("peptide_vs_hla",   "B. Peptide–HLA contacts",    "B_peptide_HLA_contacts"),
    ("hla_internal",     "C. HLA–HLA contacts",        "C_HLA_HLA_contacts"),
]


# ============================================================
# 4. Helper functions
# ============================================================

def get_system_color(system):
    return SYSTEM_COLORS[system]


def is_target_res(res):
    r = res.resSeq

    if not (RES_START <= r <= RES_END):
        return False

    if res.is_protein or res.name in TARGET_NON_PROTEIN_RESNAMES:
        return True

    return False


def classify_pair(res1, res2):
    """
    Classify a residue pair into:
    - peptide_internal
    - peptide_vs_hla
    - hla_internal
    """
    if res1 == res2:
        return None

    if not (is_target_res(res1) and is_target_res(res2)):
        return None

    r1 = res1.resSeq
    r2 = res2.resSeq

    in_pep_1 = RANGE_1_START <= r1 <= RANGE_1_END
    in_pep_2 = RANGE_1_START <= r2 <= RANGE_1_END

    in_hla_1 = RANGE_2_START <= r1 <= RANGE_2_END
    in_hla_2 = RANGE_2_START <= r2 <= RANGE_2_END

    if in_pep_1 and in_pep_2:
        return "peptide_internal"

    if (in_pep_1 and in_hla_2) or (in_pep_2 and in_hla_1):
        return "peptide_vs_hla"

    if in_hla_1 and in_hla_2:
        return "hla_internal"

    return None


def keep_this_pair(res1, res2, target_category):
    return classify_pair(res1, res2) == target_category


def is_sidechain_heavy(atom):
    if atom.element is None:
        return False

    if atom.element.symbol == "H":
        return False

    if atom.name in ("N", "CA", "C", "O", "OXT"):
        return False

    return True


def is_hydrophobic_atom(atom):
    res = atom.residue

    if res.name not in HYDROPHOBIC_RESNAMES:
        return False

    if not is_sidechain_heavy(atom):
        return False

    if atom.element is None or atom.element.symbol != "C":
        return False

    return True


def atom_label(atom):
    res = atom.residue
    return f"{res.name}{res.resSeq} atom {atom.name}"


def get_last_500ns(traj):
    """
    Keep only the last 500 ns from a 550 ns trajectory.
    """
    n_frames_total = traj.n_frames

    analysis_frames = int(round(n_frames_total * ANALYSIS_TIME_NS / TOTAL_TRAJ_TIME_NS))
    analysis_frames = min(analysis_frames, n_frames_total)

    start_frame = n_frames_total - analysis_frames
    traj_prod = traj[start_frame:]

    return traj_prod, start_frame, n_frames_total


def compute_pair_occupancy_in_batches(traj, pairs, cutoff, batch_size=20000):
    pairs = np.asarray(pairs, dtype=int)
    occ_list = []

    for start in range(0, len(pairs), batch_size):
        end = min(start + batch_size, len(pairs))
        sub_pairs = pairs[start:end]

        dist = md.compute_distances(traj, sub_pairs)
        occ = (dist < cutoff).mean(axis=0)

        occ_list.append(occ)

    if not occ_list:
        return np.array([])

    return np.concatenate(occ_list)


def mean_col(system):
    return f"Occ_{system}_mean"


def sd_col(system):
    return f"Occ_{system}_sd"


def rep_col(system, rep):
    return f"{system}_rep{rep}"


# ============================================================
# 5. Analyze one replica for one target category
# ============================================================

def analyze_one_replica(traj_file, top_file, target_category):
    traj = md.load_dcd(str(traj_file), top=str(top_file))
    top = traj.topology

    traj, start_frame, n_frames_total = get_last_500ns(traj)

    print(
        f"{traj_file} | category = {target_category} | "
        f"total frames = {n_frames_total}, "
        f"start frame = {start_frame + 1}, "
        f"analysis frames = {traj.n_frames}"
    )

    rows = []

    def append_row(res1, res2, interaction_str, interaction_type, occ_value):
        if occ_value <= OCC_THRESHOLD:
            return

        cat = classify_pair(res1, res2)

        if cat != target_category:
            return

        rows.append({
            "interaction": interaction_str,
            "interaction_type": interaction_type,
            "category": cat,
            "occupancy": float(occ_value),
        })

    # ------------------------------------------------------------
    # 5.1 Hydrogen bonds
    # ------------------------------------------------------------
    hbonds = md.baker_hubbard(
        traj,
        freq=0.0,
        distance_cutoff=HBOND_DIST_CUTOFF,
        angle_cutoff=HBOND_ANGLE_DEG
    )

    if len(hbonds) > 0:
        hbonds = np.asarray(hbonds, dtype=int)

        # baker_hubbard returns D-H-A
        ha_pairs = hbonds[:, [1, 2]]
        dha_triplets = hbonds

        ha_dist = md.compute_distances(traj, ha_pairs)
        dha_angle = md.compute_angles(traj, dha_triplets)

        angle_cutoff_rad = np.deg2rad(HBOND_ANGLE_DEG)
        present = (ha_dist < HBOND_DIST_CUTOFF) & (dha_angle > angle_cutoff_rad)

        # merge multiple H for the same donor-acceptor pair
        pair_to_indices = {}

        for i, (d_idx, h_idx, a_idx) in enumerate(hbonds):
            key = (int(d_idx), int(a_idx))
            pair_to_indices.setdefault(key, []).append(i)

        for (d_idx, a_idx), idx_list in pair_to_indices.items():
            donor_atom = top.atom(d_idx)
            acceptor_atom = top.atom(a_idx)

            res1 = donor_atom.residue
            res2 = acceptor_atom.residue

            if not keep_this_pair(res1, res2, target_category):
                continue

            sub_present = present[:, idx_list]
            occ_pair = sub_present.any(axis=1).mean()

            interaction_str = f"{atom_label(donor_atom)};{atom_label(acceptor_atom)};"

            append_row(
                res1=res1,
                res2=res2,
                interaction_str=interaction_str,
                interaction_type="hbond",
                occ_value=occ_pair
            )

    # ------------------------------------------------------------
    # 5.2 Salt bridges
    # ------------------------------------------------------------
    nterm_res_indices = set()
    cterm_res_indices = set()

    for chain in top.chains:
        chain_residues = [res for res in chain.residues if is_target_res(res)]

        if not chain_residues:
            continue

        nterm_res_indices.add(chain_residues[0].index)
        cterm_res_indices.add(chain_residues[-1].index)

    neg_atoms = []
    pos_atoms = []

    for res in top.residues:
        if not is_target_res(res):
            continue

        for atom in res.atoms:
            # negative atoms
            if res.name == "ASP" and atom.name in ("OD1", "OD2"):
                neg_atoms.append(atom)

            elif res.name == "GLU" and atom.name in ("OE1", "OE2"):
                neg_atoms.append(atom)

            elif (
                res.name in PHOSPHO_RESNAMES
                and atom.name in PHOSPHO_NEGATIVE_ATOM_NAMES
            ):
                neg_atoms.append(atom)

            if res.index in cterm_res_indices and atom.name in ("O", "OXT"):
                neg_atoms.append(atom)

            # positive atoms
            if res.name == "LYS" and atom.name == "NZ":
                pos_atoms.append(atom)

            elif res.name == "ARG" and atom.name in ("NH1", "NH2", "NE"):
                pos_atoms.append(atom)

            elif res.name in {"HIP", "HID", "HIE", "HIS"} and atom.name in ("ND1", "NE2"):
                pos_atoms.append(atom)

            if res.index in nterm_res_indices and atom.name == "N":
                pos_atoms.append(atom)

    salt_pairs = []
    salt_pair_atoms = []

    for a_neg in neg_atoms:
        for a_pos in pos_atoms:
            res1 = a_neg.residue
            res2 = a_pos.residue

            if not keep_this_pair(res1, res2, target_category):
                continue

            salt_pairs.append([a_neg.index, a_pos.index])
            salt_pair_atoms.append((a_neg, a_pos))

    if len(salt_pairs) > 0:
        occ = compute_pair_occupancy_in_batches(
            traj,
            salt_pairs,
            SALT_BRIDGE_CUTOFF,
            batch_size=PAIR_BATCH_SIZE
        )

        for i, (a1, a2) in enumerate(salt_pair_atoms):
            interaction_str = f"{atom_label(a1)};{atom_label(a2)};"

            append_row(
                res1=a1.residue,
                res2=a2.residue,
                interaction_str=interaction_str,
                interaction_type="salt_bridge",
                occ_value=occ[i]
            )

    # ------------------------------------------------------------
    # 5.3 Hydrophobic interactions
    # ------------------------------------------------------------
    hydrophobic_atoms = []

    for res in top.residues:
        if not is_target_res(res):
            continue

        for atom in res.atoms:
            if is_hydrophobic_atom(atom):
                hydrophobic_atoms.append(atom)

    hydrophobic_pairs = []
    hydrophobic_pair_atoms = []

    for i in range(len(hydrophobic_atoms)):
        a1 = hydrophobic_atoms[i]

        for j in range(i + 1, len(hydrophobic_atoms)):
            a2 = hydrophobic_atoms[j]

            if not keep_this_pair(a1.residue, a2.residue, target_category):
                continue

            hydrophobic_pairs.append([a1.index, a2.index])
            hydrophobic_pair_atoms.append((a1, a2))

    if len(hydrophobic_pairs) > 0:
        occ = compute_pair_occupancy_in_batches(
            traj,
            hydrophobic_pairs,
            HYDROPHOBIC_CUTOFF,
            batch_size=PAIR_BATCH_SIZE
        )

        for i, (a1, a2) in enumerate(hydrophobic_pair_atoms):
            interaction_str = f"{atom_label(a1)};{atom_label(a2)};"

            append_row(
                res1=a1.residue,
                res2=a2.residue,
                interaction_str=interaction_str,
                interaction_type="hydrophobic",
                occ_value=occ[i]
            )

    df = pd.DataFrame(
        rows,
        columns=["interaction", "interaction_type", "category", "occupancy"]
    )

    if not df.empty:
        df = (
            df.groupby(["interaction", "interaction_type", "category"], as_index=False)
              ["occupancy"]
              .max()
        )

        df = df.sort_values(
            ["interaction_type", "occupancy"],
            ascending=[True, False]
        )

    return df


# ============================================================
# 6. Collect replicas and summarize
# ============================================================

def collect_all_replicas(target_category):
    all_dfs = []

    for system in SYSTEMS:
        for rep in REPLICAS:
            traj_file = Path(system) / TRAJ_TEMPLATE.format(rep)
            top_file  = Path(system) / TOP_TEMPLATE.format(rep)

            if not traj_file.exists():
                raise FileNotFoundError(f"Cannot find trajectory: {traj_file}")

            if not top_file.exists():
                raise FileNotFoundError(f"Cannot find topology: {top_file}")

            df = analyze_one_replica(traj_file, top_file, target_category)

            df["system"] = system
            df["replica"] = rep

            all_dfs.append(df)

    all_df = pd.concat(all_dfs, ignore_index=True)
    return all_df


def make_summary_table(all_df):
    key_cols = ["interaction", "interaction_type", "category"]
    keys = all_df[key_cols].drop_duplicates().reset_index(drop=True)
    summary = keys.copy()

    for system in SYSTEMS:
        for rep in REPLICAS:
            sub = all_df[
                (all_df["system"] == system) &
                (all_df["replica"] == rep)
            ][key_cols + ["occupancy"]].copy()

            col_name = rep_col(system, rep)
            sub = sub.rename(columns={"occupancy": col_name})

            summary = summary.merge(
                sub,
                on=key_cols,
                how="left"
            )

            summary[col_name] = summary[col_name].fillna(0.0)

    for system in SYSTEMS:
        cols = [rep_col(system, rep) for rep in REPLICAS]

        # Bar height: mean occupancy across five replicas.
        summary[mean_col(system)] = summary[cols].mean(axis=1)

        # Error bar: standard deviation across five replicas.
        summary[sd_col(system)] = summary[cols].std(axis=1, ddof=1)

    summary["Diff_mean"] = summary[mean_col(SYSTEM_1)] - summary[mean_col(SYSTEM_2)]
    summary["Abs_Diff_mean"] = summary["Diff_mean"].abs()

    return summary


def select_top_interactions(summary, target_category):
    selected_list = []

    for t in TARGET_TYPES:
        df_type = summary[summary["interaction_type"] == t].copy()
        df_type = df_type.sort_values("Abs_Diff_mean", ascending=False)
        df_top = df_type.head(MAX_INTERACTIONS_PER_TYPE).copy()

        count = len(df_top)

        if count < MAX_INTERACTIONS_PER_TYPE:
            n_pad = MAX_INTERACTIONS_PER_TYPE - count

            pad_data = {
                "interaction": [np.nan] * n_pad,
                "interaction_type": [t] * n_pad,
                "category": [target_category] * n_pad,
                "Diff_mean": [0.0] * n_pad,
                "Abs_Diff_mean": [0.0] * n_pad,
            }

            for system in SYSTEMS:
                pad_data[mean_col(system)] = [0.0] * n_pad
                pad_data[sd_col(system)] = [0.0] * n_pad
                for rep in REPLICAS:
                    pad_data[rep_col(system, rep)] = [0.0] * n_pad

            pad = pd.DataFrame(pad_data)
            df_top = pd.concat([df_top, pad], ignore_index=True)

        selected_list.append(df_top)

    selected = pd.concat(selected_list, ignore_index=True)
    return selected


# ============================================================
# 7. Plot
# ============================================================

def plot_selected_interactions(df_top_final, panel_title, output_pdf):
    n_systems = len(SYSTEMS)
    bar_group_width = 0.80
    bar_width = bar_group_width / n_systems

    ind = np.arange(MAX_PLOTTED_INTERACTIONS)
    offsets = (np.arange(n_systems) - (n_systems - 1) / 2) * bar_width

    fig = plt.figure(figsize=(24, 8))

    gs = fig.add_gridspec(
        2, 3,
        height_ratios=[8, 0.5],
        width_ratios=[15, 0.5, 1],
        hspace=0,
        wspace=0.1
    )

    ax_bar = fig.add_subplot(gs[0, 0])
    ax_heatmap = fig.add_subplot(gs[1, 0], sharex=ax_bar)
    ax_cbar = fig.add_subplot(gs[:, 2])

    # ------------------------------------------------------------
    # Bar plot: height = mean occupancy; error bar = SD among replicas.
    # ------------------------------------------------------------
    for idx, system in enumerate(SYSTEMS):
        x_pos = ind + offsets[idx]
        means = df_top_final[mean_col(system)].values

        if SHOW_ERRORBARS:
            yerr = df_top_final[sd_col(system)].values
        else:
            yerr = None

        ax_bar.bar(
            x_pos,
            means,
            bar_width,
            yerr=yerr,
            capsize=4 if SHOW_ERRORBARS else 0,
            color=get_system_color(system),
            edgecolor="black",
            linewidth=0.9,
            alpha=BAR_ALPHA,
            label=system
        )

    ax_bar.set_ylabel("Occupancy", fontsize=AXIS_LABEL_FONTSIZE)
    ax_bar.set_ylim(0, 1.2)
    ax_bar.tick_params(axis="y", labelsize=TICK_FONTSIZE)
    ax_bar.tick_params(axis="x", length=0)
    plt.setp(ax_bar.get_xticklabels(), visible=False)

    ax_bar.set_title(panel_title, fontsize=PANEL_TITLE_FONTSIZE, pad=20)

    # Type labels and vertical dividers.
    current_x = 0
    for t in TARGET_TYPES:
        center_x = current_x + (MAX_INTERACTIONS_PER_TYPE - 1) / 2

        if current_x > 0:
            ax_bar.axvline(
                x=current_x - bar_group_width / 2,
                color="k",
                linestyle="--",
                linewidth=1
            )

        ax_bar.text(
            center_x,
            ax_bar.get_ylim()[1] * 0.95,
            t.replace("_", " ").title(),
            ha="center",
            va="top",
            fontsize=TYPE_TITLE_FONTSIZE,
            bbox=dict(facecolor="white", alpha=0.75, edgecolor="none")
        )

        current_x += MAX_INTERACTIONS_PER_TYPE

    ax_bar.legend(
        frameon=False,
        fontsize=LEGEND_FONTSIZE,
        loc="upper right",
        ncol=1,
    )

    # ------------------------------------------------------------
    # Difference heatmap: SYSTEM_1 - SYSTEM_2.
    # ------------------------------------------------------------
    heatmap_data = df_top_final["Diff_mean"].values.reshape(1, -1)

    cmap_heatmap = LinearSegmentedColormap.from_list(
        "BlueWhiteRed_Diff",
        [
            (0.0, "blue"),
            (0.5, "white"),
            (1.0, "red")
        ]
    )

    abs_limit = max(
        HEATMAP_ABS_LIMIT,
        float(np.nanmax(np.abs(heatmap_data))) if heatmap_data.size else HEATMAP_ABS_LIMIT
    )

    im = ax_heatmap.imshow(
        heatmap_data,
        cmap=cmap_heatmap,
        aspect="auto",
        vmin=-abs_limit,
        vmax=abs_limit,
        extent=[
            -bar_group_width / 2,
            MAX_PLOTTED_INTERACTIONS - bar_group_width / 2,
            0,
            1
        ]
    )

    x_labels = (
        df_top_final["interaction"]
        .fillna(" ")
        .str.replace(";", " - ", regex=False)
    )

    ax_heatmap.set_xticks(ind)
    ax_heatmap.set_xticklabels(
        x_labels,
        rotation=60,
        ha="right",
        fontsize=TICK_FONTSIZE
    )

    ax_heatmap.set_xlim(
        -bar_group_width / 2,
        MAX_PLOTTED_INTERACTIONS - bar_group_width / 2
    )

    ax_heatmap.set_yticks([])
    ax_heatmap.set_ylabel("")
    ax_heatmap.tick_params(axis="x", length=0, labelsize=TICK_FONTSIZE)

    cbar = fig.colorbar(im, cax=ax_cbar, orientation="vertical")
    cbar.ax.tick_params(labelsize=TICK_FONTSIZE)
    cbar.set_label(
        f"Mean difference\n({SYSTEM_1} - {SYSTEM_2})",
        fontsize=AXIS_LABEL_FONTSIZE
    )

    for side in ["top", "right"]:
        ax_bar.spines[side].set_visible(False)

    plt.tight_layout()
    plt.savefig(output_pdf, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Saved figure: {output_pdf}")


# ============================================================
# 8. Main
# ============================================================

def main():
    for target_category, panel_title, output_prefix in PLOT_TARGETS:
        print("\n" + "=" * 80)
        print(f"Processing category: {target_category}")
        print("=" * 80)

        output_pdf = f"{output_prefix}_{SYSTEM_TAG}.pdf"
        summary_csv = f"{target_category}_interaction_occupancy_summary_{SYSTEM_TAG}.csv"
        selected_csv = f"{target_category}_top_interactions_{SYSTEM_TAG}.csv"

        all_df = collect_all_replicas(target_category)
        summary = make_summary_table(all_df)
        selected = select_top_interactions(summary, target_category)

        summary.to_csv(summary_csv, index=False)
        selected.to_csv(selected_csv, index=False)

        print(f"Saved summary: {summary_csv}")
        print(f"Saved selected interactions: {selected_csv}")

        plot_selected_interactions(
            df_top_final=selected,
            panel_title=panel_title,
            output_pdf=output_pdf
        )


if __name__ == "__main__":
    main()
