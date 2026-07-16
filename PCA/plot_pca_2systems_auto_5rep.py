#!/usr/bin/env python3
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from pathlib import Path


# ============================================================
# 1. Basic settings
# ============================================================

# The script will automatically scan the current working directory.
# A valid system folder must be named phos, prot, or wt and contain:
#   t1.dcd, t2.dcd, t3.dcd, t4.dcd, t5.dcd
#   r1.pdb, r2.pdb, r3.pdb, r4.pdb, r5.pdb
# Exactly two valid folders are expected and will be compared in one PCA space.
SYSTEMS = None

ALLOWED_SYSTEMS = ["phos", "prot", "wt"]

REPLICAS = [1, 2, 3, 4, 5]

TRAJ_TEMPLATE = "t{}.dcd"
TOP_TEMPLATE  = "r{}.pdb"

# Total trajectory = 550 ns, only analyze the last 500 ns
TOTAL_TRAJ_TIME_NS = 550.0
ANALYSIS_TIME_NS   = 500.0

# Default chain order
# chain 0 = peptide
# chain 1 = HLA chain 2
# chain 2 = HLA chain 3
PEPTIDE_CHAIN_INDEX = 0
CHAIN2_INDEX = 1
CHAIN3_INDEX = 2

# Plot style
# Folder names determine colors, consistent with the interaction plotting script:
#   phos = purple, prot = blue, wt = gray.
# Alpha and confidence-ellipse line styles are kept from the original script.
SYSTEM_COLORS = {
    "phos": "#800080",
    "prot": "#0072B2",
    "wt":   "#666666",
}
SYSTEM_ALPHAS = [0.70, 0.35]
SYSTEM_ELLIPSE_LINESTYLES = ["-", "--"]

POINT_SIZE = 10

FONT_SIZE = 20
LABEL_FONTSIZE = 24
TICK_FONTSIZE = 18
TITLE_FONTSIZE = 26

CONF_LEVEL = 0.95
ELLIPSE_LW = 2.2

# Output figures.
# The final output filename will automatically append "{system1}_vs_{system2}".
PCA_TARGETS = [
    ("all",     "All CA atoms",         "PCA_all_CA"),
    ("peptide", "Peptide CA atoms",     "PCA_peptide_1_11_CA"),
    ("chain2",  "HLA chain 2 CA atoms", "PCA_chain2_CA"),
    ("chain3",  "HLA chain 3 CA atoms", "PCA_chain3_CA"),
]


# ============================================================
# 2. Helper functions
# ============================================================

def required_files_for_system(system_dir):
    """
    Return the required trajectory/topology files for one system folder.
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
    return Path(system_dir).is_dir() and all(p.exists() for p in required_files_for_system(system_dir))


def discover_systems(base_dir="."):
    """
    Automatically discover exactly two valid system folders.

    Rules:
    1. Only folders named phos, prot, or wt are considered.
    2. A folder is valid only if it contains all files required by REPLICAS,
       TRAJ_TEMPLATE and TOP_TEMPLATE.
    3. Exactly two valid folders are expected. This avoids silently choosing
       the wrong pair when an old result folder is left in the directory.
    4. The output order follows ALLOWED_SYSTEMS: phos, prot, wt.
    """
    base_dir = Path(base_dir)

    candidate_dirs = [base_dir / name for name in ALLOWED_SYSTEMS if (base_dir / name).is_dir()]
    valid_dirs = [p for p in candidate_dirs if is_valid_system_folder(p)]

    if len(valid_dirs) != 2:
        print("\nDetected allowed folders:")
        for name in ALLOWED_SYSTEMS:
            p = base_dir / name
            if not p.is_dir():
                print(f"  {name}: not found")
                continue

            missing = [
                str(x.relative_to(p))
                for x in required_files_for_system(p)
                if not x.exists()
            ]

            if missing:
                print(f"  {name}: incomplete, missing {missing}")
            else:
                print(f"  {name}: valid")

        raise RuntimeError(
            "\nExactly two valid system folders are required.\n"
            "Folder names must be selected from: phos, prot, wt.\n"
            "Each valid system folder must contain:\n"
            "  t1.dcd, t2.dcd, t3.dcd, t4.dcd, t5.dcd\n"
            "  r1.pdb, r2.pdb, r3.pdb, r4.pdb, r5.pdb\n"
        )

    systems = [p.name for p in valid_dirs]

    print("\nSystems used for PCA:")
    print(f"  system 1: {systems[0]}  color = {SYSTEM_COLORS[systems[0]]}")
    print(f"  system 2: {systems[1]}  color = {SYSTEM_COLORS[systems[1]]}")

    return systems


def safe_name(name):
    """
    Convert a folder name into a safe filename component.
    """
    safe_chars = []
    for ch in str(name):
        if ch.isalnum() or ch in ("-", "_"):
            safe_chars.append(ch)
        else:
            safe_chars.append("_")
    return "".join(safe_chars)


def make_output_svg(output_prefix, systems):
    """
    Build an output filename using the actual system folder names.
    """
    s1 = safe_name(systems[0])
    s2 = safe_name(systems[1])
    return f"{output_prefix}_{s1}_vs_{s2}.svg"


def get_last_500ns(traj):
    """
    Keep only the last 500 ns from a 550 ns trajectory.
    """
    n_frames_total = traj.n_frames

    analysis_frames = int(round(n_frames_total * ANALYSIS_TIME_NS / TOTAL_TRAJ_TIME_NS))
    analysis_frames = min(analysis_frames, n_frames_total)

    start_frame = n_frames_total - analysis_frames

    return traj[start_frame:], start_frame, n_frames_total


def get_ca_indices_by_chain(top, chain_index):
    """
    Select CA atoms by chain index.
    N/C caps do not have CA, so they are automatically excluded.
    """
    chains = list(top.chains)

    if chain_index >= len(chains):
        raise ValueError(
            f"chain_index={chain_index} out of range. "
            f"This topology only has {len(chains)} chains."
        )

    chain = chains[chain_index]

    indices = []
    for res in chain.residues:
        for atom in res.atoms:
            if atom.name == "CA":
                indices.append(atom.index)

    return np.array(indices, dtype=int)


def get_ca_indices(top, region):
    """
    Return CA atom indices for a specified region.
    """
    if region == "all":
        all_indices = []
        for chain in top.chains:
            for res in chain.residues:
                for atom in res.atoms:
                    if atom.name == "CA":
                        all_indices.append(atom.index)
        return np.array(all_indices, dtype=int)

    if region == "peptide":
        return get_ca_indices_by_chain(top, PEPTIDE_CHAIN_INDEX)

    if region == "chain2":
        return get_ca_indices_by_chain(top, CHAIN2_INDEX)

    if region == "chain3":
        return get_ca_indices_by_chain(top, CHAIN3_INDEX)

    raise ValueError(f"Unknown region: {region}")


def _chi2_val_2d(conf_level):
    """
    Chi-square critical value for 2D Gaussian confidence ellipse.
    """
    try:
        from scipy.stats import chi2
        return float(chi2.ppf(conf_level, df=2))
    except Exception:
        fallback = {
            0.68: 2.279,
            0.90: 4.605,
            0.95: 5.991,
            0.99: 9.210,
        }
        return fallback.get(round(conf_level, 2), 5.991)


def add_confidence_ellipse(
    x,
    y,
    ax,
    conf_level=0.95,
    edgecolor="black",
    linestyle="-",
    linewidth=2.0,
    zorder=4,
):
    """
    Add a 2D confidence ellipse.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    if x.size < 3:
        return

    cov = np.cov(x, y)

    if not np.all(np.isfinite(cov)):
        return

    vals, vecs = np.linalg.eigh(cov)

    order = vals.argsort()[::-1]
    vals = vals[order]
    vecs = vecs[:, order]

    angle = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))

    chi2_val = _chi2_val_2d(conf_level)
    scale = np.sqrt(chi2_val)

    width = 2.0 * scale * np.sqrt(max(vals[0], 0.0))
    height = 2.0 * scale * np.sqrt(max(vals[1], 0.0))

    mean_x = np.mean(x)
    mean_y = np.mean(y)

    ell = Ellipse(
        xy=(mean_x, mean_y),
        width=width,
        height=height,
        angle=angle,
        facecolor="none",
        edgecolor=edgecolor,
        linestyle=linestyle,
        linewidth=linewidth,
        zorder=zorder,
    )

    ax.add_patch(ell)


def print_chain_info(top, label):
    """
    Print chain information to verify chain indices.
    """
    print(f"\nTopology chain information for {label}:")
    for chain in top.chains:
        ca_count = 0
        residues = list(chain.residues)

        for res in residues:
            for atom in res.atoms:
                if atom.name == "CA":
                    ca_count += 1

        if residues:
            first_res = residues[0]
            last_res = residues[-1]
            print(
                f"  chain index {chain.index}: "
                f"n_res = {len(residues)}, "
                f"n_CA = {ca_count}, "
                f"first = {first_res.name}{first_res.resSeq}, "
                f"last = {last_res.name}{last_res.resSeq}"
            )
        else:
            print(f"  chain index {chain.index}: empty chain")


# ============================================================
# 3. PCA main function
# ============================================================

def run_pca_for_region(region, title, output_svg, systems):
    """
    Run PCA using all replicas from the two discovered systems combined
    into one PCA space.
    """
    loaded = []

    # Load all trajectories
    for system in systems:
        for rep in REPLICAS:
            traj_file = Path(system) / TRAJ_TEMPLATE.format(rep)
            top_file  = Path(system) / TOP_TEMPLATE.format(rep)

            if not traj_file.exists():
                raise FileNotFoundError(f"Cannot find trajectory file: {traj_file}")

            if not top_file.exists():
                raise FileNotFoundError(f"Cannot find topology file: {top_file}")

            traj = md.load(str(traj_file), top=str(top_file))
            traj, start_frame, n_frames_total = get_last_500ns(traj)

            print(
                f"{system} rep{rep} | {region}: "
                f"total frames = {n_frames_total}, "
                f"start frame = {start_frame + 1}, "
                f"analysis frames = {traj.n_frames}"
            )

            loaded.append({
                "system": system,
                "replica": rep,
                "traj": traj,
                "top": traj.topology,
                "traj_file": traj_file,
                "top_file": top_file,
            })

    # Use the first discovered system, replica 1, first production frame
    # as the common reference.
    ref_item = loaded[0]
    ref_traj = ref_item["traj"]
    ref_frame = ref_traj[0]
    ref_top = ref_item["top"]

    print_chain_info(ref_top, str(ref_item["top_file"]))

    ref_indices = get_ca_indices(ref_top, region)

    if len(ref_indices) == 0:
        raise ValueError(f"No CA atoms selected for region: {region}")

    print(f"\nRegion {region}: reference selected CA atoms = {len(ref_indices)}")

    X_list = []
    labels = []
    replica_labels = []

    # Align each trajectory to the same reference using the same region
    for item in loaded:
        traj = item["traj"]
        top = item["top"]
        system = item["system"]
        rep = item["replica"]

        atom_indices = get_ca_indices(top, region)

        if len(atom_indices) == 0:
            raise ValueError(
                f"No CA atoms selected for region {region} in {item['top_file']}"
            )

        if len(atom_indices) != len(ref_indices):
            raise ValueError(
                f"Different number of selected CA atoms for region {region}.\n"
                f"Reference: {len(ref_indices)} atoms\n"
                f"{item['top_file']}: {len(atom_indices)} atoms\n"
                f"Please check chain order or topology consistency."
            )

        traj_aligned = traj.superpose(
            ref_frame,
            atom_indices=atom_indices,
            ref_atom_indices=ref_indices
        )

        coords = traj_aligned.xyz[:, atom_indices, :]
        X = coords.reshape(coords.shape[0], -1)

        X_list.append(X)
        labels.extend([system] * X.shape[0])
        replica_labels.extend([f"{system}_rep{rep}"] * X.shape[0])

    X_all = np.vstack(X_list)
    labels = np.array(labels)
    replica_labels = np.array(replica_labels)

    # PCA using covariance matrix
    X_mean = X_all.mean(axis=0)
    X_centered = X_all - X_mean

    cov = np.dot(X_centered.T, X_centered) / (X_centered.shape[0] - 1)

    eigvals, eigvecs = np.linalg.eigh(cov)

    order = eigvals.argsort()[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    scores = np.dot(X_centered, eigvecs[:, :2])

    pc1 = scores[:, 0]
    pc2 = scores[:, 1]

    explained = eigvals / eigvals.sum()
    pc1_var = explained[0] * 100.0
    pc2_var = explained[1] * 100.0

    # ========================================================
    # Plot
    # ========================================================
    fig, ax = plt.subplots(figsize=(8, 8))

    # Plot the second system first, then the first system on top.
    # This preserves the original visual logic: comparison/background first,
    # highlighted system second.
    plot_order = [1, 0]

    for idx in plot_order:
        system = systems[idx]
        mask = labels == system

        ax.scatter(
            pc1[mask],
            pc2[mask],
            s=POINT_SIZE,
            alpha=SYSTEM_ALPHAS[idx],
            color=SYSTEM_COLORS[system],
            edgecolors="none",
            label=system,
            zorder=1 + idx
        )

    # Confidence ellipses
    for idx, system in enumerate(systems):
        mask = labels == system

        add_confidence_ellipse(
            pc1[mask],
            pc2[mask],
            ax,
            conf_level=CONF_LEVEL,
            edgecolor="black",
            linestyle=SYSTEM_ELLIPSE_LINESTYLES[idx],
            linewidth=ELLIPSE_LW,
            zorder=4
        )

    # Symmetric axis limits
    max_abs_x = max(abs(pc1.min()), abs(pc1.max()))
    max_abs_y = max(abs(pc2.min()), abs(pc2.max()))
    limit = max(max_abs_x, max_abs_y) * 1.05

    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)

    ax.set_xlabel(f"PC1 ({pc1_var:.1f}%)", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel(f"PC2 ({pc2_var:.1f}%)", fontsize=LABEL_FONTSIZE)
    ax.set_title(title, fontsize=TITLE_FONTSIZE, pad=15)

    ax.tick_params(axis="both", labelsize=TICK_FONTSIZE)

    ax.grid(True, linestyle="--", alpha=0.45)

    ax.legend(
        frameon=False,
        fontsize=FONT_SIZE,
        loc="best"
    )

    plt.tight_layout()
    plt.savefig(output_svg, dpi=300)
    plt.close()

    print(f"Saved: {output_svg}")


# ============================================================
# 4. Main
# ============================================================

def main():
    systems = discover_systems(".")

    for region, title, output_prefix in PCA_TARGETS:
        output_svg = make_output_svg(output_prefix, systems)

        print("\n" + "=" * 80)
        print(f"Running PCA for region: {region}")
        print(f"Output: {output_svg}")
        print("=" * 80)

        run_pca_for_region(region, title, output_svg, systems)


if __name__ == "__main__":
    main()
