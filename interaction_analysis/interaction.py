#!/usr/bin/env python3
import mdtraj as md
import numpy as np
import pandas as pd


HBOND_DIST_CUTOFF = 0.25    
HBOND_ANGLE_DEG   = 120.0   

SALT_BRIDGE_CUTOFF = 0.40   
HYDROPHOBIC_CUTOFF = 0.45  

RES_START, RES_END = 1, 386

RANGE_1_START, RANGE_1_END = 1, 11
RANGE_2_START, RANGE_2_END = 12, 386

HYDROPHOBIC_RESNAMES = {
    'ALA', 'VAL', 'LEU', 'ILE', 'MET',
    'PHE', 'TYR', 'TRP', 'PRO', 'CYS'
}

OCC_THRESHOLD = 0.01

TARGET_NON_PROTEIN_RESNAMES = {'SEP'}


def is_target_res(res):
    """判断一个残基是否是目标残基（在 resSeq 范围内，且是蛋白质或目标非标准残基）"""
    r = res.resSeq
    if not (RES_START <= r <= RES_END):
        return False
    
    if res.is_protein or res.name in TARGET_NON_PROTEIN_RESNAMES:
        return True
    
    return False


def in_target_ranges(res1, res2):
    """
    判断两个残基是否在统一范围内：
    - 两者都不是同一个残基
    - 两者都是目标残基（使用 is_target_res 检查）
    """
    if res1 == res2:
        return False

    if not (is_target_res(res1) and is_target_res(res2)):
        return False

    return True


def classify_pair(res1, res2):

    r1 = res1.resSeq
    r2 = res2.resSeq

    in1_1 = RANGE_1_START <= r1 <= RANGE_1_END
    in1_2 = RANGE_1_START <= r2 <= RANGE_1_END
    in2_1 = RANGE_2_START <= r1 <= RANGE_2_END
    in2_2 = RANGE_2_START <= r2 <= RANGE_2_END

    if in1_1 and in1_2:
        return "1-11_internal"
    if (in1_1 and in2_2) or (in1_2 and in2_1):
        return "1-11_vs_12-386"
    if in2_1 and in2_2:
        return "12-386_internal"
    return None


def is_sidechain_heavy(atom):

    if atom.element is None:
        return False
    if atom.element.symbol == 'H':
        return False
    if atom.name in ('N', 'CA', 'C', 'O', 'OXT'):
        return False
    return True


def is_hydrophobic_atom(atom):

    res = atom.residue
    if res.name not in HYDROPHOBIC_RESNAMES:
        return False
    if not is_sidechain_heavy(atom):
        return False
    if atom.element is None or atom.element.symbol != 'C':
        return False
    return True


def atom_label(atom):

    res = atom.residue
    return f"{res.name}{res.resSeq} atom {atom.name}"


# ------------------ 主分析函数 ------------------ #
def main():
    # 1. 读轨迹
    traj = md.load_dcd('traj.dcd', top='ref.pdb')
    top = traj.topology
    n_frames = traj.n_frames

    # 3 个文件分别对应的记录列表
    rows_1_11_internal = []
    rows_1_11_vs_12_386 = []
    rows_12_386_internal = []

    def append_row_by_category(res1, res2, interaction_str, interaction_type, occ_value):
       
        cat = classify_pair(res1, res2)
        if cat is None:
            return
        if occ_value <= OCC_THRESHOLD:
            return

        row = {
            "interaction": interaction_str,
            "interaction_type": interaction_type,
            "occupancy": float(occ_value)
        }

        if cat == "1-11_internal":
            rows_1_11_internal.append(row)
        elif cat == "1-11_vs_12-386":
            rows_1_11_vs_12_386.append(row)
        elif cat == "12-386_internal":
            rows_12_386_internal.append(row)

    # ------------------ 2. 氢键 (hbond) ------------------ #
    hbonds = md.baker_hubbard(
        traj,
        freq=0.0,
        distance_cutoff=HBOND_DIST_CUTOFF,
        angle_cutoff=HBOND_ANGLE_DEG
    )

    if len(hbonds) > 0:
        hbonds = np.asarray(hbonds, dtype=int)      # shape: (n_hbonds, 3)
        ha_pairs = hbonds[:, [1, 2]]                # (H, A)
        dha_triplets = hbonds                       # (D, H, A)

        ha_dist = md.compute_distances(traj, ha_pairs)
        dha_angle = md.compute_angles(traj, dha_triplets)
        angle_cutoff_rad = np.deg2rad(HBOND_ANGLE_DEG)

        present = (ha_dist < HBOND_DIST_CUTOFF) & (dha_angle > angle_cutoff_rad)
       

       
        pair_to_indices = {}
        for i, (d_idx, h_idx, a_idx) in enumerate(hbonds):
            key = (int(d_idx), int(a_idx))  
            pair_to_indices.setdefault(key, []).append(i)

        for (d_idx, a_idx), idx_list in pair_to_indices.items():
            donor_atom = top.atom(d_idx)
            acceptor_atom = top.atom(a_idx)

            res1 = donor_atom.residue
            res2 = acceptor_atom.residue

            
            if not in_target_ranges(res1, res2):
                continue

           
            sub_present = present[:, idx_list]       
            occ_pair = sub_present.any(axis=1).mean()  

            interaction_str = f"{atom_label(donor_atom)};{atom_label(acceptor_atom)};"

            append_row_by_category(
                res1, res2,
                interaction_str=interaction_str,
                interaction_type="hbond",
                occ_value=occ_pair
            )


    nterm_res_indices = set()
    cterm_res_indices = set()
    for chain in top.chains:
       
        chain_residues = [res for res in chain.residues
                              if is_target_res(res)] 
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
            
            if res.name == 'ASP' and atom.name in ('OD1', 'OD2'):
                neg_atoms.append(atom)
            elif res.name == 'GLU' and atom.name in ('OE1', 'OE2'):
                neg_atoms.append(atom)
            elif res.name == 'SEP' and atom.name in ('O1P', 'O2P', 'O3P'):
                neg_atoms.append(atom)

            if res.index in cterm_res_indices and atom.name in ('O', 'OXT'):
                neg_atoms.append(atom)

            
            if res.name == 'LYS' and atom.name == 'NZ':
                pos_atoms.append(atom)
            elif res.name == 'ARG' and atom.name in ('NH1', 'NH2', 'NE'):
                pos_atoms.append(atom)
            elif res.name in {'HIP', 'HID', 'HIE', 'HIS'} and atom.name in ('ND1', 'NE2'):
                pos_atoms.append(atom)

            if res.index in nterm_res_indices and atom.name == 'N':
                pos_atoms.append(atom)

    salt_pairs = []
    salt_pair_atoms = []

    for a_neg in neg_atoms:
        for a_pos in pos_atoms:
            res1 = a_neg.residue
            res2 = a_pos.residue
            
            if not in_target_ranges(res1, res2):
                continue
            salt_pairs.append([a_neg.index, a_pos.index])
            salt_pair_atoms.append((a_neg, a_pos))

    if len(salt_pairs) > 0:
        salt_pairs = np.asarray(salt_pairs, dtype=int)
        dist = md.compute_distances(traj, salt_pairs)
        occ = (dist < SALT_BRIDGE_CUTOFF).mean(axis=0)

        for i, (a1, a2) in enumerate(salt_pair_atoms):
            res1 = a1.residue
            res2 = a2.residue
            interaction_str = f"{atom_label(a1)};{atom_label(a2)};"

            append_row_by_category(
                res1, res2,
                interaction_str=interaction_str,
                interaction_type="salt_bridge",
                occ_value=occ[i]
            )

    
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
            res1 = a1.residue
            res2 = a2.residue
            
            if not in_target_ranges(res1, res2):
                continue
            hydrophobic_pairs.append([a1.index, a2.index])
            hydrophobic_pair_atoms.append((a1, a2))

    if len(hydrophobic_pairs) > 0:
        hydrophobic_pairs = np.asarray(hydrophobic_pairs, dtype=int)
        dist = md.compute_distances(traj, hydrophobic_pairs)
        occ = (dist < HYDROPHOBIC_CUTOFF).mean(axis=0)

        for i, (a1, a2) in enumerate(hydrophobic_pair_atoms):
            res1 = a1.residue
            res2 = a2.residue
            interaction_str = f"{atom_label(a1)};{atom_label(a2)};"

            append_row_by_category(
                res1, res2,
                interaction_str=interaction_str,
                interaction_type="hydrophobic",
                occ_value=occ[i]
            )

   
    df1 = pd.DataFrame(rows_1_11_internal,
                       columns=["interaction", "interaction_type", "occupancy"])
    df2 = pd.DataFrame(rows_1_11_vs_12_386,
                       columns=["interaction", "interaction_type", "occupancy"])
    df3 = pd.DataFrame(rows_12_386_internal,
                       columns=["interaction", "interaction_type", "occupancy"])

    
    if not df1.empty:
        df1 = df1.sort_values("occupancy", ascending=False)
    if not df2.empty:
        df2 = df2.sort_values("occupancy", ascending=False)
    if not df3.empty:
        df3 = df3.sort_values("occupancy", ascending=False)

    df1.to_csv("interactions_1_11_internal.csv",
               index=False,
               float_format="%.6f")
    df2.to_csv("interactions_1_11_vs_12_386.csv",
               index=False,
               float_format="%.6f")
    df3.to_csv("interactions_12_386_internal.csv",
               index=False,
               float_format="%.6f")

if __name__ == "__main__":
    main()