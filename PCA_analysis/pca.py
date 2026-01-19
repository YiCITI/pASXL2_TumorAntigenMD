import mdtraj as md
import numpy as np


ref1 = md.load_pdb('ref1.pdb')
ref2 = md.load_pdb('ref2.pdb')

traj1 = md.load('traj1.dcd', top=ref1)
traj2 = md.load('traj2.dcd', top=ref2)

ref_frame = traj1[0]


def run_pca_for_region(selection, csv_prefix):

    atom_indices_ref = ref1.topology.select(selection)   
    atom_indices1    = traj1.topology.select(selection)  
    atom_indices2    = traj2.topology.select(selection) 

    if len(atom_indices_ref) == 0:
        raise ValueError(f"No atoms selected：'{selection}'")
    if len(atom_indices1) != len(atom_indices_ref) or len(atom_indices2) != len(atom_indices_ref):
        raise ValueError(
            f"Cannot align '{selection}' , different number of atoms："
            f"ref1={len(atom_indices_ref)}, traj1={len(atom_indices1)}, traj2={len(atom_indices2)}"
        )

    traj1_aligned = traj1.superpose(
        ref_frame,
        atom_indices=atom_indices1,
        ref_atom_indices=atom_indices_ref
    )

    traj2_aligned = traj2.superpose(
        ref_frame,
        atom_indices=atom_indices2,
        ref_atom_indices=atom_indices_ref
    )

    coords1 = traj1_aligned.xyz[:, atom_indices1, :]
    coords2 = traj2_aligned.xyz[:, atom_indices2, :]

    X1 = coords1.reshape(coords1.shape[0], -1)
    X2 = coords2.reshape(coords2.shape[0], -1)

    X = np.vstack([X1, X2])

    X_mean = X.mean(axis=0)
    X_centered = X - X_mean

    U, S, Vt = np.linalg.svd(X_centered, full_matrices=False)
    scores = X_centered @ Vt.T  # (n_total_frames, n_components)

    pc12 = scores[:, :2]

    n1 = X1.shape[0]
    n2 = X2.shape[0]

    pc12_traj1 = pc12[:n1, :]
    pc12_traj2 = pc12[n1:n1 + n2, :]

    np.savetxt(f"{csv_prefix}_traj1.csv", pc12_traj1, delimiter=',')
    np.savetxt(f"{csv_prefix}_traj2.csv", pc12_traj2, delimiter=',')


run_pca_for_region(
    selection="name CA and resSeq >= 1 and resSeq <= 386",
    csv_prefix="pca_all"
)

run_pca_for_region(
    selection="name CA and resSeq >= 1 and resSeq <= 11",
    csv_prefix="pca_1_11"
)
run_pca_for_region(
    selection="name CA and resSeq >= 12 and resSeq <= 286",
    csv_prefix="pca_12_286"
)

run_pca_for_region(
    selection="name CA and resSeq >= 287 and resSeq <= 386",
    csv_prefix="pca_287_386"
)



