import mdtraj as md
import numpy as np

traj = md.load('center.dcd', top='ref.pdb')
top = traj.topology

def calc_and_save_rmsf(traj, res_start, res_end, filename):

    sel = top.select(f"name CA and resSeq {res_start} to {res_end}")
    if len(sel) == 0:
        raise ValueError(f" no CA in resSeq {res_start} to {res_end}")

    rmsf = md.rmsf(traj, traj[0], atom_indices=sel)

    res_ids = np.array([top.atom(i).residue.resSeq for i in sel], dtype=int)

    data = np.column_stack((res_ids, rmsf))

    np.savetxt(
        filename,
        data,
        delimiter=',',
        header='residue_index,RMSF_nm',
        comments='',          
        fmt=['%d', '%.6f']    
    )

calc_and_save_rmsf(traj, 1,   11,  'rmsf_1_11.csv')
calc_and_save_rmsf(traj, 12,  286, 'rmsf_12_286.csv')
calc_and_save_rmsf(traj, 287, 386, 'rmsf_287_386.csv')

