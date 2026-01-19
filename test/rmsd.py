import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

TRAJ_FILE = 'center.dcd'
TOP_FILE  = 'ref.pdb'

TOTAL_TIME_NS = 500.0  

COLOR_PURPLE = '#808080'
MARKER_SIZE  = 3
MARKER_ALPHA = 0.6
HIST_BINS    = 30


LABEL_FONTSIZE = 32   
TICK_FONTSIZE  = 26   
# TITLE_FONTSIZE = 20  


traj = md.load(TRAJ_FILE, top=TOP_FILE)
atom_indices = traj.topology.select('name CA')
if len(atom_indices) == 0:
    raise ValueError("canot find CA in topology")

rmsd = md.rmsd(traj, traj, frame=0, atom_indices=atom_indices)  
n_frames = traj.n_frames

time_ns = np.linspace(0.0, TOTAL_TIME_NS, n_frames)
print(f"总帧数: {n_frames}, 时间轴范围: 0 ~ {TOTAL_TIME_NS:.1f} ns")

fig = plt.figure(figsize=(11, 6))
gs = gridspec.GridSpec(1, 2, width_ratios=[4, 1], wspace=0.0)

ax_main = plt.subplot(gs[0, 0])
ax_dist = plt.subplot(gs[0, 1])

ax_main.plot(
    time_ns, rmsd,
    marker='o',
    linestyle='None',
    markersize=MARKER_SIZE,
    alpha=MARKER_ALPHA,
    color=COLOR_PURPLE
)

ax_main.set_xlabel('Time (ns)', fontsize=LABEL_FONTSIZE)
ax_main.set_ylabel('RMSD (nm)', fontsize=LABEL_FONTSIZE)

ax_main.set_xlim(0.0, TOTAL_TIME_NS)

ax_main.set_ylim(0.0, 0.5)

ax_main.tick_params(axis='both', labelsize=TICK_FONTSIZE)

ax_dist.hist(
    rmsd,
    bins=HIST_BINS,
    orientation='horizontal',
    density=False,            
    color=COLOR_PURPLE,
    alpha=0.5,
    edgecolor='black'
)

ax_dist.set_ylim(ax_main.get_ylim())

ax_dist.set_yticks([])

ax_dist.set_xlabel('')
ax_dist.set_ylabel('')
ax_dist.set_xticks([])
ax_dist.set_yticks([])
ax_dist.tick_params(
    axis='both', which='both',
    bottom=False, top=False, left=False, right=False,
    labelbottom=False, labelleft=False
)

ax_dist.set_xlim(0.0, 1000.0)

for side in ['top', 'bottom', 'left', 'right']:
    ax_dist.spines[side].set_visible(False)
ax_dist.set_frame_on(False)

ax_main.spines['right'].set_visible(True)

plt.tight_layout()
plt.savefig('new_rmsd.svg', dpi=300)
