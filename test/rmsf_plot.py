import numpy as np
import matplotlib.pyplot as plt


COLOR_SITE_4 = '#800080'    
COLOR_NO_PHOS = '#808080' 

TICK_FONTSIZE = 18     
LABEL_FONTSIZE = 22    
LEGEND_FONTSIZE = 18   

def load_rmsf_csv(filename):
   
    data = np.loadtxt(filename, delimiter=',', skiprows=1)
    residues = data[:, 0]
    rmsf = data[:, 1]
    return residues, rmsf

try:
    
    res1, rmsf1 = load_rmsf_csv('1_1.csv')  
    res2, rmsf2 = load_rmsf_csv('2_1.csv')  
    res3, rmsf3 = load_rmsf_csv('3_1.csv')  

    plt.figure(figsize=(8, 5))

   
    plt.plot(res1, rmsf1, marker='o', linestyle='-', color=COLOR_SITE_4,
             label='Phosphorylated', zorder=3)

   
    plt.plot(res2, rmsf2, marker='o', linestyle='-', color=COLOR_NO_PHOS,
             label='No phosphorylated', zorder=2)

    
    plt.plot(res3, rmsf3, marker='s', linestyle='--', color=COLOR_SITE_4,
             label='Protonated-phosphorylated', zorder=3)

    
    plt.xlabel('Residue index', fontsize=LABEL_FONTSIZE)
    plt.ylabel('RMSF (nm)', fontsize=LABEL_FONTSIZE)

    plt.ylim(0, 0.5)

    all_residues = np.unique(np.concatenate([res1, res2, res3]))
    plt.xlim(all_residues.min(), all_residues.max())

    if len(all_residues) > 50:
        step = max(1, len(all_residues) // 20)
        tick_positions = all_residues[::step]
        plt.xticks(tick_positions.astype(int), rotation=0, fontsize=TICK_FONTSIZE)
    else:
        plt.xticks(all_residues.astype(int), rotation=0, fontsize=TICK_FONTSIZE)

    plt.tick_params(axis='both', which='major', labelsize=TICK_FONTSIZE)

    plt.legend(fontsize=LEGEND_FONTSIZE)

    plt.tight_layout()
    plt.savefig('1.svg', dpi=300)
    plt.show()


