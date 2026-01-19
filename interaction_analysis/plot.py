import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

FILE_PATH = '3diff.csv'

COLOR_PHOS = '#800080'    
COLOR_NONPHOS = '#808080' #


AXIS_LABEL_FONTSIZE = 32   
TICK_FONTSIZE = 26         
TYPE_TITLE_FONTSIZE = 32   


N_SAMPLES_PER_BAR = 5000      
SIG_FONTSIZE = 26             
SIG_LINEWIDTH = 1.2          
SIG_Y_PAD = 0.03             

from math import erf, sqrt

def _norm_cdf(x: float) -> float:
    return 0.5 * (1.0 + erf(x / sqrt(2.0)))

def two_proportion_ztest_pvalue(p1: float, p2: float, n1: int, n2: int) -> float:

    k1 = int(round(p1 * n1))
    k2 = int(round(p2 * n2))
    k1 = max(0, min(n1, k1))
    k2 = max(0, min(n2, k2))
    p_pool = (k1 + k2) / (n1 + n2) if (n1 + n2) > 0 else 0.0
    if p_pool <= 0.0 or p_pool >= 1.0:
        return 1.0 if (k1 / n1) == (k2 / n2) else 0.0
    se = sqrt(p_pool * (1.0 - p_pool) * (1.0 / n1 + 1.0 / n2))
    if se == 0:
        return 1.0 if (k1 / n1) == (k2 / n2) else 0.0
    z = ((k1 / n1) - (k2 / n2)) / se
    pval = 2.0 * (1.0 - _norm_cdf(abs(z)))
    return max(0.0, min(1.0, pval))

def p_to_star(p: float) -> str:
    if p < 1e-4:
        return '****'
    if p < 1e-3:
        return '***'
    if p < 1e-2:
        return '**'
    if p < 5e-2:
        return '*'
    return 'ns'

MAX_INTERACTIONS_PER_TYPE = 5  
MAX_PLOTTED_INTERACTIONS = 15  


try:
    df = pd.read_csv(
        FILE_PATH,
        header=None,
        names=['Interaction_String', 'Type', 'Occ_Phospho', 'Occ_NonPhospho', 'Abs_Diff']
    )
except FileNotFoundError:
    print(f"Error: File not found at: {FILE_PATH}")
    exit()

df['Diff'] = df['Occ_Phospho'] - df['Occ_NonPhospho']

TARGET_TYPES = ['hbond', 'salt_bridge', 'hydrophobic']
df_list_padded = []

for t in TARGET_TYPES:
    df_type = df[df['Type'] == t].copy()
    df_top_n = df_type.nlargest(MAX_INTERACTIONS_PER_TYPE, 'Abs_Diff')

    df_non_zero_phos_sub = df_top_n[df_top_n['Occ_Phospho'] > 0].sort_values(by='Occ_Phospho', ascending=False)
    df_zero_phos_sub = df_top_n[df_top_n['Occ_Phospho'] == 0].sort_values(by='Occ_NonPhospho', ascending=True)
    df_sorted_sub = pd.concat([df_non_zero_phos_sub, df_zero_phos_sub])

    count = len(df_sorted_sub)

    if count < MAX_INTERACTIONS_PER_TYPE:
        n_padding = MAX_INTERACTIONS_PER_TYPE - count
        padding_data = {
            'Interaction_String': [np.nan] * n_padding,
            'Type': [t] * n_padding,
            'Occ_Phospho': [0.0] * n_padding,
            'Occ_NonPhospho': [0.0] * n_padding,
            'Abs_Diff': [0.0] * n_padding,
            'Diff': [0.0] * n_padding,
        }
        df_padding = pd.DataFrame(padding_data)
        df_padded_sub = pd.concat([df_sorted_sub, df_padding], ignore_index=True)
    else:
        df_padded_sub = df_sorted_sub

    df_list_padded.append(df_padded_sub)

df_top_final = pd.concat(df_list_padded, ignore_index=True)

bar_group_width = 0.8
bar_width = bar_group_width / 2
ind = np.arange(MAX_PLOTTED_INTERACTIONS)

phos_bar_x = ind - bar_width / 2
nonphos_bar_x = ind + bar_width / 2

colors_heatmap_data = [
    (0.0, 'blue'),
    (0.5, 'white'),
    (1.0, 'red')
]
cmap_heatmap = LinearSegmentedColormap.from_list("BlueWhiteRed_Diff", colors_heatmap_data)
vmin_heatmap = -1.0
vmax_heatmap = 1.0

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

ax_bar.bar(
    phos_bar_x, df_top_final['Occ_Phospho'], bar_width,
    color=COLOR_PHOS, edgecolor='black', label='Phosphorylated'
)
ax_bar.bar(
    nonphos_bar_x, df_top_final['Occ_NonPhospho'], bar_width,
    color=COLOR_NONPHOS, edgecolor='black', alpha=0.8, label='Non-phosphorylated'
)

for i in range(len(df_top_final)):

    if pd.isna(df_top_final.loc[i, 'Interaction_String']):
        continue

    p1 = float(df_top_final.loc[i, 'Occ_Phospho'])
    p2 = float(df_top_final.loc[i, 'Occ_NonPhospho'])
    pval = two_proportion_ztest_pvalue(p1, p2, N_SAMPLES_PER_BAR, N_SAMPLES_PER_BAR)
    star = p_to_star(pval)

    x1 = phos_bar_x[i]
    x2 = nonphos_bar_x[i]
    y = max(p1, p2) + SIG_Y_PAD
    h = SIG_Y_PAD * 0.6
    ax_bar.plot([x1, x1, x2, x2], [y, y + h, y + h, y], color='black', linewidth=SIG_LINEWIDTH)
    ax_bar.text((x1 + x2) / 2.0, y + h + SIG_Y_PAD * 0.2, star,
                ha='center', va='bottom', fontsize=SIG_FONTSIZE)


ax_bar.set_ylabel('Occupancy', fontsize=AXIS_LABEL_FONTSIZE)
ax_bar.set_ylim(0, 1.2)
ax_bar.tick_params(axis='y', labelsize=TICK_FONTSIZE)
ax_bar.tick_params(axis='x', length=0)
plt.setp(ax_bar.get_xticklabels(), visible=False)


current_x = 0
for t in TARGET_TYPES:
    center_x = current_x + 2  

    if current_x > 0:
        ax_bar.axvline(x=current_x - bar_group_width / 2, color='k', linestyle='--', linewidth=1)

    ax_bar.text(
        center_x, ax_bar.get_ylim()[1] * 0.95,
        t.replace('_', ' ').title(),
        ha='center', va='top',
        fontsize=TYPE_TITLE_FONTSIZE,
        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none')
    )

    current_x += MAX_INTERACTIONS_PER_TYPE


heatmap_data = df_top_final['Diff'].values.reshape(1, -1)
im = ax_heatmap.imshow(
    heatmap_data,
    cmap=cmap_heatmap,
    aspect="auto",
    vmin=vmin_heatmap, vmax=vmax_heatmap,
    extent=[-bar_group_width / 2, MAX_PLOTTED_INTERACTIONS - bar_group_width / 2, 0, 1]
)


x_labels = df_top_final['Interaction_String'].fillna(' ').str.replace(';', ' - ')
ax_heatmap.set_xticks(ind)
ax_heatmap.set_xticklabels(x_labels, rotation=60, ha='right', fontsize=TICK_FONTSIZE)

ax_heatmap.set_xlim(-bar_group_width / 2, MAX_PLOTTED_INTERACTIONS - bar_group_width / 2)
ax_heatmap.tick_params(axis='x', length=0, labelsize=TICK_FONTSIZE)


ax_heatmap.set_yticks([])
ax_heatmap.set_ylabel('') 

cbar = fig.colorbar(im, cax=ax_cbar, orientation='vertical')
cbar.ax.tick_params(labelsize=TICK_FONTSIZE)
cbar.set_label('Difference', fontsize=AXIS_LABEL_FONTSIZE)

plt.tight_layout()
plt.savefig('3.svg', dpi=300, bbox_inches='tight')
# plt.show()
