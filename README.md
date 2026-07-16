# pASXL2_TumorAntigenMD
The detailed scripts we used in the article "Phosphorylation of a tumor-derived ASXL2 epitope remodels the HLA-bound peptide conformational ensemble and interaction network of the peptide-HLA complex".

## ASXL2 analysis
ASXL2_analysis: script for differential expression and pathway enrichment analysis. 
### ASXL2_Code.R
Analyze the bioinformatic information for ASXL2 and its phosphorylation.


## MD
MD: Amber scripts we used in this project.
### make.leap
Make the top file and inpcrd file.  
### fix_opt.in
Optimizing the coordinates with a fixed protein complex.  
### move_opt.in
Optimizing the coordinates.  
### warp_up.in
Warm up the system. 
### nvt.in
An NVT run.  
### v_opt.in
A CPU NPT run to optimize the volume of the box.  
### npt.in
Production run. 
### npt_pep.in
Peptide-only production run.



## RMS
RMS: root-mean-square analysis. 
### rmsd_plot.py
Analysis RMSD and plot RMSD figures.
### rmsf_plot.py
Analysis RMSF and plot RMSD figures.



## Interaction analysis
interaction_analysis: the non-bonding interaction analysis of the pHLA complex. The second part of the results.
### plot_interactions_2systems_auto_5rep_mean_diff_heatmap_pdf.py
Extract non-bonding interactions from the trajectories.


## FE analysis
FE_analysis: the GB and PB scripts for this project. The third part of the results.
### make_top.sh
Extract ligand and complex topology.
### mmpbsa_gb.in
Process GB calculation.
### mmpbsa_pb.in
Process PB calculation.
### plot_free_energy_phos_prot_wt_5rep_meanline_pvalue.py
Plot the free energy change figure. 

## PCA analysis
PCA_analysis: the PCA analysis scripts for this project. The fourth part of the results. 
### pca_entropy_2systems_auto.py
PCA-derived divergence analysis.
### plot_pca_2systems_auto_5rep.py
Plot the PC map. 
