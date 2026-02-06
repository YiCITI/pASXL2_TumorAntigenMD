# pASXL2_TumorAntigenMD
The detailed scripts we used in the article "The phosphorylation of an antigen peptide from ASXL2 protein alters the peptide-HLA binding affinity and the interaction dynamics of the pHLA complex".

## ASXL2 analysis
ASXL2_analysis: script for differential expression and pathway enrichment analysis. 
### ASXL2_Code.R
Analyze the bioinformatic information for ASXL2 and its phosphorylation. Processed by Lin Lv.

## MD
MD: Amber scripts we used in this project.
### make.leap
Make the top file and inpcrd file.  Processed by Jiahui Zhang.
### fix_opt.in
Optimizing the coordinates with a fixed protein complex.  Processed by Jiahui Zhang.
### move_opt.in
Optimizing the coordinates.  Processed by Jiahui Zhang.
### warp_up.in
Warm up the system.  Processed by Jiahui Zhang.
### nvt.in
An NVT run.  Processed by Jiahui Zhang.
### v_opt.in
A CPU NPT run to optimize the volume of the box.  Processed by Jiahui Zhang.
### npt.in
Production run. 



## Test
test: Force field test, RSMD test, and RMSF test. The first part of the results. 
### ff_test_md.py
The script to run MD for testing the force field. Processed by Jiahui Zhang.
### RMSF_plot.R
Plot the RMSF statistical box plot for the force field test. Processed by Lin Lv.
### rmsd.py
Analysis RMSD and plot RMSD figures. Processed by Jiahui Zhang.
### rmsf_csv.py
Extract RMSF information into a CSV file. Processed by Jiahui Zhang.
### rmsf_plot.py
Plot the RMSF figure. Processed by Jiahui Zhang.


## Interaction analysis
interaction_analysis: the non-bonding interaction analysis of the pHLA complex. The second part of the results.
### interaction.py
Extract non-bonding interactions from the trajectories.  Processed by Jiahui Zhang.
### plot.py
Plot the figure. Processed by Jiahui Zhang.


## FE analysis
FE_analysis: the GB and PB scripts for this project. The third part of the results.
### make_top.sh
Extract ligand and complex topology. Processed by Jiahui Zhang/
### mmpbsa_gb.in
Process GB calculation.Processed by Jiahui Zhang.
### mmpbsa_pb.in
Process PB calculation. Processed by Jiahui Zhang
### plot.R
Plot the free energy change figure. Processed by Lin Lv.

## PCA analysis
PCA_analysis: the PCA analysis scripts for this project. The fourth part of the results. 
### pca.py
Extract the PCs from trajectories. Processed by Jiahui Zhang.
### plot.py
Plot the PC map. Processed by Jiahui Zhang.

