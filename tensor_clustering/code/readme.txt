###### Software and reproducing code for ``Multiway Spherical Clustering via Degree-Corrected Tensor Block Models"  ########

Jan 23, 2022, Jiaxin

############### Our software: /software/dTBM.tar.gz with manual /software/dTBM_manual.pdf

############### Simulation procedure ############

# Files:

- 8 R scripts to implement the simulations and plot the figures in the main text;

- 1 folder named ``software" stores the softwares of different tensor methods;
- 1 folder named ``presaved" stores the simulation results from prior simulations;
- 1 folder named ``figure" stores the figures with prior simulation results.

** IMPORTANT: The software for tensor-SCORE (Ke et. al. 2019), named ``tensor_score.R", is translated from original Python codes ``Tensor-SCORE-main.zip" provided by the author. **

# Dependencies:

1. R version 4.0.2 or higher.
2. R packages: RSKC, ggplot2, patchwork, RColorBrewer.

# Pipeline:

1. To reproduce the Figure 4, 5, 6, 7, 8, and Table 2 with prior simulation results, run the first part in the scripts Figure4.R, Figure5.R, Figure6.R, Figure7.R, Figure8.R, and Table2.R. 

2. To reproduce the simulations in Figure 4, 5, 6, 7, 8, and Table 2, run the second part in the scripts Figure4.R, Figure5.R, Figure6.R, Figure7.R, Figure8.R, and Table2.R. ** Note that the results for HOSVD, HOSVD+, SCORE may be slightly different with prior simulations due to the randomness of kmeans(). **

3. To reproduce the Peru data analysis in Table 3, run the script Table3.R.

4. To reproduce the HCP data analysis in Figure 9, 10, and 11, run the script Figure9_10_11.R. Particularly, Figure 9 and 11 are plotted by an external software, BrainNetViewer, which is not included here. Necessary files for plotting Figure 9 and 11 are saved in the folder ``presaved/for_Figure9_11".