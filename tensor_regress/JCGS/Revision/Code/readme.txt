##### Our software and reproducing codes for ``Generalized tensor decomposition with features on multiple modes" ######

### Our software: /software/tensorregress.tar.gz, also available in https://CRAN.R-project.org/package=tensorregress.

### Simulation Procedure: 

- 9 R scripts and 3 Matlab scripts to implement the simulations and plot the figures in the main text; 

- 1 folder named ``rawdata" stores the data for real data analyses;
- 1 folder named ``software" stores the softwares of different tensor methods;
- 2 folders named ``mat_output" and ``presaved" store the outputs from Matlab scripts and the simulation results from prior simulations, respectively.

** IMPORTANT: The software for GLSNet (Zhang et. al. 2018), named ``netglm_code_adjust.R", is NOT PUBLIC yet, and the version we submit is modified from personal correspondence with the author. **

# Simulation Dependencies:

1. R version 4.0.2 or higher.
2. Matlab version R2020b or higher. 
3. R packages: rTensor, TRES, rrpack, rmatio, MASS, pracma, speedglm, lattice, ggplot2, patchwork, RColorBrewer.

** Note that ``rTensor" and ``TRES" are not available in RCRAN currently. Try install_version() in package ``remotes" to install the previous version of these two packages.  **

# Simulation Pipeline #

1. To reproduce the Figures 3, 4, 5, 6, and Table 3 with prior simulation results, run the first part in scripts Figure3.R, Figure4.R, Figure5.R, and Table3.R.

2. To reproduce the simulations and analyses in Figure 2, 3, 8, and Table 3, run the second part in scripts Figure2.R, Figure3.R, Figure8.R, and Table3.R.

3. To reproduce the simulations in Figure 4, 5, and 6, run the second part of scripts Figure4.R, Figure5.R, and Figure6.R to obtain the results for the methods implemented in R and the input data for the method implemented in Matlab. Then, run the scripts Figure4.m, Figure5.m, Figure6.m to obtain the output from Matlab. Last, go back to R and generate the figures with the first part in Figure4.R, Figure5.R, Figure6.R.

4. To reproduce the analysis in Figure 7, run the script Figure7.R. Due to the high-dimensionality of the data, we do not suggest to run on a personal laptop. Figure 7 is plotted by an external software, which is not included here. The matrix for plotting the connectivity is saved in the folder presaved/data_for_Figure7.