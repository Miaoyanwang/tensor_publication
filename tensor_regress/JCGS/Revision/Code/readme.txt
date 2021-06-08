# Revision of ``Supervised Tensor Decomposition with Side Information"

# Create: 06/08/2021 Last edition: 06/08/2021 (Jiaxin)

- Figures: simulation Figures 4,5,6,7 in the main paper.
- function: useful functions and softwares required by the simulations.
- mat_output: output by implementing SupCP method in MatLab.
- presaved: simulation results for plotting.
- Rscripts: run the simulations corresponding to Figures 4,5,6,7.
- supcp.m: implement SupCP method in Matlab.


### Outline ###

To reproduce the Figure 4,5,6 in the main text, follow the steps below:

1. Run the second part of Figure4.R, Figure5.R, Figure5.R to get the input for ssupcp.m.
2. Run supcp.m to get the results of SupCP in mat_output.
3. Read SupCP results in mat_output and plot the simulation figure with the first part of Figure4.R, Figure5.R, Figure5.R.