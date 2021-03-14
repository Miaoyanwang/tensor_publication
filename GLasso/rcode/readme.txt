# Brain tissue analysis
# Origin: 03/07/21 Lastest Update: 03/12/21

### Scripts ###
To go through the brain tissue analysis, five scripts are needed:

pipeline.m : Matlab code to obtain the results using SCENT method; 

###  Input/Dependency ### 
Algorithm: SCENT_tucker1_v6_2.m, logdet.m
ss_data.mat: a collection of sample covariance matrices for 13 brain tissues 
nvector.csv: a csv file recording the sample size for each brain tissue.  

Output: U.csv(membership matrix), Theta0.csv(global connection), Theta_l.csv(l = 1,...,r, component in the prediction matrix), Obj_curve.csv(trajectory of objective value); 

### Tuning parameters ###

The default tuning parameters in pipeline.m are:
 - r = 3
 - rho (penalty parameter) = 1500
 - option.Niter (max iteration) = 50
 - option.tol (stopping threshold) = 1
 - option.initial (initialization) = 1 (HOSVD)

The default setting is the best setting for ss_data.mat is r = 3, rho = 1500. 

### Results ###

Results obtained by the algorithm are stored in folder results, where
 
 - r3rho1500: has the results with ss_data.mat under r = 3, rho = 1500.

###  plot ### 
barplot.R: the R code to make the vertical bar plots in network.pdf. 

## instruction for plotting network?