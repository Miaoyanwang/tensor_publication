%%% Brain tissue analysis %%%
%%%  Jiaxin Hu 03/01/21   %%%

% Here is the pipeline to use the Joint_tucker_v6.m function.
% Remember to revise the paths!


userpath(func_path); % REVISED!
% func_path should be the path involves the main function and the logdet.m

% Read ss data 
K = 13;
gSS = cell(K,1);

files = dir('ss_data_path/*.csv'); % REVISED!
% ss_data_path should be the path of the folder ONLY includes ss_data

for k=1:K
    name = strcat('ss_data_path/', files(k).name);% REVISED!
    read = readtable(name, 'ReadRowNames',true);
    read = read{:,:};
    gSS{k} = read;
end

% read n_vector(sample size vector)
gnn = readtable("nvec_path/nvector.csv");% REVISED!
% nvec_path should be the path of the sample size vector
% nvec_path should be different with ss_data_path, 
% unless sample size vector is not stored as .csv
gnn = gnn{:,:};
gnn(:,1) = [];
gnn = gnn';

% if we use balanced sample size, try
% gnn = 193*ones(1,K); 


% step up parameters
option.initial = 1; % HOSVD, relatively quick
option.Niter = 30;
option.tol =  1;

r = 3; option.rho = 1000;
[Omega,U,Theta0,Theta,convg,rec_obj] = Joint_tucker_v6(gSS,gnn,r,option);

% write the results

% result_path should be the path you store the results

for i = 1:r % if the rank degenerates, change r to the degenerated rank.
    csvwrite("result_path/Theta_"+i+"_r"+r+"_rho"+option.rho+".csv",Theta{i});% REVISED!
end

csvwrite("result_path/U"+"_r"+r+"_rho"+option.rho+".csv",U);% REVISED!
csvwrite("result_path/Theta0"+"_r"+r+"_rho"+option.rho+".csv",Theta0);% REVISED!
csvwrite("result_path/Obj_curve"+"_r"+r+"_rho"+option.rho+"_iter" + option.Niter+".csv",rec_obj);% REVISED!



