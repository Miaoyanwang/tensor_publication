%%% Brain tissue analysis %%%
%%%  Jiaxin Hu 03/07/21   %%%

% lastest update 03/11/21

% Here is the pipeline to use the Joint_tucker_v6.m function.
% Remember to revise the paths!

userpath("/Users/March/Desktop/myresearch/graphical_lasso/data_analysis/multi_layer");
userpath(func_path); % REVISED!
% func_path should be the path involves the main function and the logdet.m

% Read ss data 
load("/Users/March/Desktop/myresearch/graphical_lasso/data_analysis/multi_layer/test1/varvar_withname/ss_data.mat");
gss = struct2cell(ss_data);  

% read n_vector(sample size vector)
gnn = readtable("/Users/March/Desktop/myresearch/graphical_lasso/data_analysis/multi_layer/gene_cov/gtex_brain/pre_gene/nvector.csv");% REVISED!
% nvec_path should be the path of the sample size vector
% nvec_path should be different with ss_data_path, 
% unless sample size vector is not stored as .csv
gnn = gnn{:,:};
gnn(:,1) = [];
gnn = gnn';

% if we use balanced sample size, try
% gnn = 193*ones(1,13); 


% step up parameters
option.initial = 1; % HOSVD, relatively quick
option.Niter = 30;
option.tol =  1;

r = 3; option.rho = 1500;
[Omega,U,Theta0,Theta,convg,rec_obj] = Joint_tucker_v6(gss,gnn,r,option);

% write the results

% result_path should be the path you store the results
result_path = "/Users/March/Desktop/myresearch/graphical_lasso/data_analysis/multi_layer/test1/varvar_withname/result/";
col = cell(1,length(gene_name)+1);
col(:,2:end) = gene_name;
col(:,1) = {'name'};

for i = 1:r % if the rank degenerates, change r to the degenerated rank.
    A = Theta{i};
    A=mat2cell(A,ones(1,size(A,1)),ones(1,size(A,2)));
    A(:,2:end+1)=A;
    A(:,1)=gene_name;
    %A(1,2:end)=gene_name; 
    
    T = cell2table(A);
    T.Properties.VariableNames = col;
    writetable(T, result_path+"Theta_"+i+"_r"+r+"_rho"+option.rho+".csv"); % REVISED!
end

%write Theta0
A = Theta0;
A=mat2cell(A,ones(1,size(A,1)),ones(1,size(A,2)));
A(:,2:end+1)=A(:,1:end);
A(:,1)=gene_name;

T = cell2table(A);
T.Properties.VariableNames = col;
writetable(T, result_path+"Theta0"+"_r"+r+"_rho"+option.rho+".csv"); % REVISED!

%write U
A = U;
A=mat2cell(A,ones(1,size(A,1)),ones(1,size(A,2)));
A(:,2:end+1)=A(:,1:end);
A(:,1)=tissue_name;
T = cell2table(A);
writetable(T, result_path+"U"+"_r"+r+"_rho"+option.rho+".csv"); % REVISED!

csvwrite(result_path+"Obj_curve"+"_r"+r+"_rho"+option.rho+"_iter" + option.Niter+".csv",rec_obj);% REVISED!



