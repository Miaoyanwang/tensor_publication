%%% Brain tissue analysis %%%

% lastest update. Miaoyan Wang 03/13/21

% Here is the pipeline to use the SCENT_tucker1_v6_2.m function.
% Remember to revise the paths!

% userpath("");
% userpath(func_path); % REVISED!
% func_path should be the path involves the main function and the logdet.m

% Read ss data 
load("input/ss_data.mat");
gss = ss_data;


% read n_vector(sample size vector)
gnn = readtable("input/nvector.csv");% REVISED!
% nvec_path should be the path of the sample size vector
% nvec_path should be different with ss_data_path, 
% unless sample size vector is not stored as .csv
gnn = gnn(:,:);
gnn(:,1) = [];
gnn = table2array(gnn);
gnn=gnn';

% if we use balanced sample size, try
% gnn = 193*ones(1,13); 


% step up parameters
option.initial = 1; % HOSVD, relatively quick
option.Niter = 30;
option.tol =  1;
option.fig = 1;
option.TOL = 1;

r = 3; option.rho = 1500;
[Omega,U,Theta0,Theta,convg,rec_obj] = SCENT_tucker1_v6_2(gss,gnn,r,option);

% write the results

gene_name=strrep(gene_name,'.','_');
gene_name=strrep(gene_name,'-','_');


% result_path should be the path you store the results
result_path = "output_r3rho1500/";
col = cell(1,length(gene_name)+1);
col(:,2:end) = gene_name;
col(:,1) = {'name'};

for i = 1:r % if the rank degenerates, change r to the degenerated rank.
    A = Theta{i};
    A=mat2cell(A,ones(1,size(A,1)),ones(1,size(A,2)));
    A(:,2:end+1)=A;


    A(:,1)=gene_name;
    
    T = cell2table(A);
    T.Properties.VariableNames= col;
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



