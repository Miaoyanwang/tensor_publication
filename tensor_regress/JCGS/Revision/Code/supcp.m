% supCP vs STD comparison

addpath("function/SupCP-master");

% PMSE/Cor vs info number
rng(0);
%load("/Volumes/GoogleDrive/My Drive/research/tensor_regress/jcgs_revision/jcgs_code/test/s_2_r_2_i_5_n_30_fig4_sample.mat");
load("mat_data/mode_fig5_v1.mat");

dup = 30;
% store the errors
err_sup = zeros(2,2,3,dup);
cor_sup = zeros(2,2,3,dup);

final_mode_sup = zeros(2,2,3);
final_sd_mode_sup = zeros(2,2,3);
final_mode_cor_sup = zeros(2,2,3);
final_sd_mode_cor_sup = zeros(2,2,3);

rank_range = [2,4,6,8,10];
rank_choose = zeros(2,2,3,dup);
rank_choose_cor = zeros(2,2,3,dup);

args = struct('max_niter',2000,'AnnealIters',1000);
for s=1:2 % signal alpha = 3,6
    for r=1:2 % core shape = (3,3,3) or (4,5,6)
        for i=1:3 % info number
            for n=1:dup % dup
                
                %s = 2; r = 2; i = 3; n = 30;
                disp([para{s,r}{i,n}'])
                Y = y{s,r}{i,n};
                X = tsr{s,r}{i,n};
                X_true = x_true{s,r}{i,n};
                
                test = zeros(1,length(rank_range)); test_cor = zeros(1,length(rank_range));
                for(k=1:length(rank_range))
                    rank = rank_range(k);
                    disp(['choose rank as ', num2str(rank)]);
                    
                    [B,V,U,se2,Sf,rec]=SupParafacEM(Y,X,rank,args);
                    
                    % calculate fitted value
                    x_fit = TensProd([{U}, V]);
                    test(k) = mean((x_fit - X_true).^2, 'all');
                    correlation = corrcoef(x_fit,X_true);
                    test_cor(k) = correlation(1,2);
                end % rank choose
                
                [M,I] = min(test); [M_cor,I_cor] = max(test_cor);
                rank_choose(s,r,i,n) = rank_range(I);
                rank_choose_cor(s,r,i,n) = rank_range(I_cor);
                
                err_sup(s,r,i,n) = M;
                cor_sup(s,r,i,n) = M_cor;
            end % end dup
            
            final_mode_sup(s,r,i) = mean(err_sup(s,r,i,:));
            final_sd_mode_sup(s,r,i) = std(err_sup(s,r,i,:));
            
            final_mode_cor_sup(s,r,i) = mean(cor_sup(s,r,i,:));
            final_sd_mode_cor_sup(s,r,i) = std(cor_sup(s,r,i,:));
        end
    end
end


save("mat_output/fig5_mode.mat",'err_sup','cor_sup', 'final_mode_sup', 'final_sd_mode_sup',...
    'final_mode_cor_sup','final_sd_mode_cor_sup','rank_choose','rank_choose_cor');

% PMSE/Cor vs dim
rng(2357);

load("mat_data/sample_fig5_v1.mat");
dup = 30;
% store the errors
err_sup = zeros(2,2,5,dup);% signal, rank, dim, dup
cor_sup = zeros(2,2,5,dup);

final_sample_sup = zeros(2,2,5); % signal, rank, dim
final_sd_sample_sup = zeros(2,2,5);
final_sample_cor_sup = zeros(2,2,5);
final_sd_sample_cor_sup = zeros(2,2,5);

rank_range = [2,4,6,8,10];
rank_choose = zeros(2,2,5,dup);% signal, rank, dim, dup
rank_choose_cor = zeros(2,2,5,dup);

args = struct('max_niter',2500,'AnnealIters',1000);
for(s=1:2) % signal alpha = 3,6
    for(r=1:2) % core shape = (3,3,3) or (4,5,6)
        for(i=1:5) % dim 
            for(n=1:dup) % dup
                
                %s = 2; r = 2; i = 5; n = 30;
                disp([para{s,r}{i,n}'])
                Y = y{s,r}{i,n};
                X = tsr{s,r}{i,n};
                X_true = x_true{s,r}{i,n};
                
                test = zeros(1,length(rank_range)); test_cor = zeros(1,length(rank_range));
                for(k=1:length(rank_range))
                    rank = rank_range(k);
                    disp(['choose rank as ', num2str(rank)]);
                    
                    [B,V,U,se2,Sf,rec]=SupParafacEM(Y,X,rank,args);
                    
                    % calculate fitted value
                    x_fit = TensProd([{U}, V]);
                    test(k) = mean((x_fit - X_true).^2, 'all');
                    correlation = corrcoef(x_fit,X_true);
                    test_cor(k) = correlation(1,2);
                end % rank choose
                
                [M,I] = min(test); [M_cor,I_cor] = max(test_cor);
                rank_choose(s,r,i,n) = rank_range(I);
                rank_choose_cor(s,r,i,n) = rank_range(I_cor);
                
                err_sup(s,r,i,n) = M;
                cor_sup(s,r,i,n) = M_cor;
            end % end dup
            
            final_sample_sup(s,r,i) = mean(err_sup(s,r,i,:));
            final_sd_sample_sup(s,r,i) = std(err_sup(s,r,i,:));
            
            final_sample_cor_sup(s,r,i) = mean(cor_sup(s,r,i,:));
            final_sd_sample_cor_sup(s,r,i) = std(cor_sup(s,r,i,:));
        end
    end
end

save("mat_output/fig5_sample.mat",'err_sup','cor_sup', 'final_sample_sup', 'final_sd_sample_sup',...
    'final_sample_cor_sup','final_sd_sample_cor_sup','rank_choose','rank_choose_cor');

%load("/Volumes/GoogleDrive/My Drive/research/tensor_regress/jcgs_revision/code/test/i_1_j_1_k_1_m_3_fig6.mat")

% robustness to non-i.i.d.
rng(0);
%load("/Volumes/GoogleDrive/My Drive/research/tensor_regress/jcgs_revision/jcgs_code/test/s_2_r_2_i_5_n_30_fig4_sample.mat");
load("mat_data/fig6.mat");

dup = 30;
% store the errors
err_sup = zeros(2,2,3,dup);
cor_sup = zeros(2,2,3,dup);

final_c_sup = zeros(2,2,3);
final_sd_c_sup = zeros(2,2,3);
final_c_cor_sup = zeros(2,2,3);
final_sd_c_cor_sup = zeros(2,2,3);

rank_range = [2,4,6,8,10];
rank_choose = zeros(2,2,3,dup);
rank_choose_cor = zeros(2,2,3,dup);

args = struct('max_niter',2000,'AnnealIters',1000);
for s=1:2 % rank
    for r=1:2 % signal
        for i=1:3 % cor_level
            for n=1:dup % dup
                
                %s = 2; r = 2; i = 2; n = 30;
                disp([para{s,r}{i,n}'])
                Y = y{s,r}{i,n};
                X = tsr{s,r}{i,n};
                X_true = x_true{s,r}{i,n};
                
                test = zeros(1,length(rank_range)); test_cor = zeros(1,length(rank_range));
                for(k=1:length(rank_range))
                    rank = rank_range(k);
                    disp(['choose rank as ', num2str(rank)]);
                    
                    [B,V,U,se2,Sf,rec]=SupParafacEM(Y,X,rank,args);
                    
                    % calculate fitted value
                    x_fit = TensProd([{U}, V]);
                    test(k) = mean((x_fit - X_true).^2, 'all');
                    correlation = corrcoef(x_fit,X_true);
                    test_cor(k) = correlation(1,2);
                end % rank choose
                
                [M,I] = min(test); [M_cor,I_cor] = max(test_cor);
                rank_choose(s,r,i,n) = rank_range(I);
                rank_choose_cor(s,r,i,n) = rank_range(I_cor);
                
                err_sup(s,r,i,n) = M;
                cor_sup(s,r,i,n) = M_cor;
            end % end dup
            
            final_c_sup(s,r,i) = mean(err_sup(s,r,i,:));
            final_sd_c_sup(s,r,i) = std(err_sup(s,r,i,:));
            
            final_c_cor_sup(s,r,i) = mean(cor_sup(s,r,i,:));
            final_sd_c_cor_sup(s,r,i) = std(cor_sup(s,r,i,:));
        end
    end
end


save("mat_output/fig6.mat",'err_sup','cor_sup', 'final_c_sup', 'final_sd_c_sup',...
    'final_c_cor_sup','final_sd_c_cor_sup','rank_choose','rank_choose_cor');



% sanity check 



% test
%%Simple simulation for supervised parafac (cp) method

Utrue=20*randn(30,2);
Vtrue=randn(10,2);
Vtrue=normc(Vtrue); 
Wtrue=randn(10,2);
Wtrue=normc(Wtrue);
Ztrue=randn(10,2);
Ztrue=normc(Ztrue);
Xtrue=TensProd({Utrue,Vtrue,Wtrue,Ztrue});
X=  Xtrue+5*randn(30,10,10,10); %%add noise
%Y=Utrue(:,1)+randn(30,1); %Y is vector related to one of the signals in U
Y=Utrue(:,2)+randn(30,2); %Y is vector related to one of the signals in U

args = struct('convg_thres',10^(-5),'max_niter',2000,'AnnealIters',1000); 
args = struct('AnnealIters',1000);
% random intialization leads to different convergence results
[B,V,U,se2,Sf,rec]=SupParafacEM(Y,X,2,args);


% calculate fitted value
x_fit = TensProd([{U}, V]);
err = mean(x_fit - Xtrue, 'all');

% load mat

load("test.mat");
load("/Volumes/GoogleDrive/My Drive/research/tensor_regress/jcgs_revision/jcgs_code/test/s_1_r_1_i_1_n_1_fig5_mode.mat");
load("/Volumes/GoogleDrive/My Drive/research/tensor_regress/jcgs_revision/jcgs_code/test/s_2_r_2_i_5_n_29_fig5_sample.mat");

% use r data 
args = struct('max_niter',2000,'AnnealIters',1000); 
[B,V,U,se2,Sf,rec]=SupParafacEM(y,tsr,10,args);

x_fit = TensProd( [{U}, V]);
err = mean((x_fit - x_true).^2, 'all'); % mse

a = corrcoef(x_fit,x_true); % corre

a = [2,4,6];

save("info_output.mat",'a');
