% Code for Figure 6: ME comparison of STD and other methods under model misspecifications

% load software
addpath("software/SupCP-master");

% model misspecification with noniid noise --------------------------------
% load data
rng(10);
load("Figure6_noniid_data.mat");

dup = 30;

err_sup = zeros(2,2,4,dup);
cor_sup = zeros(2,2,4,dup);

final_c_sup = zeros(2,2,4);
final_sd_c_sup = zeros(2,2,4);
final_c_cor_sup = zeros(2,2,4);
final_sd_c_cor_sup = zeros(2,2,4);

rank_range = [2,4,6,8,10];
rank_choose = zeros(2,2,4,dup);
rank_choose_cor = zeros(2,2,4,dup);

args = struct('max_niter',2000,'AnnealIters',1000);
for s = 1:2 % rank
    for r = 1:2 % signal
        for i = 1:4 % cor_level
            for n = 1:dup % dup
                
                disp([para{s,r}{i,n}'])
                Y = y{s,r}{i,n};
                X = tsr{s,r}{i,n};
                X_true = x_true{s,r}{i,n};
                
                test = zeros(1,length(rank_range)); 
                test_cor = zeros(1,length(rank_range));
                
                for k = 1:length(rank_range)
                    rank = rank_range(k);
                    disp(['choose rank as ', num2str(rank)]);
                    
                    [B,V,U,se2,Sf,rec] = SupParafacEM(Y,X,rank,args);
                    
                    x_fit = TensProd([{U}, V]);
                    test(k) = mean((x_fit - X_true).^2, 'all');
                    correlation = corrcoef(x_fit,X_true);
                    test_cor(k) = correlation(1,2);
                end
                
                [M,I] = min(test); 
                [M_cor,I_cor] = max(test_cor);
                rank_choose(s,r,i,n) = rank_range(I);
                rank_choose_cor(s,r,i,n) = rank_range(I_cor);
                
                err_sup(s,r,i,n) = M;
                cor_sup(s,r,i,n) = M_cor;
            end
            
            final_c_sup(s,r,i) = mean(err_sup(s,r,i,:));
            final_sd_c_sup(s,r,i) = std(err_sup(s,r,i,:));
            
            final_c_cor_sup(s,r,i) = mean(cor_sup(s,r,i,:));
            final_sd_c_cor_sup(s,r,i) = std(cor_sup(s,r,i,:));
        end
    end
end

save("mat_output/Figure6_noniid.mat",'err_sup','cor_sup', 'final_c_sup', 'final_sd_c_sup',...
    'final_c_cor_sup','final_sd_c_cor_sup','rank_choose','rank_choose_cor');

