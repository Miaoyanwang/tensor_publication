% Code for Figure 4: ME comparison among STD, Envelope, GLSNet, mRRR, and SupCP

%load software
addpath("software/SupCP-master");

% Error vs Info Number ----------------------------------------------------
% load data
rng(0);
load("Figure4_mode_data.mat");

dup = 30;

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
for s = 1:2 % signal alpha = 3,6
    for r = 1:2 % core shape = (3,3,3) or (4,5,6)
        for i = 1:3 % info number
            for n = 1:dup % dup
                
                disp([para{s,r}{i,n}'])
                Y = y{s,r}{i,n};
                X = tsr{s,r}{i,n};
                X_true = x_true{s,r}{i,n};
                
                test = zeros(1,length(rank_range)); test_cor = zeros(1,length(rank_range));
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
            
            final_mode_sup(s,r,i) = mean(err_sup(s,r,i,:));
            final_sd_mode_sup(s,r,i) = std(err_sup(s,r,i,:));
            
            final_mode_cor_sup(s,r,i) = mean(cor_sup(s,r,i,:));
            final_sd_mode_cor_sup(s,r,i) = std(cor_sup(s,r,i,:));
        end
    end
end

save("mat_output/Figure4_mode.mat",'err_sup','cor_sup', 'final_mode_sup', 'final_sd_mode_sup',...
    'final_mode_cor_sup','final_sd_mode_cor_sup','rank_choose','rank_choose_cor');

% Error vs Dimension ------------------------------------------------------
rng(2357);

load("Figure4_sample_data.mat");
dup = 30;

err_sup = zeros(2,2,5,dup);
cor_sup = zeros(2,2,5,dup);

final_sample_sup = zeros(2,2,5); 
final_sd_sample_sup = zeros(2,2,5);
final_sample_cor_sup = zeros(2,2,5);
final_sd_sample_cor_sup = zeros(2,2,5);

rank_range = [2,4,6,8,10];
rank_choose = zeros(2,2,5,dup);
rank_choose_cor = zeros(2,2,5,dup);

args = struct('max_niter',2500,'AnnealIters',1000);
for s = 1:2 % signal alpha = 3,6
    for r = 1:2 % core shape = (3,3,3) or (4,5,6)
        for i = 1:5 % dim 
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
            
            final_sample_sup(s,r,i) = mean(err_sup(s,r,i,:));
            final_sd_sample_sup(s,r,i) = std(err_sup(s,r,i,:));
            
            final_sample_cor_sup(s,r,i) = mean(cor_sup(s,r,i,:));
            final_sd_sample_cor_sup(s,r,i) = std(cor_sup(s,r,i,:));
        end
    end
end

save("mat_output/Figure4_sample.mat",'err_sup','cor_sup', 'final_sample_sup', 'final_sd_sample_sup',...
    'final_sample_cor_sup','final_sd_sample_cor_sup','rank_choose','rank_choose_cor');

