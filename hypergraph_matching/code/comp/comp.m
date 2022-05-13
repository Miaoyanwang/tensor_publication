%%%%% distance statistics comparison

clear all;
addpath './degree_profile-master';
addpath './degree_profile-master/Slashdot';

sigma_vec= [0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6];
n = 50;

dup = 10; 
r_mat = zeros(length(sigma_vec), 2); % sigma, 2 methods
r_mat_sd = zeros(length(sigma_vec), 2);

for ind_s = 1:1:length(sigma_vec)

    r_WASS_vec = zeros(dup,1);
    r_TV_vec = zeros(dup,1);
    
    sigma = sigma_vec(ind_s);

    for ind_dup = 1:1:dup

        disp(['sigma = ', num2str(sigma), ', duplication = ', num2str(ind_dup)]);

        seed = 100*ind_s + 10*ind_dup;
        rng(seed);

        % generate data

        p=1/2;
        %parent graph
        A=binornd(1,p,n,n);
        A=triu(A,1);
        A=A+A';
        
        s=1-sigma^2;      %sigma^2 = delta = 1-s
        %subsample graph
        Z1=binornd(1,s,n,n);
        Z1=triu(Z1,1);
        Z1=Z1+Z1';
        A1=A.*Z1; 
        
        perp_rnd=randperm(n);
        %perp_rnd=[1:1:n];
        S_rnd=zeros(n,n);
        S_rnd(1:1:n,perp_rnd)=eye(n);
        A_permuted=S_rnd*A*S_rnd';
        
        Z2=binornd(1,s,n,n);
        Z2=triu(Z2,1);
        Z2=Z2+Z2';
        A2=A_permuted.*Z2;
        
        W1=A1;
        W2=A2;

        % calculate the distributions

        deg1=sum(W1);       % degree sequence % colSum
        deg2=sum(W2);
        [deg1_sort,ind1]=sort(deg1); % increasing 
        [deg2_sort,ind2]=sort(deg2);

         % degree profile (sorted)        
         N_deg1=cell(n,1);
         F_deg1=cell(n,1);
         N_deg2=N_deg1;
         F_deg2=F_deg1;
         for ind=1:1:n  
            %ind = 1;
            temp1=deg1_sort(logical(W1(ind,ind1)));
            temp2=deg2_sort(logical(W2(ind,ind2))); 
            [temp1_unique,temp1_a,temp1_c]= unique(temp1);
            N_deg1{ind}=temp1_unique;
            temp1_counts = accumarray(temp1_c,1);
            F_deg1{ind}=temp1_counts/deg1(ind);

            [temp2_unique,temp2_a,temp2_c]= unique(temp2);
            N_deg2{ind}=temp2_unique;
            temp2_counts = accumarray(temp2_c,1);
            F_deg2{ind}=temp2_counts/deg2(ind); 
         end 


         % empirical W1

         D=n*ones(n,n);
                for ind_i=1:1:n
                  if deg1(ind_i)>0
                    for ind_j=1:1:n 
                        if deg2(ind_j)>0    
                        D(ind_i,ind_j) =  dwass_discrete2(N_deg1{ind_i},N_deg2{ind_j},F_deg1{ind_i},F_deg2{ind_j});  
                            %Wasserstein distance
                        end
                    end
                   end
                end

        S_WASS= greedy_match((-1)*D');
        r_WASS_vec(ind_dup)=full(sum(dot(S_rnd,S_WASS))/n);
        
        % smooth tv 
        L_seq = 2:2:n;
        r_tv_seq = zeros(length(L_seq),1);

        for int_L = 1:1:length(L_seq)
            L = L_seq(int_L);
            disp(['smoothed TV norm, smooth parameter L = ', num2str(L)]);
            D_tv=n*ones(n,n);
            for ind_i=1:1:n
              if deg1(ind_i)>0
                for ind_j=1:1:n 
                    if deg2(ind_j)>0    
                    D_tv(ind_i,ind_j) =  smooth_tv_b(N_deg1{ind_i},N_deg2{ind_j},F_deg1{ind_i},F_deg2{ind_j},n,p,L); 
                    end
                end
               end
            end

            S_tv= greedy_match((-1)*D_tv');
            r_tv_seq(int_L)=full(sum(dot(S_rnd,S_tv))/n);
        end

        r_TV_vec(ind_dup) = max(r_tv_seq);

    end % end ind_dup

    r_mat(ind_s, 1) = mean(r_WASS_vec);
    r_mat(ind_s, 2) = mean(r_TV_vec);

    r_mat_sd(ind_s, 1) = std(r_WASS_vec);
    r_mat_sd(ind_s, 2) = std(r_TV_vec);

end % end ind_s


%%% plotting
line_width=1.5;
Marker_size=6;
plot_spec={'k-+','m-*'};
leng_spec={'Z dist','W1 dist'};
i=1;
figure;

errorbar(sigma_vec, r_mat(:,2), r_mat_sd(:,2),plot_spec{1},'LineWidth', line_width, 'MarkerSize', Marker_size );
hold on;
errorbar(sigma_vec, r_mat(:,1), r_mat_sd(:,1),plot_spec{2},'LineWidth', line_width, 'MarkerSize', Marker_size );


legend(leng_spec,'location', 'best', 'FontSize', 20,'Interpreter','latex');
xlabel('$\sigma$','FontSize',20,'Interpreter','latex');
ylabel ('fraction of correctly matched pairs','FontSize',20,'Interpreter','latex');

