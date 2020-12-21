
% following the setup in Guo et al (2011)
% Focus on:
% 1. effect of n, p, K, r
% 2. tissue clustering (vs clustering methods)
% 3. network recovery (vs separate graphical model and some joint graphical models)
% 3'. rank-1 case (ROC with varying tuning parameter) 


%% Basic Setup
% simulate low-rank precision tensor
K=10; % num of layers
p=100; % dimension
r=3; % number of low-rank structure
n=100*ones(1,K); % sample size per tissue


% common network: none
Theta0=eye(p);
% generate rank-3 group networks
% chain network
temp_s=cumsum(rand(1,p)*0.5+0.5);
temp_Sigma1=zeros(p);
for i=1:p
    for j=1:p
        temp_Sigma1(i,j)=exp(-abs(temp_s(i)-temp_s(j))/2);
    end
end
temp_Omega1=2*inv(temp_Sigma1); 
temp1=temp_Omega1-diag(diag(temp_Omega1));  % remove diagonal
[~,Theta1]=cov2corr(temp1+diag(sum(abs(temp1),2))); % add diagonal to make psd, then covert to diag=1
Theta1=Theta1-eye(p); % no diagonal
a=diag(Theta1,1);
Theta1=diag(a,1)+diag(a,1)';
heatmap(Theta1);


% nearest-neighbor network
tempx=rand(p);
tempy=rand(p);
dist=inf(p);
temp_Omega2=zeros(p);
for i=1:p
    for j=(i+1):p
        dist(i,j)=sqrt((tempx(i)-tempx(j))^2+(tempy(i)-tempy(j))^2);
    end
end
dist=triu(dist)+triu(dist)';
for i=1:p
    [~,ind]=sort(dist(i,:));
    ind=ind(1:3);% 3-nearest neighbor
    index=ind(ind>i); % to avoid repeatence
    for j=index  
        temp_rn=sign(randn(1))*(rand(1)*0.5+0.5);
        temp_Omega2(i,j)=temp_rn;
    end
end
temp1=temp_Omega2+temp_Omega2';
[~,Theta2]=cov2corr(temp1+diag(sum(abs(temp1),2))); % add diagonal to make psd, then covert to diag=1
Theta2=Theta2-eye(p); % no diagonal
heatmap(Theta2);


% random network
temp_Omega3=zeros(p);
for i=1:p
    for j=(i+1):p
        temp_rn=sign(randn(1))*(rand(1)*0.5+0.5);
        if rand(1)<0.1 % only 10% edges
            temp_Omega3(i,j)=temp_rn;
        end
    end
end
temp1=temp_Omega3+temp_Omega3';
[~,Theta3]=cov2corr(temp1+diag(sum(abs(temp1),2))); % add diagonal to make psd, then covert to diag=1
Theta3=Theta3-eye(p); % no diagonal
heatmap(Theta3);



% true Theta tensor
trueTheta={Theta1,Theta2,Theta3};
trueTTheta=zeros(r,p,p);
for k=1:r
    trueTTheta(k,:,:)=trueTheta{k};
end


U=normalize([[1,-1,zeros(1,K-2)]',[0,0,1,-1,zeros(1,K-4)]',[0,0,0,0,1,-1,zeros(1,K-6)]'],'norm'); % orthogonal
Omega=cell(1,K); % precision matrix (desired output)
TOmega=zeros(K,p,p);
SS=cell(1,K); % covariance matrix (observed, input, no noise)
for k=1:K
    Omega{k}=Theta0+U(k,1)*Theta1+U(k,2)*Theta2+U(k,3)*Theta3;
    if min(eig(Omega{k}))<0
    error([num2str(k),'th Precision matrix not positive definite!'])
             end
    SS{k}=inv(Omega{k});
    SS{k}=(SS{k}+SS{k}')/2;
    TOmega(k,:,:)=Omega{k};
end

                    
% now we have U, Theta0/Theta1/Theta2/Theta3/trueTheta/trueTTheta, Omega/TOmega, SS
           
%% simulate raw data from SS
sSS=cell(1,K); % sample covariance cell
for k=1:K
    X=mvnrnd(zeros(1,p),SS{k},n(k));
    sSS{k}=cov(X);
end
% % check
figure(1);clf;
subplot(1,2,1)
mesh(SS{3}) % nearest neighbor network
subplot(1,2,2)
mesh(sSS{3})
[~,out]=cov2corr(SS{4});
out(1:15,1:15) % we hope the nonzeros are distinct from 0 (what we call, high signal)
           
           
% noisy case
k=2; %% try k =1, 2, 3,
Omega_glasso = GLasso(sSS{k}, n(k)); 
% check
figure(4);clf;
subplot(1,2,1)
heatmap(double(Omega{k}~=0));
title(['sparsity=',num2str(1-sum(sum(triu(Omega{k})~=0))/(p*(p+1)/2))]);
subplot(1,2,2)
heatmap(double(Omega_glasso~=0))
title(['sparsity=',num2str(1-sum(sum(triu(Omega_glasso)~=0))/(p*(p+1)/2))]);
%  why the estimation is so bad compared to the truth?? can I conclude this is because of matlab GLasso? 
           
           
%% Evaluation of Results
% clustering results
U_est
% network recovery
[FP1,FN1]=falsePN(Theta1,Theta_est{1})
figure(1);clf;
subplot(1,2,1)
heatmap(Theta1);
subplot(1,2,2)
heatmap(Theta_est{1})

[FP2,FN2]=falsePN(Theta2,Theta_est{2})
figure(2);clf;
subplot(1,2,1)
heatmap(Theta2);
subplot(1,2,2)
heatmap(Theta_est{2})

[FP3,FN3]=falsePN(Theta3,Theta_est{3})
figure(3);clf;
subplot(1,2,1)
heatmap(Theta3);
subplot(1,2,2)
heatmap(Theta_est{3})
