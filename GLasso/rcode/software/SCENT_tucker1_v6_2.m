function [Omega,U,Theta0,Theta,convg,rec_obj]=SCENT_tucker1_v6_2(SS,nn,r,option)
% This function solves Simultaneous Clustering and Estimation of Networks
% via sparse Tucker1 decomposition (SCENT). The input is a collection of sample
% covariance matrices on the same variables (e.g., genes) but under 
% different conditions/populations (e.g., tissues, or cancer types); the
% output is two-fold: precision matrices corresponding to different conditions, capturing 
% an undirected graph among variables (i.e. Omega), and multiple clusters of
% conditions (NON-overlapping and maybe nonexhaustive) in U and corresponding
% subgraph structure in Theta.
% 
%
% Model: 
%       C=1\circ Theta_0 + \sum_{l=1}^r u_l\circ Theta_l
% where C is K*p*p precision tensor, each u_l gives one cluster of
% conditions, and each Theta_l give one subgraph defining that cluster.
% Theta_0 is a global precision structure (acting as an intercept)
% right-hand-side is essentially a sparse Tucker1 decomposition.
%
% Optimization:
% Penalized likelihood framework where likelihood is on C and sparsity
% penalty is on u and Theta. We use ADMM.
% 
%
% input: 
%
%   SS     K*1 cell array, each cell is a p*p sample covariance matrix (PSD) 
%
%   nn     K*1 vector, sample sizes for K different populations
% 
%   r      positive integer, desired number of clusters for K populations 
%          although the actual number of clusters may be smaller than r
%
%   option
%          initial  1 = HOSVD initial (default); 
%                   0 = zero initial
%                   2 = random initial
%                   3 = ground truth
%
%          threshold  1=thresholding (default); 
%                     0=no thresholding 
%
%          varyrho  0=fixed rho (default); 1=adaptive rho
%          maxrho   5 (default): max rho. Unused if varyrho==0 (with
%                   initial rho=0.1 and increase of 1.1, it takes 41
%                   iterations to reach maxrho)
%
%          rho      initial step size, default rho=0.1
%
%          Tol      default 1E-4, 
%
%          Niter	default 500
%
%          fig      1 show checking figures; 0 (default) no show
%
% Output: 
%
%   Omega   K*1 cell array, each being a pd precision matrix for one population
%
%   U       K*r matrix, with unit-norm, strictly orthogonal (nonoverlapping) sparse columns, centered,
%           each column gives a cluster of the K populations
%
%   Theta0  p*p matrix, global structure of precision 
%
%   Theta   r*1 cell array, each being a p*p sparse, symmetric matrix,
%           corresponding to one column of U, capturing a subgraph that
%           defines the cluster of populations. 
%
%   convg   an indicator of whether convergence is achieved in the end
%           1=converging, 0=not converging
%
%
%
% updated BIC for more accurate selection for subsequent ranks
% Modified stopping criterion (from primal res<1E-3 to relative primal res<1E-4) 
% 
% Add 1U=0 constraint (3/12/2021)

% default parameters
K=length(SS); % number of populations
p=size(SS{1},1); % number of variables

if isfield(option, 'initial')
    initial=option.initial;
else
    initial=1;
end
if isfield(option, 'threshold')
    threshold=option.threshold;
else
    threshold=1;
end
if isfield(option, 'varyrho')
    varyrho=option.varyrho;
else
    varyrho=0;
end
if isfield(option, 'maxrho')
    maxrho=option.maxrho;
else
    maxrho=5;
end
if isfield(option, 'rho')
    rho=option.rho;
else
    rho=0.1;
end
if isfield(option, 'Tol')
    Tol=option.Tol;
else
    Tol=1E-4;
end
if isfield(option, 'Niter')
    Niter=option.Niter;
else
    Niter=500;
end
if isfield(option, 'fig')
    fig=option.fig;
else
    fig=0;
end






% initial values for Omega, TOmega, Lambda, TLambda; Theta0, TTheta, U
Omega=cell(1,K); % precision matrices
Theta=cell(1,r); % subgraph matrices
Lambda=cell(1,K); % Lagrange multiplier matrices
U=zeros(K,r);
for k=1:K
    if ~issymmetric(SS{k}) 
        error('Input sample covariance matrix is not symmetric!')
    end
    if min(eig(SS{k}))<-1E-5
        error('Input sample covariance matrix is not positive semi-definite!')
    end
    Omega{k}=pinv(SS{k}+(0.1*mean(diag(SS{k})))*eye(p)); % make SS pd, so initial Omega is not too off
    Lambda{k}=zeros(p);
end
TOmega=cell2tensor(Omega); % K*p*p
TLambda=cell2tensor(Lambda); % K*p*p
if initial==1
    % HOSVD initial
    Theta0=permute(mean(TOmega,1),[2,3,1]); % p*p mean
    deTOmega=bsxfun(@minus,matricize1(TOmega),Theta0(:)'); % demeaned Omega tensor -> matrix
    [U,~,~]=svds(deTOmega,r);
    U=normalize(normalize(U,'center'),'norm'); % center and unit norm
    TTheta=product1(U',TOmega-product1(ones(K,1),mean(TOmega,1))); % r*p*p
    if r==1 % TTheta is a p*p matrix
        TTheta=permute(TTheta,[3,1,2]);
    end
elseif initial==0
    % zero initial
    Theta0=zeros(p); 
    U=zeros(K,r);
    TTheta=zeros(r,p,p);
elseif initial==2 
    % random initial
    temp=rand(p);
    Theta0=eye(p)+(temp+temp')/2;
    [U,~,~]=svds(randn(K),r);
    U=normalize(normalize(U,'center'),'norm'); % center and unit norm
    for l=1:r
        temp=randn(p);
        Theta{l}=0.1*(temp+temp'); 
    end
    TTheta= cell2tensor(Theta);
elseif initial==3
    Theta0=option.trueTheta0;
    TTheta=option.trueTTheta;
    U=option.trueU;
end



% initial objective function value
obj_ini=ObjVal(SS,nn,Omega); % without penalty






%%%%%%%%%%%%%%%
% ADMM
niter=0;
diff=inf;
rec_obj=obj_ini; % record obj value
rec_primal=[]; % record total primal residual (relative)
rec_dual=[]; % record total dual residual (relative)
rec_Lambda=[]; % record total norm of Lagrange multipliers
while niter<Niter  && abs(diff)>Tol
    niter=niter+1;
    TOmega_old=TOmega;
    niter
    
    %%%%%%%%%%%%% Primal Steps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % estimate Omega
    for k=1:K
        Mat=-nn(k)*SS{k}+rho*(Theta0+product1(U(k,:),TTheta))-Lambda{k};
        Mat=(Mat+Mat')/2;
        [V,D]=eig(Mat);
        tempd=diag(D);
        tempd_new=-tempd/nn(k); % miaoyan's alternative (more robust???)
        newd_new=2./(tempd_new+sqrt(tempd_new.^2+4*rho/nn(k)));    
        Omega{k}=V*diag(newd_new)*V';
    end
    TOmega=cell2tensor(Omega);
    


    % estimate Theta0 (no penalty)
    Theta0=permute(mean(TOmega+TLambda/rho,1),[2,3,1]);
    
    
    % estimate U (with hard thresholding) and Theta (with soft thresholding) , together (not sequentially)
    % **unit-rank deflation, adaptive BIC (for both tuning in each layer) 
    curr_tens=TOmega+TLambda/rho-product1(ones(K,1),permute(Theta0,[3,1,2]));
    keep=1:K; % nonzero U space
    if threshold==0 % no thresholding
        [U,~,~]=svds(matricize1(curr_tens),r);
        U=normalize(normalize(U,'center'),'norm'); % center and unit norm (this should be redundant since curr_tens has been demeaned)
        rawTheta=product1(U',curr_tens);  
        for l=1:r
            Theta{l}=permute(rawTheta(l,:,:),[2,3,1]);
            Theta{l}=(Theta{l}+Theta{l}')/2;
        end
    else
    for l=1:r
        curr_mat=matricize1(curr_tens(keep,:,:));
        [rawu,~,~]=svds(curr_mat,1);
        rawTheta=product1(rawu',curr_tens(keep,:,:));  
        rawTheta=(rawTheta+rawTheta')/2;

        % thresholding
        minlam=min(abs(rawTheta(:)));
        maxlam=max(abs(rawTheta(:)));
        seg=100; % number of L1 tuning candidates
        candlam=minlam+(1:seg)*(maxlam-minlam)/(seg+1);        
        BIC=inf(length(keep),seg);
        for ind_u=2:length(keep) % at least 2 nonzeros, in order to have mean zero
            for ind_v=1:seg % larger=sparser
                % est u with L0
                [truncated_u,nonzero_index]=thres(rawu,ind_u);
                nonzero_u=truncated_u(nonzero_index);
                curr_u=zeros(size(rawu));
                curr_u(nonzero_index)=normalize(normalize(nonzero_u,'center'),'norm'); % L0 hard thresholding, center and unit norm
                
                % est Theta with L1
                curr_Theta=SoftThres(rawTheta,candlam(ind_v)); % L1 soft thresholding
        
                % BIC
%                 BIC(ind_u,ind_v)=log(norm(curr_mat-curr_u*curr_Theta(:)','fro')^2/numel(curr_mat))+log(numel(curr_mat))*(ind_u*sum(curr_Theta(:)~=0))/numel(curr_mat);

                % modified BIC
                BIC(ind_u,ind_v)=log(norm(curr_mat-curr_u*curr_Theta(:)','fro')^2/numel(curr_mat))+log(K*p^2)*(ind_u*sum(curr_Theta(:)~=0))/(K*p^2); % i.e. keep totalSampleSize unchanged over ranks; or equiv, treating fitted part as zero
            end
        end     
        [~,index]=min(BIC(:));
        [opt_u,opt_v]=ind2sub(size(BIC),index);
        % best u and theta
        [truncated_u,nonzero_index]=thres(rawu,opt_u);
        nonzero_u=truncated_u(nonzero_index);
        curr_u=zeros(size(rawu));
        curr_u(nonzero_index)=normalize(normalize(nonzero_u,'center'),'norm'); % L0 hard thresholding, center and unit norm
        curr_Theta=SoftThres(rawTheta,candlam(opt_v));
            % check
            if fig==1 && l<=4 && size(BIC,1)>1
            figure(99);
            subplot(2,2,l)
%             mesh(BIC);
            plot(BIC(opt_u,:));
            title(['BIC: Layer ',num2str(l),' (',num2str(opt_u),' nonzeros)']);
            xlabel('L1 sparsity level in Theta')
            ylabel('BIC score')
            end
      
        % deflation, get residual tensor for next rank
        U(:,l)=zeros(K,1);
        U(keep,l)=curr_u;
        Theta{l}=curr_Theta;
        curr_tens=curr_tens-product1(U(:,l),permute(Theta{l},[3,1,2]));
        keep=keep(curr_u==0);
        if length(keep)<2 % cannot seek for more layers
            U=U(:,1:l);
            Theta=Theta(1:l); % compress subsequent ranks
            break
        end
    end    
    end
    TTheta=cell2tensor(Theta);
    

%     U
%     mean(U,1)
    

    
   
    %%%%%%%%%%%%% Dual Steps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for k=1:K
        Lambda{k}=Lambda{k}+rho*(Omega{k}-Theta0-product1(U(k,:),TTheta));
    end
    TLambda=cell2tensor(Lambda);
    rec_Lambda=[rec_Lambda,norm(TLambda(:),'fro')];
    
   

    
    
    %%%%%%%%%%%%% Stopping Rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % primal and dual residuals
    res=TOmega-product1(ones(K,1),permute(Theta0,[3,1,2]))-product1(U,TTheta);
    primal=norm(res(:),'fro');
    relprimal=primal/norm(TOmega(:),'fro');
    rec_primal=[rec_primal,relprimal];
    dual=norm(TOmega(:)-TOmega_old(:),'fro'); 
    reldual=dual/norm(TOmega(:),'fro');
    rec_dual=[rec_dual,reldual];  
    
    % update rho
    if varyrho && niter>100 % only start varying rho after 100 iterations (to expedite the convergence)
        if primal>5*dual % if primal residual too large, increase rho to make Omega similar to U*Theta
            rho=min(maxrho,1.1*rho); % steadily increasing rho
        end
    end
    
    % objective function value
    obj=ObjVal(SS,nn,Omega);  % just calc by Omega, irrelevant to U and Theta
    obj
    rec_obj=[rec_obj,obj];
    
    % stopping rule
    diff=relprimal; % relative diff
%     diff=reldual;
%     diff=rec_obj(1,end-1)-rec_obj(1,end);


    % Check Figures
    if fig==1
        % obj fcn values
        figure(100);clf; 
        subplot(2,2,1)
        plot(0:niter,rec_obj,'bo-');
        title('-loglik (w/o penalty)');
        subplot(2,2,2)
        plot(1:niter,rec_Lambda,'o-'); 
        title('Norm of Multipliers');
        % primal and dual residuals
        subplot(2,2,3)
        plot(1:niter,rec_primal,'o-');
        title(['Relative Primal: ',num2str(primal)]);
        subplot(2,2,4)
        plot(1:niter,rec_dual,'o-');
        title(['Relative Dual: ',num2str(dual)]);
        drawnow
    end


end


if niter==Niter
    disp(['SCENT does NOT converge after ',num2str(Niter),' iterations!']);
    convg=0;
else
    disp(['SCENT converges after ',num2str(niter),' iterations.']);      
    convg=1;
end
 
end



function out=SoftThres(Theta,lam)
% this function soft thresholds every entry in Theta by lambda
pos=max(Theta-lam,0);
neg=min(Theta+lam,0);
out=pos+neg;
end


function [out,nonzero_index]=thres(u,lam)
% this function truncate u to only keep the largest lam absolute values
% lam: # of non-zero values
% nonzero_index: index vector of length lam
[~,ind]=sort(abs(u),'descend');
nonzero_index=ind(1:lam);
out=u;
out(ind((lam+1):end))=0;
end

function obj=ObjVal(SS,nn,Omega)
% Calc the negative log likelihood function  
obj=0;
K=length(nn);
for k=1:K
    obj=obj+nn(k)*(trace(SS{k}*Omega{k})-logdet(Omega{k},'chol'));
end
end


function T=cell2tensor(C)
% Convert a length-K cell array of p*q matrices to a K*p*q tensor
K=length(C);
p=size(C{1},1);
q=size(C{1},2);
T=zeros(K,p,q);
for k=1:K
    T(k,:,:)=C{k};
end
end

function T1=matricize1(T)
% Matricize a K*p*p tensor along the 1st mode to a K*p^2 matrix
T1=reshape(T,size(T,1),size(T,2)*size(T,3));
end

function out=product1(mat,tensor)
% the mode-1 product of k*r matrix mat and r*p*q tensor 
% output is a k*p*q tensor (if k=1, then a p*q matrix)
[k,r]=size(mat);
[r1,p,q]=size(tensor);
if r~=r1
    error('Dimension does not match to perform mode-1 product');
end
out=reshape(mat*matricize1(tensor),k,p,q);
if k==1
    out=permute(out,[2,3,1]);
end
end


