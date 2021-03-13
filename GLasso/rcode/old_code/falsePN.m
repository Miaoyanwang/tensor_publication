function [FP,FN]=falsePN(Omega,Omega_est)
% the function calculates false positive and false negative rates based on
% true network Omega and estimated network Omega_est (ignoring diagonal)

falsepositive= (Omega_est~=0 & Omega==0);
num_FP=sum(sum(falsepositive-diag(diag(falsepositive))))/2;
negative=(Omega==0);
num_N=sum(sum(negative-diag(diag(negative))))/2;

FP=num_FP/num_N;


falsenegative = (Omega_est==0 & Omega~=0);
num_FN=sum(sum(falsenegative-diag(diag(falsenegative))))/2;
positive=(Omega~=0);
num_P=sum(sum(positive-diag(diag(positive))))/2;

FN=num_FN/num_P;