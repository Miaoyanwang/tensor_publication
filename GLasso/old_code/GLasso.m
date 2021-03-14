function Omega = GLasso(S,n)
% input:
% S         sample covariance
% n         sample size
% output: 
% Omega     sparse precision matrix
% 
% tuning was selected by BIC


tmp = max(max(abs(S)));
tune=(tmp / 10):(tmp / 10):tmp;
BIC=zeros(size(tune));
for i=1:length(tune)
    rho=tune(i);
    Omega=graphicalLasso(S,rho);
    num=sum(sum(triu(Omega)~=0));
    BIC(i)=2*n*(trace(S*Omega)-log(det(Omega)))+log(n)*num;
end

figure(1);clf;
plot(tune,BIC,'o-');
xlabel('Rho')
ylabel('BIC')
[~,ind]=min(BIC);
Omega=graphicalLasso(S,tune(ind));

end
