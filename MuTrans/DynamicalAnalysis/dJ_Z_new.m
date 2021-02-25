% The derivative of the objective function J with respect to Z.
% p_hat = diag(1./sum(exp(Y),2))*exp(Y);
% rho = exp(Z)*diag(1./sum(exp(Z),1));
%-----------------------------------------------

function dJZ = dJ_Z_new(num, k, P, mu, Y, Z)

YY = Y - repmat(max(Y,[],2),1,k);
ZZ = Z - repmat(max(Z,[],1),k,1);

p_hat = diag(1./sum(exp(YY),2))*exp(YY);
rho = exp(ZZ)*diag(1./sum(exp(ZZ),1));

%rho_2 = exp(Z)*diag(1./(sum(exp(Z),1).^2));

mu_hat = rho*mu;

C = rho'*p_hat*diag(1./mu_hat)*rho -P*diag(1./mu);
A = p_hat*diag(1./mu_hat)*rho;
B = p_hat'*rho;
%----------------------------------

dJZ1 = A*diag(mu)*C'.*rho;

dJZ2 = rho*diag(rho'*A.*C*mu);

dJZ3 = diag(1./mu_hat)*(B*diag(mu)*C).*rho;

dJZ4 = rho*diag(mu'*(rho'*A.*C));

dJZ5 = diag(B*diag(mu)*C.*rho *mu./(mu_hat.^2))*rho;

dJZ6 = rho*diag(mu'* (C*diag(mu)*rho'.*B' *diag(1./(mu_hat.^2))*rho) );

%-----------------------------------

dJZ = 2.*(dJZ1 - dJZ2 + dJZ3 - dJZ4 - dJZ5 + dJZ6)*diag(mu); 

end



