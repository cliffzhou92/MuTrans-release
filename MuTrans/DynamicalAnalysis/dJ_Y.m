% The derivative of the objective function J with respect to Y.
% p_hat = diag(1./sum(exp(Y),2))*exp(Y);
% rho = exp(Z)*diag(1./sum(exp(Z),1));
%---------------------------------------------------

function dJY = dJ_Y(num, k, P, mu, Y, Z)
global lamda;

YY = Y - repmat(max(Y,[],2),1,k);
ZZ = Z - repmat(max(Z,[],1),k,1);

p_hat = diag(1./sum(exp(YY),2))*exp(YY);
rho = exp(ZZ)*diag(1./sum(exp(ZZ),1));

%p_hat_2 = diag(1./(sum(exp(Y),2).^2))*exp(Y);

mu_hat = rho*mu;

A = rho'*p_hat*diag(1./mu_hat)*rho - P*diag(1./mu);

dJY1 = rho*diag(mu)*A*diag(mu)*rho'.*p_hat *diag(1./mu_hat);

dJY2 = diag( (rho*diag(mu)*A.*( p_hat*diag(1./mu_hat)*rho ))*mu)*p_hat;

Ent = -lamda * ( p_hat.*( log(p_hat)+1 ) - p_hat.* repmat( sum( p_hat.*( log(p_hat)+1 ),2 ),1,k) );
%Ent =0;

dJY = 2.*(dJY1 - dJY2) + Ent;



