function [J_out,grad] = J_Z(num, k, P, mu, P_hat, Z)
Z_new = Z;
ZZ_new = Z_new - repmat(max(Z_new,[],1),k,1);
rho_new= exp(ZZ_new)*diag(1./sum(exp(ZZ_new),1));
J_out =  J(num, k, rho_new, P, P_hat, mu);
Y = log(P_hat);
grad = dJ_Z_new(num, k, P, mu, Y, Z);
end