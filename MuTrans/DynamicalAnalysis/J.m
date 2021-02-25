%The objective function J_E is:

function J0 = J(num, k, rho, P, p_hat, mu)

mu_hat = rho*mu;

Z1 = P./repmat(mu',num,1);

Z2 = rho'*p_hat*diag(1./mu_hat)*rho;

Z = Z1-Z2;

J0 = mu'*(Z.^2)*mu;
