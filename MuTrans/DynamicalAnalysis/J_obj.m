%The objective function J_E is:

function J0 = J_obj(rho_reshape,class_order,P_perm,mu_perm,P_hat)

N_cell = length(class_order);
N_cluster = max(class_order);
rho = reshape(rho_reshape,N_cluster,N_cell);

mu_hat = rho*mu_perm;

Z1 = P_perm./repmat(mu_perm',N_cell,1);

Z2 = rho'*P_hat*diag(1./mu_hat)*rho;

Z = Z1-Z2;

J0 = mu_perm'*(Z.^2)*mu_perm;
