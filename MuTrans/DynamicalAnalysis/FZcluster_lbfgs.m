function [rho, p_hat] = FZcluster_lbfgs(class,P,mu,P_hat,par)


k = max(class);  % The number of partitions.
num = size(P,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization of rho by indicator function


Z = -5* ones(k,num);

for i = 1:num
    Z(class(i),i) = 0;
end

options = struct('GradObj','on','Display','iter','LargeScale','on','HessUpdate','lbfgs','InitialHessType','identity','GoalsExactAchieve',0);
J_obj_fun = @(Z)J_Z(num, k, P, mu, P_hat, Z);
Z_new = fminlbfgs(J_obj_fun,Z,options);
ZZ_new = Z_new - repmat(max(Z_new,[],1),k,1);
rho = exp(ZZ_new)*diag(1./sum(exp(ZZ_new),1));
p_hat = P_hat;

end