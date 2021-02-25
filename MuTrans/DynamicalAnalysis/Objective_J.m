% Given the transition matrix P and invariant distribution mu, we compute
% the error for the partition "classes" with k partitions.

function E_star = Objective_J( P, mu, P_hat, classes )

num = length(mu);
k = size(P_hat,1);

if k == 1
    E_star = mu'*(P.^2)*(1./mu) - 1.0;

else
    temp_class = repmat(classes',1,k); % size : n * k
    temp_mat = ones(num,1)*(1:k); % size : n * k
    class_mat = temp_class==temp_mat; % size : n * k, Indicator matrix for classes

    [rInd, cInd, vals] = find(class_mat);
    clear class_mat;
    class_mat = sparse(rInd, cInd, vals); % Convert class_mat into a sparse matrix.

    mu_hat = (mu'*class_mat)'; % size : k * 1

    E_star = mu'*(P.^2)*(1./mu) - mu_hat'*(P_hat.^2)*(1./mu_hat);
end







