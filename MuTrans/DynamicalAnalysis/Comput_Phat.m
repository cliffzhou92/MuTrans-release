function P_hat = Comput_Phat( P, mu, classes, k )

num = length(mu);
if k == 1
    P_hat = 1.0;

else
    % Comput P_hat :
    temp_class = repmat(classes',1,k); % size : n * k
    temp_mat = ones(num,1)*(1:k); % size : n * k
    class_mat = temp_class==temp_mat; % size : n * k, Indicator matrix for classes

    [rInd, cInd, vals] = find(class_mat);
    clear class_mat;
    class_mat = sparse(rInd, cInd, vals); % Convert class_mat into a sparse matrix.

    mu_hat = (mu'*class_mat)'; % size : k * 1
    mu_hat_ext = class_mat*mu_hat; % size : n * 1
    mu_i = mu./mu_hat_ext; % size : n * 1
    mu_i_ext = repmat(mu_i,1,k); % size : n * k

    P_hat = (mu_i_ext.*class_mat)'*P*class_mat; % size : k * k
end