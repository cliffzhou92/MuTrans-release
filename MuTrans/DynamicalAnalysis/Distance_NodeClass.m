function [E_bar classes_bar k_bar] = Distance_NodeClass( P, mu, P_hat, classes )

k = size(P_hat,1);
num = length(mu);


if k == 1
    k_bar = k;
    classes_bar = ones(1,num);
    E_bar = 1.0;

else
    temp_class = repmat(classes',1,k); % size : n * k
    temp_mat = ones(num,1)*(1:k); % size : n * k
    class_mat = temp_class==temp_mat; % size : n * k, Indicator matrix for classes

    [rInd, cInd, vals] = find(class_mat);
    clear class_mat;
    class_mat = sparse(rInd, cInd, vals); % Convert class_mat into a sparse matrix.

    mu_hat = (mu'*class_mat)'; % size : k * 1

    clear temp_mat;
    temp_mat = P_hat.*(1./mu_hat*ones(1,k))'; % size : k * k

    E_bar = mu.*((P.^2)*(1./mu))*ones(1,k) + mu*((P_hat.^2)*(1./mu_hat))' - 2*(mu*ones(1,k)).*(P*class_mat*temp_mat');

    [value_bar, classes_bar] = min(E_bar');

    %-----------------------------------------------------------

    mm = 0;
    for i = 1:k
        if length(find( classes_bar ==i ))==0
            mm = mm+1;

        else
            classes_bar( find( classes_bar ==i ) ) = i - mm;
        end
    end
    k_bar = k-mm;
end