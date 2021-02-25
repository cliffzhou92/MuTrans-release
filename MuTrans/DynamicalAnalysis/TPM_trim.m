function [P_trim,theta_matrix] = TPM_trim (P_hat,class_order, rho_class, par)
k = size(P_hat,1);
theta_matrix = zeros(k,k);
P_trim = P_hat;
    for i = 1:(k-1)
        for j = (i+1):k
        theta_matrix (i,j) = 0.5*(logistic_transition_score(i, j ,class_order, rho_class,par)+logistic_transition_score(j, i ,class_order, rho_class, par));
        theta_matrix (j,i) = theta_matrix (i,j);
        end
    end
trim_matrix = theta_matrix < par.thresh_sharp;
P_trim(~trim_matrix) = 0;
figure;
hist(theta_matrix);
end