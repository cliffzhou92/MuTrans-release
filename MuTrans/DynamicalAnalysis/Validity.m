% Validity Index1
function V = Validity( P, mu, P_hat, classes )

k = size(P_hat,1);

%P_hat = full(P_hat);

%S = E_star/num/sum(eig(P_hat_star).^2);

%S =  E_star/num/(min(abs(eig(P_hat_star))).^2);

%S =  E_star/num/min(abs(eig(P_hat_star)));

%S =  E_star/num/sum( sum(P_hat_star.^2,2) );

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%(1) max~
% max_p_ii = max(diag(p_hat));
% %min_p_ii = min(diag(P_hat));
% for i = 1:k
%     p_hat(i,i)=0;
% end
% max_p_ij = max(p_hat(:));
%
% if k==1
%     S = E_star;
% else
%    S = E_star/(max_p_ii/max_p_ij);
%   %S = E_star/(min_p_ii/max_p_ij);
% end


%(2) average~


if k==1
    V = mu'*(P.^2)*(1./mu) - 1;

else
    E_star = Objective_J( P, mu, P_hat, classes );

    p_ii = sum( diag( P_hat ) )/k;
    p_ij = ( sum( P_hat(:) )- sum( diag( P_hat ) ) )/(k*(k-1));

    V = E_star * p_ij / p_ii;
end





% The following small part is wrong...
% sum_p_ii = sum(diag(P_hat_star));
% for i = 1:k
%     P_hat_star(i,i)=0;
% end
% sum_p_ij = sum(sum(P_hat_star,2));
%
% if k==1
%     S = E_star;
% else
%     S = E_star/(sum_p_ii/sum_p_ij);
% end


%S = E_star/(sum(eig(P_hat_star).^2)/k);

