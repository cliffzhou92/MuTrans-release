% Split Centre is computated by :

function [ classes_s  k_s] =  Split_class( classes, p_hat, mu )

%global K_max p_ref;
global K_max;

k = size(p_hat,1);

[split_value, split_id] = max(diag(p_hat));

node_id = find(classes==split_id);

%mu_min = min( mu(node_id) );
sN = length( mu(node_id) );

%mu_average = sum( mu(node_id) )/mu_N;
%---------------------------------------------------------------


%if mu_average*p_ref < mu_min | k==K_max | mu_N<4
 
if  k==K_max | sN<=2
    k_s = k;
    classes_s = classes;
else

    k_s = k+1;
    choose_id = fix( (sN-1)*rand( 1, fix(sN/2) ) )+1;
    
    classes_s = classes;
    classes_s( node_id(choose_id) ) = k+1;

end
