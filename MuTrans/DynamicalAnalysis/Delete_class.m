% Delete Centre is computated by :

function [classes_d  k_d] = Delete_class( classes, p_hat )

k = size(p_hat,1);

if k==1
    k_d = k;
    classes_d = classes;

else

    k_d = k-1;

    [delete_value, delete_id] = min(diag(p_hat));

    for i = 1:k
        p_hat(i,i)=0;
    end

    [adjacent_value, adjacent_id] = max(p_hat(delete_id,:));

    classes(find(classes==delete_id)) = adjacent_id;


    classes_d = classes;
    if delete_id == 1
        classes_d = classes-1;
    elseif delete_id ==k
        classes_d = classes;
    else

        for i = 1:(delete_id-1)
            classes_d( find(classes==i) ) = i;
        end

        for i = (delete_id+1) : k
            classes_d( find(classes==i) ) = i-1;
        end
    end

end
