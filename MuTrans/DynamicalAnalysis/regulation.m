% To make sure rho >=0 and sum(rho,1) = (1,...,1)

function re_rho = regulation(num,k,rho)

for j = 1:num

    delta = (1.0-sum(rho(:,j)))/k;
    rho(:,j) = rho(:,j) + repmat(delta,k,1);
    i=0;
    clear idd;
    while 1
        [min_value,minid] = min( rho(:,j) );
        [max_value,maxid] = max( rho(:,j) );
        if  min_value < 0.0
            i = i+1;
            idd(i) = minid;
            rho(minid,j)=0.0;

            num_zero = length(idd);
            num_nonzero = k - num_zero;

            delta = (1.0-sum(rho(:,j)))/num_nonzero;

            rInd = idd;

            cInd = ones(num_zero,1);

            vals = delta.*ones(num_zero,1);

            rho(:,j) = rho(:,j) + repmat(delta,k,1)-sparse(rInd,cInd,vals,k,1);

        else
            break;
        end
    end



end

re_rho =rho;