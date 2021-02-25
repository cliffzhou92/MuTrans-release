function out = EstClusterNum(data,par)
k = par.order;

[N,~] = size(data);
perplex = N/2.5;

Dist = squareform(pdist(data,par.choice_distance));
[E, ~] = d2p(Dist, perplex, 1e-4);
E = (E+E')/2;
D = diag(sum(E));
W = D^(-0.5)*E*D^(-0.5);
lambda = sort(abs(eig(W)),'descend');
lambda_k = lambda.^k;
sum_all = sum(lambda_k);
sum_cul = 0;
prop = zeros(N,1);
ratio = zeros(N-1,1);

    for i = 1:N
        sum_cul = sum_cul + lambda_k(i);
        prop(i) = sum_cul/sum_all;
        if i<N
            ratio(i) = lambda_k(i)/lambda_k(i+1);
        end
    end

out.lambda_k = lambda_k;
out.prop = prop;
out.ratio = ratio;
end