function Output = DynamicalAnalysis (data, par)
% MuTrans Dynamical Analysis of Single cell Data
% Input : 
%   data : Single-cell expression matrix
%   par.K_cluster : Numeric paramter. Number of cell clusters.
%   par.options : Whether to allow an automatic selection of
%   K_cluster (might be slow, not recommended). Default is 'fixed'.
%   par.perplex : Numeric paramter in constructing cell-cell rwTPM. The perplexity determines the local variance gaussian-kernel similarity between cells. Default is N_cell/2.5.
%   par.choice_distance : Numeric paramter in constructing cell-cell rwTPM.
%   The metric in gaussian-kernel similarity between cells. Default is
%   Euclidean. Can also be Cosine or Correlation.
%   par.weight_scale: Logic parameter in constructing cell-cell rwTPM. Wheter
%   to scale the constructed rwTPM. Default is true.
%   par.alpha_scale: Logic parameter in constructing cell-cell rwTPM. The power appear in the denominator to scale rwTPM as in diffusion map. Default is 0.
%   par.initial: The initial cell assignment in coarse-graining of cell-cell rwTPM. Default is tSNE+
%   kemans. Can also be 'random' or 'pca'.
%   par.trials: Number of trials in coarse-graining of cell-cell rwTPM.
%   Default is 20.
%   par.softclustering: Whether performs refinement optimization to obtain cell-cluster soft clustering membership.
%   Default is true.
%   par.opt_tool: Optimization tool in soft-clustering. Default is using fminunc function in matlab. Can also
%   be 'lbfgs'.
%   par.tol: Optimality Tolerance parameter in optimization. Default is 1e-6.
%   par.step: Maximum iteration steps allowed in optimization. Default is
%   1000.
%   par.print_results: Whether to output the clustering results in tSNE.
%   Default is false.
%   par.true_labs: The true experimental or original labels of the cells.

% Output :
%   class_order : The cluster label of each individual cell through
%   coarse-graining, with cells re-ordered by their label number.
%   P_hat : The cluster-cluster rwTPM through coarse-graining.
%   rho_class : The cell-cluster membership matrix through dynamics refinement,with cells re-ordered by their label number.
%   mu_perm : The re-ordered stationary distribution of rwTPM at cell-cell
%   level.
%   P_perm : The re-ordered rwTPM at cell-cell level.
%   P_rho: The refined rwTPM by considering uncertainty in cell-cluster
%   membership.
%   H: The Shannon entropy of membership function.
%   data_perm: The re-ordered input gene expression matrix.
%   labs_perm: The re-ordered input experimental or original labels.

data = double(data);
[N_cell, ~] = size(data);

if ~isfield(par,'perplex')
    par.perplex = N_cell/2.5;
end

if ~isfield(par,'choice_distance')
    par.choice_distance = 'euclid';
end

if ~isfield(par,'initial')
    par.initial = 'tsne';
end

if ~isfield(par,'trials')
    par.trials = 20;
end

if ~isfield(par,'opt_tool')
    par.opt_tool = 'matlab';
end

if ~isfield(par,'softclustering')
par.softclustering = true;
end

if ~isfield(par,'weight_scale')
par.weight_scale = true;
end

if ~isfield(par,'print_results')
    par.print_results = false;
end

if ~isfield(par,'alpha_scale')
    par.alpha_scale = 0;
end

if ~isfield(par,'options')
    par.options = 'fixed';
end

if ~isfield(par,'tol')
    par.tol = 1e-6;
end

if ~isfield(par,'step')
    par.step = 1000;
end

if ~isfield(par,'true_labs')
    par.true_labs = ones(N_cell,1);
end

if ~isfield(par,'reduce_large_scale')
    par.reduce_large_scale = false;
end

if ~isfield(par,'reduce_num_meta_cell')
    par.reduce_num_meta_cell = 1000;
end

par.data = data;


%% Construction of the cell-to-cell random walk, E is the weight matrix

perplex = par.perplex; % in general one may choose N_cell/2.5;200

Dist = squareform(pdist (data,par.choice_distance));

if par.weight_scale
    [E, ~] = d2p(Dist, perplex, 1e-4);
else
    [E, ~] = d2p_new(Dist, perplex, 1e-4);
    sum_col = sum(E,1);
    sum_row = sum(E,2);
    alpha_scale= par.alpha_scale;
    scale_col = diag(sum_col.^(-alpha_scale));
    scale_row = diag(sum_row.^(-alpha_scale));
    E = scale_row*E*scale_col;
end

E = 0.5*(E+E');



%% Coarse Graining
    
% parameter for simulated anealling if the number of cluster is not
% specified
    if ~isfield(par,'K_cluster')
    par.options = 'optimal';
    par.K_max = floor(N_cell/100);
    par.K_min = 3;
    par.T_max = 3;
    par.T_min = 1e-3;
    par.Rloop = 50;  % The iterative times at each temperature
    par.alpha = 0.9;  % The cooling rate
    end


out_cg = CoarseGraining(E,par);

class_out = out_cg.class_out;
Output.P_hat = out_cg.P_hat;
network_out = out_cg.network_out;

if par.reduce_large_scale
    Output.data = out_cg.data;
    Output.reduce_class =out_cg.reduce_class;
end

P = network_out.P;
mu = network_out.mu;
k = max(class_out);
mu_hat = zeros(k,1);
for i =1:1:k
    mu_hat(i) = sum(mu(class_out == i));
end
Output.mu_hat = mu_hat;
Output.k = k ;
% Reorder of the cells accoring to cluster label
N_cell = size(P,1);
P_appr = zeros(N_cell,N_cell);
for i = 1: N_cell
    for j = 1: N_cell
    cluster_i = class_out(i);
    cluster_j = class_out(j);
    P_appr(i,j) = Output.P_hat(cluster_i,cluster_j)*mu(j)/mu_hat(cluster_j);
    end
end

[Output.class_order,perm_class] = sort(class_out);
Output.perm_class = perm_class;
Output.P_appr_perm = P_appr(perm_class,perm_class);
Output.P_perm = P(perm_class,perm_class);
Output.mu_perm = mu(perm_class);
Output.data_perm = data(perm_class,:);
if isfield(par,'true_labs')
    Output.labs_perm = par.true_labs(perm_class);
end
%% Soft (Fuzzy) Clustering of the data
if par.softclustering 
    tic;
    
    if ~isfield(par,'tol') || isempty(par.tol)
    par.tol = 1e-6;
    end
    if ~isfield(par,'step') || isempty(par.step)
    par.step = 1000;
    end
    
    if strcmp(par.opt_tool , 'matlab')
    [rho, ~] = FZcluster_matlab(Output.class_order,Output.P_perm,Output.mu_perm,Output.P_hat,par);
    else 
    [rho, ~] = FZcluster_lbfgs(Output.class_order,Output.P_perm,Output.mu_perm,Output.P_hat,par);    
    end
    Output.rho_class = rho';
    Output.H = -sum(Output.rho_class' .* log2((Output.rho_class+eps)'));
    Output.P_rho = rho'*Output.P_hat*diag(1./Output.mu_hat)*rho*diag(Output.mu_perm);
    toc;
end

if par.print_results
   
   if  ~isfield(par,'score')
   Dist = squareform(pdist (data,par.choice_distance));
   par.score = tsne_d(Dist);
   end
   
   score = par.score;
   figure;
   gscatter(score(:,1),score(:,2),par.true_labs,[],[],22); 
   figure;
   gscatter(score(:,1),score(:,2),class_out,[],[],22); 
end

end