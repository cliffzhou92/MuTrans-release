function out = bifur_genes_within (cluster_ids, data_perm,class_order,rho_class,par)
genes = par.genes;
cell_id = par.cells;
true_labs = par.true_labs;

% whole data, select genes, 
id_whole = ismember(class_order,cluster_ids);
id_bifur = ismember(class_order,cluster_ids(1));
data_whole = data_perm(id_whole,:);
data_bifur = data_perm(id_bifur,:);
class_order_ini_whole = class_order(id_whole);
id_type_whole = ismember(1:max(class_order),cluster_ids);
id_ntype_whole = ~id_type_whole;
rho_other_whole = rho_class(id_whole,id_ntype_whole);
rho_other_bifur = rho_class(id_bifur,id_ntype_whole);  % must be id_data_keep_whole


% only considers the cells in three states (useful if g1 is the bifurcation state)
id_data_keep_whole = max(rho_other_whole,[],2)<par.otherkeep;
id_data_keep_bifur = max(rho_other_bifur,[],2)<par.otherkeep;
data_whole = data_whole(id_data_keep_whole,:);
data_bifur = data_bifur(id_data_keep_bifur,:);
rho_ini_bifur = rho_class(id_bifur,cluster_ids);
class_order_whole = class_order_ini_whole(id_data_keep_whole);
rho_bifur = rho_ini_bifur(id_data_keep_bifur,:);
labs_ini = true_labs(id_bifur);
labs_bifur = labs_ini(id_data_keep_bifur);


% initial sort of DE genes  
[~,N_genes_whole] = size(data_whole);
p_anova = nan(N_genes_whole,1);

for k = 1:N_genes_whole
    gene_temp = data_whole(:,k);
    p_anova(k) = anova1(gene_temp,class_order_whole,'off');
end

gene_keep_id = p_anova < par.de_pvalues;
genes_keep = genes(gene_keep_id);
data_whole_left = data_whole(class_order_whole== cluster_ids(2),gene_keep_id);
data_whole_right = data_whole(class_order_whole== cluster_ids(3),gene_keep_id);
% sort the cells in bifurcation state
score_type = [0,-1,1];
[~,type_shuffle] = sort(cluster_ids);
score_type  = score_type(type_shuffle);
cell_score = rho_bifur*score_type';
[cell_score_sort,cell_id_sort] = sort(cell_score);
data_bifur_sort = data_bifur(cell_id_sort,gene_keep_id);
labs_bifur_sort = labs_bifur(cell_id_sort);
id_within_stable= abs(cell_score_sort) < 0.02;
id_within_left = cell_score_sort < -0.02;
id_within_right = cell_score_sort > 0.02;
data_within_left = data_bifur_sort(id_within_left,:);
data_within_right = data_bifur_sort(id_within_right,:);
data_within_stable = data_bifur_sort(id_within_stable,:);


[p_in_left,t_in_left] = mattest(data_within_left',data_within_stable');
[p_in_right,t_in_right] = mattest(data_within_right',data_within_stable');
[p_out_left,t_out_left] = mattest(data_whole_left',data_bifur_sort');
[p_out_right,t_out_right] = mattest(data_whole_right',data_bifur_sort');

significant_in_left = p_in_left < par.de_pvalues_inbifur;
significant_in_right = p_in_right < par.de_pvalues_inbifur;
significant_out_left = p_out_left < par.de_pvalues_inbifur;
significant_out_right = p_out_right < par.de_pvalues_inbifur;
consist_left = significant_in_left-significant_out_left;
consist_right = significant_in_right-significant_out_right;
delete_left = consist_left ==1;
delete_right = consist_right ==1;
delete_all = logical(delete_left+delete_right);
gene_keep_id = ~delete_all;
genes_keep = genes_keep(gene_keep_id);
data_bifur_sort = data_bifur_sort(:,gene_keep_id); 
significant_in_left = significant_in_left (gene_keep_id);
t_in_left =  t_in_left(gene_keep_id);
significant_in_right = significant_in_right(gene_keep_id);
t_in_right = t_in_right(gene_keep_id);

[~,N_genes_bifur] = size(data_bifur_sort);


%marker_gene_type \in [1,2,3,4];
%1: 100 011 2:010 101 3:001 110 4: 000 111

marker_type = nan(N_genes_bifur,1);
t_value = nan(N_genes_bifur,1);
for k = 1:N_genes_bifur
    sig_in_diff = significant_in_left(k)-significant_in_right(k);
    switch sig_in_diff
        case 1
            marker_type (k) = 1;
            t_value(k) = t_in_left(k);
        case 0
            if significant_in_left(k)*significant_in_right(k)==1
                marker_type (k) = 2;
                t_value(k) = 0.5*(t_in_left(k)+t_in_right(k));
            else
                marker_type (k) = 4;
                t_value(k) = 0.5*(t_in_left(k)+t_in_right(k));
            end
        case -1
            marker_type (k) = 3;
            t_value(k) = t_in_right(k);
    end
end

[~,gene_sort_id] = sortrows([marker_type, t_value]);
data_bifur_sort = normalize_data(data_bifur_sort(:,gene_sort_id));
data_smooth = smooth_data(data_bifur_sort,par.smooth_span,par.smooth_method);

genes_sort = genes_keep(gene_sort_id);
marker_type_sort = marker_type(gene_sort_id);
figure('rend','painters','pos',[10 10 900 N_genes_bifur*par.unit_length])
colormap redbluecmap;
clims = [-2 2];
imagesc(data_smooth',clims);
hold on 
%plot(get(gca, 'xlim'),[id_sig_up+0.5 id_sig_up+0.5]);
set(gca,'ytick',[]);
yticks(1:length(genes_sort));
set(gca,'xtick',[]);
yticklabels(genes_sort)
set(gca,'TickLength',[0 0])

out.labs = labs_bifur_sort;
out.marker_type = marker_type_sort;
end