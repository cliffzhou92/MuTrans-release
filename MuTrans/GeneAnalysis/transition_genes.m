function out = transition_genes (cluster_id1, cluster_id2, data_perm,class_order,rho_class,par)
% cluster_id1, cluster_id2 start and target

genes = par.genes;
cell_id = par.cells;
true_labs = par.true_labs;

g1 = cluster_id1;
g2 = cluster_id2;


%% pre-processing the data
id_data = ismember(class_order,[g1 g2]);
data_ini = data_perm(id_data,:);
class_order_ini = class_order(id_data);
id_type = ismember(1:max(class_order),[g1 g2]);
id_ntype = ~id_type;
rho_other = rho_class(id_data,id_ntype);
true_labs_ini = true_labs(id_data);
cell_id_ini = cell_id(id_data);

% only considers the cells transit between two states (useful if g1 is the bifurcation state)
id_data_keep = max(rho_other,[],2)<par.otherkeep;
data = data_ini(id_data_keep,:);
rho_ini = rho_class(id_data,:);
rho = rho_ini(id_data_keep,:);
class_order_keep = class_order_ini(id_data_keep);
true_labs_keep = true_labs_ini(id_data_keep);
cell_id_keep = cell_id_ini(id_data_keep);
% delete the genes with low variance
std_genes = std(data);
ndelete_gene = std_genes > 0.01;
genes = genes(ndelete_gene);
data = data(:,ndelete_gene);

% two groups of data
bifurid_g1 = class_order_keep == g1;
bifurid_g2 = class_order_keep == g2;
data_g1 = data(bifurid_g1,:);
data_g2 = data(bifurid_g2,:);
rho_g1 = rho(bifurid_g1,:);
rho_g2 = rho(bifurid_g2,:);


score_g1 = rho_g1(:,g1)./(rho_g1(:,g1)+rho_g1(:,g2));
score_g2 = (rho_g2(:,g1))./(rho_g2(:,g1)+rho_g2(:,g2));
[score1,sort1] = sort(score_g1,'descend');
[score2,sort2] = sort(score_g2,'descend');

lab_g1 = true_labs_keep(bifurid_g1);
lab_g2 = true_labs_keep(bifurid_g2);
cell_id_g1 = cell_id_keep(bifurid_g1);
cell_id_g2 = cell_id_keep(bifurid_g2);
lab_sort = [lab_g1(sort1);lab_g2(sort2)];
cell_id_sort = [cell_id_g1(sort1);cell_id_g2(sort2)];


%sort the cells
[PValues, TValues] = mattest(data_g1',data_g2');
select_id = PValues < par.de_pvalues;%0.05
genes = genes(select_id);
data_g1 = data_g1(:,select_id);
data_g2 = data_g2(:,select_id);
data_g1_sort = data_g1(sort1,:);
data_g2_sort = data_g2(sort2,:);
t_keep = TValues(select_id);
data_all = [data_g1_sort;data_g2_sort];
score_logistic = [score1;score2];
id_sep = length(sort1);


%% logistic transition
% data normalization

%sort the genes
[N_cell, N_genes] = size(data_all);

L = length(score_logistic);
h = 1/(L-1);
x_fit = 0:h:1;
[score_pre, p_logistic, ~, ~] = fit_logistic(x_fit',score_logistic);


% the transition score of genes
id_trans = abs(score_pre-0.5)< par.trans_thresh;
data_trans = data_all(id_trans,:);
data_trans = normalize_data(data_trans);
data_trans_smooth = smooth_data(data_trans,par.smooth_span_trans,par.smooth_method);
gene_score_trans = zeros(N_genes,1);

for k = 1: N_genes
gene_temp = data_trans_smooth(:,k);
gene_max = max(gene_temp);
gene_min = min(gene_temp);
gene_temp = (gene_temp-gene_min)/(gene_max-gene_min);
gene_score_trans(k) = corr(gene_temp,score_pre(id_trans),'Type',par.corr_choice);
    if gene_score_trans(k)* t_keep(k) < 0
    gene_score_trans(k) = 0;
    end
end
%% trans

gene_keep_id = abs(gene_score_trans)> par.score_thresh;
gene_score_trans = gene_score_trans(gene_keep_id);

data_trans_keep = data_trans_smooth(:,gene_keep_id);
genes_trans = genes(gene_keep_id);

[gene_score_sort,gene_sort]= sort(gene_score_trans);
data_all_sort_trans = data_trans_keep(:,gene_sort);
gene_name_sort = genes_trans(gene_sort);


L_genes_trans = length(gene_name_sort);
l = par.unit_length;
L_tot = 900;
trans_start =  sum((score_pre-0.5)> par.trans_thresh);
l_cell = L_tot/N_cell;


figure('rend','painters','pos',[10 10 L_tot  L_genes_trans*l])
score_pre_plot = score_pre;
if par.flip
    score_pre_plot = flip(score_pre_plot);
    trans_start = N_cell - trans_start;
end
plot(x_fit,score_pre_plot,'linewidth',2.0)
hold on
xlim([0 1])
ylim([0 max(score_pre_plot)])
plot(h*[id_sep id_sep], ylim)
plot(h*[trans_start trans_start], ylim)
plot(h*[trans_start+sum(id_trans) trans_start+sum(id_trans)], ylim)

set(gca,'xtick',[]);
set(findall(gcf,'-property','FontSize'),'FontSize',18)
box off
set(gca,'color','none')

figure('OuterPosition',[10+l_cell*trans_start 10 sum(id_trans)*l_cell L_genes_trans*l])
colormap redbluecmap;
clims = [-2 2];
if par.flip
    data_all_sort_trans = flip(data_all_sort_trans);
end
imagesc(data_all_sort_trans',clims);
set(gca,'ytick',[]);
set(gca,'xtick',[]);
if(par.genes_label)
yticks(1:length(gene_name_sort));
yticklabels(gene_name_sort)
set(gca,'TickLength',[0 0])
end
%{
data_scale_org = data_scale(:,gene_keep_id);
data_all_sort_org = data_scale_org(:,gene_sort);
figure('rend','painters','pos',[10 10 400 l*L])
colormap redbluecmap;
clims = [-3 3];
imagesc(data_all_sort_org',clims);
set(gca,'ytick',[]);
yticks(1:length(gene_name_sort));
set(gca,'xtick',[]);
yticklabels(gene_name_sort)
%}
trans_up_id = gene_score_sort < 0 ;
trans_down_id = gene_score_sort > 0 ;
out.genes_trans_up = gene_name_sort(trans_up_id);
out.genes_trans_down = gene_name_sort(trans_down_id);
%% differential expression
%
id_g1_stable = (score_pre-0.5)> par.trans_thresh;
id_g2_stable = (0.5-score_pre)> par.trans_thresh;

n_keep = ~gene_keep_id;
genes_n_trans = genes(n_keep);
data_g1_sort = data_g1_sort(:,n_keep);
data_g2_sort = data_g2_sort(:,n_keep);
data_all = [data_g1_sort;data_g2_sort];

data_all = normalize_data(data_all);
[PValues, TScores] = mattest(data_g1_sort',data_g2_sort');
gene_keep_id = PValues < par.de_pvalues;
t_keep = TScores(gene_keep_id);
p_keep = PValues(gene_keep_id);
data_keep = data_all(:,gene_keep_id);
genes_keep = genes_n_trans(gene_keep_id);


id_up = t_keep>0;
data_up = data_keep(:,id_up);
p_up = p_keep(id_up);
genes_up = genes_keep(id_up);

data_up_stable = data_up(id_g1_stable,:);
data_up_trans = data_up(id_trans,:);
[p_marker_up,~] = mattest(data_up_stable',data_up_trans');
[p_marker_up_sort,id_up_sort] = sort(p_marker_up,'ascend');
id_sig_up = p_marker_up_sort > par.de_significant;
data_up_sort = data_up(:,id_up_sort);
genes_up_sort = genes_up(id_up_sort);
genes_up_diff = genes_up_sort(~id_sig_up);
genes_up_undif = genes_up_sort(id_sig_up);
data_smooth = smooth_data(data_up_sort,par.smooth_span,par.smooth_method);

L_genes_up = length(genes_up);
figure('rend','painters','pos',[10 10 900 L_genes_up*l])
colormap redbluecmap;
clims = [-2 2];
if par.flip
    data_smooth = flip(data_smooth);
end
imagesc(data_smooth',clims);
hold on 
plot(get(gca, 'xlim'),[sum(~id_sig_up)+0.5 sum(~id_sig_up)+0.5]);
set(gca,'ytick',[]);
set(gca,'xtick',[]);
if(par.genes_label)
yticks(1:length(genes_up_sort));
yticklabels(genes_up_sort)
end
set(gca,'TickLength',[0 0])


id_down = t_keep<0;
data_down = data_keep(:,id_down);
p_down = p_keep(id_down);
genes_down = genes_keep(id_down);


data_down_stable = data_down(id_g2_stable,:);
data_down_trans = data_down(id_trans,:);
[p_marker_down,~] = mattest(data_down_stable',data_down_trans');
[p_marker_down_sort,id_down_sort] = sort(p_marker_down,'descend');
id_sig_down = p_marker_down_sort < par.de_significant;
data_down_sort = data_down(:,id_down_sort);
genes_down_sort = genes_down(id_down_sort);
genes_down_diff= genes_down_sort(id_sig_down);
genes_down_undiff= genes_down_sort(~id_sig_down);

data_smooth = smooth_data(data_down_sort,par.smooth_span,par.smooth_method);

L_genes_down = length(genes_down);
figure('rend','painters','pos',[10 10 900 L_genes_down*l])
colormap redbluecmap;
clims = [-2 2];
if par.flip
    data_smooth = flip(data_smooth);
end
imagesc(data_smooth',clims);
hold on 
plot(get(gca, 'xlim'),[sum(~id_sig_down)+0.5 sum(~id_sig_down)+0.5]);


set(gca,'ytick',[]);
set(gca,'xtick',[]);
if(par.genes_label)
yticks(1:length(genes_down_sort));
yticklabels(genes_down_sort)
set(gca,'TickLength',[0 0])
end
out.genes_ms_g1 = genes_up_diff;
out.genes_ms_g2 = genes_down_diff;
out.genes_ih_g1 = genes_up_undif;
out.genes_ih_g2 = genes_down_undiff;

%}
%% output
out.score_cell = score_logistic;
out.score_gene = gene_score_sort;   
out.cell_true_lab = lab_sort;
out.cell_id_sort = cell_id_sort;
out.p = p_logistic;
end