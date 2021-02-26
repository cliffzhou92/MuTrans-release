function out = ConstructLandscape(DA_results, par)
% Construct Cell-fate Landscape
% Input:
%   DA_results: The output structure generated by DynamicalAnalysis
%   function. It must contains a field called embedding_2d, which is the user-provided 2D dimension reduction of data.
%   par.label: The input clustering results for landscape construction. Default is the cluster by MuTrans coarse-graining
%   par.thresh_calc_center: threshold of minimum membership function to
%   include for estimating local means of each cluster.
%   par.thresh_calc_cov: threshold of minimum membership function to
%   include for estimating local covariance of each cluster.
%   par.N_mesh: the number of grids to calculate and plot smooth landscape.
%   par.scaleaxis: the scale of x-y axis to plot landscape.
%   par.scalevalue: the scale of z axis to plot landscape.
%   par.colors: colormap of the clusters.
%   par.color_mixing: whether display the cell colors according to the
%   mixutre of membership. Default is true.
%   par.display_legend: whether display the legend of clusters.
%   par.legend_text: text of the legends.
%   par.plot_label: labels to plot on the landscape. Useful when need to
%   play the experimental labels.

%% default parameters

rho_class = DA_results.rho_class;
score_2d = DA_results.embedding_2d;
mu_hat = DA_results.mu_hat;


if ~isfield(par,'label')
    par.label = DA_results.class_order;
end

[N_cell, k] = size(rho_class);

if ~isfield(par,'colors')
    par.colors = brewermap(k,'set1');
end

if ~isfield(par,'thresh_calc_center')
    par.thresh_calc_center = 0.6;
end

if ~isfield(par,'thresh_calc_cov')
    par.thresh_calc_cov = 0.3;
end

if ~isfield(par,'N_mesh')
    par.N_mesh = 500;
end

if ~isfield(par,'scaleaxis')
    par.scaleaxis = 1.2;
end

if ~isfield(par,'scalevalue')
    par.scalevalue = 1.5;
end

if ~isfield(par,'fontsize')
    par.fontsize = 15;
end

if ~isfield(par,'visual_dim')
    par.visual_dim = 3;
end

if ~isfield(par,'alpha')
    par.alpha = 0.3;
end

if ~isfield(par,'mksize')
    par.mksize = 30;
end

if ~isfield(par,'display_legend')
par.display_legend = false;
end


if ~isfield(par,'color_mixing')
par.color_mixing = true;
end

if ~isfield(par,'plot_label')
par.plot_label = par.label;
end

if ~isfield(par,'plot_landscape')
par.plot_landscape = true;
end

if ~isfield(par,'reduce_large_scale')
    par.reduce_large_scale = false;
end

colors = par.colors;
label = par.label;
level = unique(par.plot_label);



%%
    if par.reduce_large_scale
        classes = DA_results.reduce_class;
        k_meta_cell = par.reduce_num_meta_cell;
        score_2d_reduced = zeros(k_meta_cell,2);
        for i_g = 1:k_meta_cell
            score_2d_reduced(i_g,:) = mean(score_2d(classes==i_g,:));
        end
        score_2d = score_2d_reduced(DA_results.perm_class,:);
    end
    
%%
score_center = zeros(k,2);
cov_cluster = zeros(2,2,k);
thresh_calc_center = par.thresh_calc_center; 
thresh_calc_cov = par.thresh_calc_cov;
N_mesh = par.N_mesh;
mksize = par.mksize;
scaleaxis = par.scaleaxis;
scalevalue = par.scalevalue;
for cluster_id = 1:k
    member_id = label==cluster_id;
    stable_id = rho_class(:,cluster_id)> thresh_calc_center;
    select_id = logical(max(member_id',stable_id))';
    score_center(cluster_id,:) = mean(score_2d(select_id,:));
end

    score_aver = rho_class* score_center; 

for cluster_id = 1:k
    member_id = label ==cluster_id;
    stable_id = rho_class(:,cluster_id)> thresh_calc_cov;
    select_id = logical(max(member_id',stable_id))';
    cov_cluster(:,:,cluster_id) = cov(score_aver(select_id,:));
end


GMModel = gmdistribution(score_center,cov_cluster,mu_hat);
score_aver_x = score_aver(:,1);
score_aver_y = score_aver(:,2);
xlim = scaleaxis*[min(score_aver_x),max(score_aver_x)];
ylim = scaleaxis*[min(score_aver_y),max(score_aver_y)];

d_x = (xlim(2)-xlim(1))/N_mesh;
d_y = (ylim(2)-ylim(1))/N_mesh;
land_cell = zeros(N_cell,1);
for cell_id = 1:N_cell
land_cell(cell_id) = -log(pdf(GMModel,[score_aver_x(cell_id),score_aver_y(cell_id)]));
end
max_land_cell = max(land_cell);
max_value_lim = scalevalue * max_land_cell; 

x = xlim(1):d_x:xlim(2);
y = ylim(1):d_y:ylim(2);
L_x = length(x);
L_y = length(y);
Land = zeros(L_x,L_y);
for i = 1:L_x
    for j = 1:L_y
    Land_value = -log(pdf(GMModel,[x(i),y(j)]));
        if Land_value > max_value_lim
        Land(i,j) = nan;
        else
        Land(i,j) = Land_value;
        end
    end
end



if(par.plot_landscape)
    figure;
    gscatter(score_aver_x,score_aver_y,label,colors,[],mksize);
    set(gca,'xtick',[],'ytick',[]);
    
    figure;
    [X1Grid,X2Grid] = meshgrid(x,y);
    colormap(flipud(brewermap([],'greys')));

    if ~isfield(par,'visual_dim')
    par.visual_dim = 3;
    end


    if par.visual_dim == 2
        pcolor(X1Grid,X2Grid,Land');
        hold on
        alpha(0.5)
        grid off
        shading interp    

        colors_cell = rho_class*colors;

        for cell_id = 1:size(rho_class,1)
            plot(score_aver_x(cell_id),score_aver_y(cell_id),'.','color',colors_cell(cell_id,:),'markersize',mksize); 
            hold on 
            set(gca,'layer','top')
        end

    else
        ld = surf(X1Grid,X2Grid,Land');
        hold on
        alpha(ld,par.alpha)
        grid off
        shading interp    


        if par.color_mixing 

            colors_cell = rho_class*colors;
            for cell_id = 1:size(rho_class,1)
            h = scatter3(score_aver_x(cell_id),score_aver_y(cell_id),land_cell(cell_id),mksize,colors_cell(cell_id,:),'filled'); 
            h.MarkerFaceAlpha = 0.8;
            hold on 
            set(gca,'layer','top')
            end

        else
            for cluster_id = 1:length(level)
             id = par.plot_label == level(cluster_id) ;
             plot3(score_aver_x(id),score_aver_y(id),land_cell(id),'.','color',colors(cluster_id,:),'markersize',mksize); 
             hold on 
             set(gca,'layer','top')
            end
        end
    end

    %% figure legend
    if par.display_legend
                h = zeros(length(par.legend_text), 1);
                for legend_id = 1:length(par.legend_text)
                h(legend_id) = plot(NaN,NaN,'.','color',colors(legend_id,:),'markersize',mksize);
                end
                legend(h,par.legend_text);
    end

    axis off
    tickCell = {'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{}};
    set(gca,tickCell{:});
end
%% output
out.trans_coord = score_aver;
out.land_cell = land_cell;
out.land = Land;
out.model = GMModel;
out.grid_x = x;
out.grid_y = y;
end