function out = landscape_gaussian (score_2d, rho_class,mu_hat,label,colors, par)
%landscape function
[N_cell, k] = size(rho_class);
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

figure;
gscatter(score_aver(:,1),score_aver(:,2),label,colors,[],25);
set(gca,'xtick',[],'ytick',[]);

GMModel = gmdistribution(score_center,cov_cluster,mu_hat);
haxis = gca;
xlim = haxis.XLim * scaleaxis;
ylim = haxis.YLim * scaleaxis;
d = (max([xlim ylim])-min([xlim ylim]))/N_mesh;

land_cell = zeros(N_cell,1);
for cell_id = 1:N_cell
land_cell(cell_id) = -log(pdf(GMModel,[score_aver(cell_id,1),score_aver(cell_id,2)]));
end
max_land_cell = max(land_cell);
max_value_lim = scalevalue * max_land_cell; 

x = xlim(1):d:xlim(2);
y = ylim(1):d:ylim(2);
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
    
    level = unique(par.plot_label);
    colors_cell = rho_class*colors;
%     for cluster_id = 1:length(level)
%         id = par.plot_label == level(cluster_id) ;
%         plot(score_aver(id,1),score_aver(id,2),'.','color',colors(cluster_id,:),'markersize',mksize); 
%         hold on 
%         set(gca,'layer','top')
%     end

    for cell_id = 1:size(rho_class,1)
        plot(score_aver(cell_id,1),score_aver(cell_id,2),'.','color',colors_cell(cell_id,:),'markersize',mksize); 
        hold on 
        set(gca,'layer','top')
    end
    
else
    ld = surf(X1Grid,X2Grid,Land');
    hold on
    alpha(ld,par.alpha)
    grid off
    shading interp    
    
    level = unique(par.plot_label);
%     for cluster_id = 1:length(level)
%         id = par.plot_label == level(cluster_id) ;
%         plot3(score_aver(id,1),score_aver(id,2),land_cell(id),'.','color',colors(cluster_id,:),'markersize',mksize); 
%         hold on 
%         set(gca,'layer','top')
%     end

        colors_cell = rho_class*colors;
        for cell_id = 1:size(rho_class,1)
        plot3(score_aver(cell_id,1),score_aver(cell_id,2),land_cell(cell_id),'.','color',colors_cell(cell_id,:),'markersize',mksize); 
        hold on 
        set(gca,'layer','top')
        end
    
    
end

if par.legend
legend_item = get(gca,'Children');
lg = legend(flip(legend_item(1:end-1)),par.clustername,'fontsize',par.fontsize); 
%legend boxoff
end
axis off
tickCell = {'XTickLabel',{},'YTickLabel',{},'ZTickLabel',{}};
set(gca,tickCell{:});
out.land_cell = land_cell;
out.land = Land;
out.model = GMModel;
out.colors_cell = colors_cell;
end