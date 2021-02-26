function FigH = PlotGeneTrend (gene,par)
if ~any(isnan(gene))
    L_cells = length(gene);
    L_window = min(30,L_cells/10);
    cell_id = 1:L_cells;
    gene_smooth = smoothdata(gene,'gaussian',L_window);
    f = fit(cell_id',gene_smooth,'smoothingspline','SmoothingParam',0.1);
    Fig = plot(cell_id,f(cell_id),'Linewidth',5.0);
    FigH  = ancestor(Fig, 'figure');
    hold on
    fig_gene = plot(cell_id',gene,'Color',Fig.Color);
    set(gca,'ytick',[]);
    set(gca,'xtick',[]);
    xlabel('Ordered Cells', 'FontSize', 24);
    ylabel('Gene Expression', 'FontSize', 24);
end
end