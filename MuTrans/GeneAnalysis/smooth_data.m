function data_smooth = smooth_data(data, span , method)
[N_cell,N_genes] = size(data);
data_smooth = nan(N_cell,N_genes);

for k = 1: N_cell
    data_smooth(k,:) = smooth(data(k,:),span, method);
end

end
