function data_scale = normalize_data (data)
idata = data;
center = mean(idata,1);
scale = std(idata, 0,1);
tscale = scale;
%=Check for zeros and set them to 1 so not to scale them.
scale(tscale == 0) = 1;
%== Center and scale the data
idata = bsxfun(@minus, idata, center);
data_scale = bsxfun(@rdivide, idata, scale);
end