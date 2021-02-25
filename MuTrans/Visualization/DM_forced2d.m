function mappedX = DM_forced2d (data,par)
dims = par.dims;
sigma = par.sigma;
k = par.isomap_kn;
[phi,~] = diffusionmap(data,'loc',sigma,dims);
mappedX = isomap(phi, 2, k);
end