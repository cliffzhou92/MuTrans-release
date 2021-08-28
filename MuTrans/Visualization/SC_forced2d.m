function mappedX = SC_forced2d (P,par)
dims = par.dims;
k = par.isomap_kn;
[phi,~] = eigs(P,dims+1);
mappedX = isomap(phi(:,2:end), 2, k);
end