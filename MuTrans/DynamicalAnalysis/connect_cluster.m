function out = connect_cluster(id1, id2 ,class_order, rho_class, thresh)
    rho_12 = rho_class(class_order == id1, id2);
    rho_21 = rho_class(class_order == id2, id1);
    rho_combo = [rho_12;rho_21];
    out = max(rho_combo)> thresh;
end