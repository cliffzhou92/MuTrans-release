function theta = logistic_transition_score (cluster_id1, cluster_id2, class_order, rho_class, par)

g1 = cluster_id1;
g2 = cluster_id2;

%% pre-processing the data
id_data = ismember(class_order,[g1 g2]);
class_order_ini = class_order(id_data);
id_type = ismember(1:max(class_order),[g1 g2]);
id_ntype = ~id_type;
rho_other = rho_class(id_data,id_ntype);


% only considers the cells transit between two states (useful if g1 is the bifurcation state)
id_data_keep = max(rho_other,[],2)<par.otherkeep;
rho_ini = rho_class(id_data,:);
rho = rho_ini(id_data_keep,:);
class_order_keep = class_order_ini(id_data_keep);

% two groups of data
bifurid_g1 = class_order_keep == g1;
bifurid_g2 = class_order_keep == g2;
rho_g1 = rho(bifurid_g1,:);
rho_g2 = rho(bifurid_g2,:);

score_g1 = rho_g1(:,g1)./(rho_g1(:,g1)+rho_g1(:,g2));
score_g2 = (rho_g2(:,g1))./(rho_g2(:,g1)+rho_g2(:,g2));
[score1,~] = sort(score_g1,'descend');
[score2,~] = sort(score_g2,'descend');
score_logistic = [score1;score2];

%% logistic transition
% data normalization

%sort the genes

L = length(score_logistic);
h = 1/(L-1);
x_fit = 0:h:1;
[~, p_logistic, ~, ~] = fit_logistic(x_fit',score_logistic);
theta = abs(p_logistic(3));

end