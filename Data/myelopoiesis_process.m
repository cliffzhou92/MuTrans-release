clear all;
load Olsson_data.mat % the raw count matrix
text = data1.textdata;
gene_name_all = string(text(2:end,1));
data_all = data1.data';

load olsson.mat % ICGS matrix
gene_id = ismember(gene_name_all,gene_name);
data_olsson = data_all(:,gene_id);
gene_name_select = gene_name_all(gene_id);


load Olsson_origin_labs.txt %original annotations
labs_cluster = Olsson_origin_labs;
load Olsson_SP_labs.txt
labs_type = Olsson_SP_labs;


data = log(data_olsson+1);
%% tsne plot
rng(1)
figure;
%[~,score] = pca(data);
score = tsne(data,'Algorithm','exact','Distance','cosine');
gscatter(score(:,1),score(:,2),labs_cluster)

%% remove apparent outliers

rm_cell = [94,104,124,97,129,107,221];
N_cell = size(data_olsson,1);
all_id = 1:N_cell;
rm_id = ismember(all_id,rm_cell);
keep_id = ~rm_id;

gscatter(score(:,1),score(:,2),keep_id)


%%
data_olsson_clean = data_olsson(keep_id,:);
labs_cluster_clean = labs_cluster(keep_id);
labs_type_clean = labs_type(keep_id);


data = log(data_olsson_clean+1);
figure;
%[~,score] = pca(data);
score = tsne(data,'Algorithm','exact','Distance','correlation');
gscatter(score(:,1),score(:,2),labs_type_clean)

figure;
gscatter(score(:,1),score(:,2),labs_cluster_clean)

save('olsson_new.mat','data_olsson_clean','labs_cluster_clean','labs_type_clean','gene_name_select');
