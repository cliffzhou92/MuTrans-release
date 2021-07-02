# MuTrans (The MUltiscale method for TRANSient cells )

**Dissecting Transition Cells from Single-cell Transcriptome through Multi-Scale Stochastic Dynamics**

Peijie Zhou, Shuxiong Wang, Tiejun Li, Qing Nie


We proposed a method based on the multiscale technique for stochastic dynamical systems to analyze single-cell transcriptome data and identify transition cells. MuTrans naturally connects the languages of dynamical system with single-cell data analysis to describe cell-fate transitions, i.e. Attractor Basins with Meta-stable States, Saddle points with Transient States, and Most Probable Paths with Cell Lineages.

<img src="https://github.com/cliffzhou92/MuTrans-release/blob/main/img/figure1-01.png" width="700">

## Installation
### In Matlab:
The code has been tested in Matlab R2019b and R2020a. To display gene expression matrix by Transcendental, one need the 'smooth' function (https://www.mathworks.com/help/curvefit/smoothing.html) from Curve Fitting toolbox of Matlab.

### In Python:
1. Install the matlab engine API for python, [instructions here](https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html)
2. Install the dependency package [PyEMMA >=2.5.6](http://emma-project.org/latest/INSTALL.html), [Scanpy](https://scanpy.readthedocs.io/en/stable/installation.html), Numpy, Pandas and Seaborn.
3. cd to the ``./Example/`` folder and analysis in Jupyter notebook

## Basic Usage
### In Matlab:

```
Output = DynamicalAnalysis(data,par)

%% Run the multi-scale analysis, returns the attractor basins, coarse-grained transition probabilities and membership assignment

%% data: the pre-processed single-cell gene expression, with N_cells x N_genes

%% par: the adjustable parameters
```

```
Lineage = InferLineage(Output,par)

%% Infer the cell lineage based on multi-scale analysis results

%% Output: the output object from DynamicalAnalysis function

%% par: the adjustable parameters

```

```
Land = ConstructLandscape(Output,par)

%% Construct and visualize the dynamical manifold

%% Output: the output object from DynamicalAnalysis function

%% par: the adjustable parameters

```

```
Genes = GeneAnalysis(i,j,Output,par)

%% Transition cell and gene analysis from state transition i to j.

%% i and j: the starting and targeting states for analysis

%% Output: the output object from DynamicalAnalysis function

%% par: the adjustable parameters

```

### In Python:

```
import pyMuTrans as pm
out = pm.plot_cluster_num(adata, par, k_plot= 10) # use EPI to determine number of clusters
par["reduce_large_scale"] = True # optional, to use DECLARE module speeding-up the calculation
par["reduce_num_meta_cell"] = 1500 # optional, to set the number of microsopic meta-stable states in DECLARE
adata = pm.dynamical_analysis(adata,par) # MuTrans Analysis on Anndata obejct
pm.infer_lineage(adata,si=2,sf=0,method = "MPPT",size_point =40, size_text = 10,alpha_point = 0.5) # plot the transition trajectory on dynamical manifold
```


## Example Notebooks
**System** | **Data Source** | **Notebook File** | **PDF File**
------------| -------------- | ------------|------------
Saddle-Node Bifurcation | Simulation Data in this study | [example_saddle_node.mlx](https://github.com/cliffzhou92/MuTrans-release/blob/main/Example/example_saddle_node.mlx) | [example_saddle_node.pdf](https://github.com/cliffzhou92/MuTrans-release/blob/main/Example/example_saddle_node.pdf)
Epithelial-Mesenchymal Transition |[Pastushenko et al.](https://www.nature.com/articles/s41586-018-0040-3)|[Matlab:example_emt.mlx](https://github.com/cliffzhou92/MuTrans-release/blob/main/Example/example_emt.mlx); [Python:example-emt.ipynb](https://github.com/cliffzhou92/MuTrans-release/blob/main/Example/example-emt.ipynb) |[Matlab: example_emt.pdf](https://github.com/cliffzhou92/MuTrans-release/blob/main/Example/example_emt.pdf); [Python: example-emt.pdf](https://github.com/cliffzhou92/MuTrans-release/blob/main/Example/example-emt.pdf)
iPSC Reprogramming |[Bargaje et al.](https://www.pnas.org/content/early/2017/02/03/1621412114)|[example_ipsc.mlx](https://github.com/cliffzhou92/MuTrans-release/blob/main/Example/example_ipsc_early.mlx) | [example_ipsc.pdf](https://github.com/cliffzhou92/MuTrans-release/blob/main/Example/example_ipsc_early.pdf)
Myelopoiesis |[Olsson et al.](https://www.nature.com/articles/nature19348) |[example_olsson.mlx](https://github.com/cliffzhou92/MuTrans-release/blob/main/Example/example_olsson.mlx) | [example_olsson.pdf](https://github.com/cliffzhou92/MuTrans-release/blob/main/Example/example_olsson.pdf)
Lymphoid Lineage Blood Differentiation | [Herman et al.](https://www.nature.com/articles/nmeth.4662) |[example_mpp.mlx](https://github.com/cliffzhou92/MuTrans-release/blob/main/Example/example_mpp.mlx) | [example_mpp.pdf](https://github.com/cliffzhou92/MuTrans-release/blob/main/Example/example_mpp.pdf)
Human Bone Marrow | [Setty et al.](https://www.nature.com/articles/s41587-019-0068-4) | [example_bone_marrow-new.ipynb](https://github.com/cliffzhou92/MuTrans-release/blob/main/Example/example_bone_marrow-new.ipynb) | [example_bone_marrow.pdf](https://github.com/cliffzhou92/MuTrans-release/blob/main/Example/example_bone_marrow-new.pdf)
Blood Development in Mouse Gastrulation | [Pijuan-Sala et al.](https://www.nature.com/articles/s41586-019-0933-9) | [example_haem_development_15K.ipynb](https://github.com/cliffzhou92/MuTrans-release/blob/main/Example/example_haem_development_15K.ipynb) | [example_haem_development_15K.pdf](https://github.com/cliffzhou92/MuTrans-release/blob/main/Example/example_haem_development_15K.pdf)
