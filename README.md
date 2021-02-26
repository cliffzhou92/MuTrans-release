# MuTrans (The MUltiscale method for TRANSient cells )

**Dissecting Transition Cells from Single-cell Transcriptome through Multi-Scale Stochastic Dynamics**

Peijie Zhou, Shuxiong Wang, Tiejun Li, Qing Nie


We proposed a method based on the multiscale technique for stochastic dynamical systems to analyze single-cell transcriptome data and identify transition cells. MuTrans naturally connects the languages of dynamical system with single-cell data analysis to describe cell-fate transitions, i.e. Attractor Basins with Meta-stable States, Saddle points with Transient States, and Most Probable Paths with Cell Lineages.

<img src="https://github.com/cliffzhou92/MuTrans-release/tree/main/img/algorithm.png" width="600">

## Installation
### In Matlab:
To display gene expression matrix by Transcendental, one need the 'smooth' function (https://www.mathworks.com/help/curvefit/smoothing.html) from Curve Fitting toolbox of Matlab.

### In Python:
1. install the matlab engine API for python, [instructions here](https://www.mathworks.com/help/matlab/matlab_external/install-the-matlab-engine-for-python.html)
2. install the dependency package [PyEMMA](http://emma-project.org/latest/INSTALL.html).
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
Lineage = InferLineage(Output,root,par)

%% Infer the cell lineage based on multi-scale analysis results

%% Output: the output object from DynamicalAnalysis function

%% root: the root state of Lineage

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



## Example Notebooks
**System** | **Data Source** | **Notebook File** | **PDF File**
------------| -------------- | ------------|------------
Saddle-Node Bifurcation | Simulation Data in this study | here | here
Epithelial-Mesenchymal Transition |[Pastushenko et al.](https://www.nature.com/articles/s41586-018-0040-3)|here |here
iPSC Reprogramming |[Bargaje et al.](https://www.pnas.org/content/early/2017/02/03/1621412114)|here|here
Myelopoiesis |[Olsson et al.](https://www.nature.com/articles/nature19348) |here |here
Lymphoid Lineage Blood Differentiation | [Herman et al.](https://www.nature.com/articles/nmeth.4662) |here|here
