#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python wrapper functions for MuTrans

@author: cliffzhou
"""

import matlab.engine
import os, sys
from pathlib import Path

eng = matlab.engine.start_matlab()
sys.path.append(os.path.join(Path(os.getcwd()).parent, "MuTrans"))
eng.eval("addpath(genpath('../'))")

import numpy as np
import pyemma.msm as msm
import pyemma.plots as mpl
import matplotlib.pyplot as plt
import scipy
import networks as nw
import seaborn as sns

def dynamical_analysis(sc_object,par):
    """ MuTrans dynamical system analysis by multi-scale reduction.
    
    Arguments: 
        
        sc_object: the anndata object storing single-cell datasets, including the expression matrix and the low-dimensional coordinates (umap or tsne) 
        
        par: parameters in MuTrans, stored in dictionary, with important keys:
        
        ### must provide:
        "K_cluster": the number of clusters/attractors. The eigen-peak can serve as the reference, and it's suggested to do Louvain clustering or other biological analysis for the prior knowledge.
        
        ### optional with defaults
        
        [parameters in construcing cell-cell random walk]
        "perplex": the parameter controlling the "diffusiveness" of cell-cell random walk, similar to the parameter in tSNE. Control the local bandwith of Gaussian kernels. Larger perplex represents larger bandwith and more connectiveness between cells, and smaller perplex may reduce "shortcuts" in lineage analysis. Default is N/2.5. Suggest to change from N/5 to N/2.
        "choice_distance": the choice of distance when constructing the Gaussian kernel. Can be "euclid"(default) or "cosine". This is similar to the "metric" parameter in UMAP.
        "weight_scale": When construcing the cell-cell weights by Gaussian kernels, whether adopting the strategy in tSNE to scale the Gaussian kernels by row sums. Default is "True" and suggested for large-scale 10X data to stabilize. 
        "alpha_scale": When construcing the cell-cell weights by Gaussian kernels, the power row sums to scale the Gaussian kernels by adopting the strategy in diffusion map. Default is 0.
        
        [parameters in multiscale reduction]
        "reduce_large_scale": whether to boost the analysis by the DECLARE dynamical-preserving pre-processing step. Default is False, recommended True for large-scale datasets with more than 5000 cells or more than 10 attractors, whenever the default analysis is slow.
        "reduce_num_meta_cell": number of meta-cells or microscopic states in the DECLARE reduction. Only valid when "reduce_large_scale" is True.
        "trials": Because MuTrans uses the EM algorithm for attractor assignment, the results are local minimums and subject to initial assignments. The number of trials allow run the assignments for multiple times, and select the assignment with minimum objective functions. Increasing the trails will increase the running time, while increase the robustness of results.
        "initial": The initial cluster assignment in each trial. Can be "random" (random assignment -- fastest),"pca" (PCA+kmeans) or tsne (TSNE+kmeans). The default is "tsne".
        
        [parameters in dynamical manifold]
        "reduction_coord": The coordinates provided as the basis to construct dynamical manifold. Can be "tsne" or "umap", and the coordinates should be stored in the "obsm" attribute of anndata object. 
        
    Output: return the anndata object, with MuTrans outputs stored in the "da_out" and "land" key in the "uns" attribute of anndata object. Can be loaded back into matlab for gene analysis.
        
    """
    data = sc_object.X
    data_convert = matlab.double(data.tolist())
    
    if par.__contains__('perplex'):
        par['perplex'] = float(par['perplex'])
    
    if par.__contains__('K_cluster'):
        par['K_cluster'] = float(par['K_cluster'])
        
    out = eng.DynamicalAnalysis(data_convert,par)
    perm = np.asarray(out['perm_class']).ravel()
    perm = perm.astype(int)-1
    
    if (par.__contains__('reduction_coord') == False) or par['reduction_coord']== 'tsne':
        X_embedding = sc_object.obsm['X_tsne']
    elif par['reduction_coord']== 'umap':
        X_embedding = sc_object.obsm['X_umap']
    
    X_embedding = X_embedding- np.mean(X_embedding, axis =0)
    
    if (par.__contains__('reduce_large_scale') == False) or par['reduce_large_scale']== False :
        embedding_r = X_embedding[perm]
        out["embedding_2d"] = matlab.double(embedding_r.tolist())
    else:
        out["embedding_2d"] = matlab.double(X_embedding.tolist())
    
    land = eng.ConstructLandscape(out,par)
    sc_object.uns['da_out']=out
    sc_object.uns['land']=land
    
    if par.__contains__('write_anndata') and par['write_anndata']:
        ind = np.argsort(perm)
        sc_object.obs['land'] = np.asarray(land['land_cell'])[ind]
        sc_object.obs['entropy'] = np.asarray(out['H']).ravel()[ind]
        label = np.asarray(out['class_order']).astype(int)-1
        sc_object.obs['attractor'] = label.astype(str).ravel()[ind]
        sc_object.obsm['trans_coord'] = np.asarray(land['trans_coord'])[ind]
        sc_object.obsm['membership'] = np.asarray(out['rho_class'])[ind]
    
    
    return sc_object

def plot_landscape(sc_object,alpha_land = 0.5, show_colorbar = False):
    land = sc_object.uns['land']
    land_value = np.asarray(land['land']).T
    x_grid = np.asarray(land['grid_x']).ravel()
    y_grid = np.asarray(land['grid_y']).ravel()
    plt.contourf(x_grid, y_grid, land_value, levels=14, cmap="Greys_r",zorder=-100, alpha = alpha_land)
    if show_colorbar:
        plt.colorbar()

def plot_cluster_num (sc_object, par = None, k_plot = 20, order =2.0):
    data = sc_object.X
    data_convert = matlab.double(data.tolist())
    if par == None:
        par = {}
        par ['choice_distance'] = 'euclid'
    
    par['order'] = order
    cluster_num_est = eng.EstClusterNum(data_convert,par)
    plt.plot(np.arange(k_plot)+1,np.asarray(cluster_num_est['ratio'])[:k_plot])
    
    return cluster_num_est

    
    
def infer_lineage(sc_object,si=0,sf=1,method = 'MPFT',flux_fraction = 0.9, size_state = 0.1, size_point = 50, alpha_land = 0.5, alpha_point = 0.5, size_text=20, show_colorbar = False, color_palette = None):
    """ Infer the lineage by MPFT or MPPT approach.
    
    Arguments: 
        
        sc_object: the anndata object storing single-cell datasets and MuTrans output (after running the dynamical_analysis function).
        
        method: method to infer the lineage, can be "MPFT" (Maximum Probability Flow Tree, global structure) or "MPPT" (Most Probable Path Tree, local transitions).
        
        si and sf: starting and targeting attractors in MPPT method for transition paths analysis.
        
        flux_fraction: the cumulative percentage of top transition paths probability flux shown on the figure.
        
    """
    projection = np.asarray(sc_object.uns['land']['trans_coord'])
    K = sc_object.uns['da_out']['k']
    K = int(K)
    P_hat = np.asarray(sc_object.uns['da_out']['P_hat'])
    mu_hat = np.asarray(sc_object.uns['da_out']['mu_hat'])
    
    centers = []
    labels = np.asarray(sc_object.uns['da_out']['class_order']).ravel()
    labels = labels.astype(int)-1
    
    if color_palette is None:
        color_palette = sns.color_palette('Set1', K)
    cluster_colors = [color_palette[x] for x in labels]
    
    
    for i in range(K):
        index = labels==i
        p = np.mean(projection[index], axis=0)
        centers.append(p)
    centers = np.array(centers)

    if method == 'MPFT':
        Flux_cg = np.diag(mu_hat.reshape(-1)).dot(P_hat)
        max_flux_tree = scipy.sparse.csgraph.minimum_spanning_tree(-Flux_cg)
        max_flux_tree = -max_flux_tree.toarray()
        for i in range(K):
            for j in range(i+1,K):   
                max_flux_tree[i,j]= max(max_flux_tree[i,j],max_flux_tree[j,i])
                max_flux_tree[j,i] = max_flux_tree[i,j]
         
        nw.plot_network(max_flux_tree, pos=centers, state_scale=size_state, state_sizes=mu_hat, arrow_scale=2.0,arrow_labels= None, arrow_curvature = 0.2, ax=plt.gca(),max_width=1000, max_height=1000)
        plot_landscape(sc_object, alpha_land=alpha_land,show_colorbar = show_colorbar)
        plt.scatter(*projection.T, s=size_point, linewidth=0, c=cluster_colors, alpha=alpha_point)
        
    if method == 'MPPT':
        M = msm.markov_model(P_hat)
        
        #state_reorder = np.array(range(K))
        #state_reorder[0] = si
        #state_reorder[-1] = sf
        #state_reorder[sf+1:-1]=state_reorder[sf+1:-1]+1
        #state_reorder[1:si]=state_reorder[1:si]-1
        tpt = msm.tpt(M, [si], [sf])
        Fsub = tpt.major_flux(fraction=flux_fraction)
        Fsubpercent = 100.0 * Fsub / tpt.total_flux
        
       
        plot_landscape(sc_object,alpha_land=alpha_land, show_colorbar = show_colorbar)
        nw.plot_network(Fsubpercent, state_scale=size_state,pos=centers, arrow_label_format="%3.1f",arrow_label_size = size_text,ax=plt.gca(), max_width=1000, max_height=1000)
        plt.scatter(*projection.T, s=size_point, linewidth=0, c=cluster_colors, alpha=alpha_point)
        plt.axis('off')        
                
                
                
                