#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python wrapper functions for MuTrans

@author: cliffzhou
"""

import matlab.engine
eng = matlab.engine.start_matlab()
import os, sys
sys.path.append(os.path.join(os.getcwd(), "DynamicalAnalysis"))
eng.eval("addpath(genpath('./DynamicalAnalysis'))")
eng.eval("addpath(genpath('./tsne'))")
eng.eval("addpath(genpath('./Visualization'))")
eng.eval("addpath(genpath('./Brewermap'))")

import numpy as np
import pyemma.msm as msm
import pyemma.plots as mpl
import matplotlib.pyplot as plt
import scipy
import networks as nw

def dynamical_analysis(sc_object,par):
    
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
    
    
    if (par.__contains__('reduce_large_scale') == False) or par['reduce_large_scale']== False :
        embedding_r = X_embedding[perm]
        out["embedding_2d"] = matlab.double(embedding_r.tolist())
    else:
        out["embedding_2d"] = matlab.double(X_embedding.tolist())
    
    land = eng.ConstructLandscape(out,par)
    sc_object.uns['da_out']=out
    sc_object.uns['land']=land
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

    
    
def infer_lineage(sc_object,si=0,sf=1,method = 'MPFT',flux_fraction = 0.9, size_state = 0.1, size_point = 50, alpha_land = 0.5, alpha_point = 0.5, size_text=20, show_colorbar = False):
    projection = np.asarray(sc_object.uns['land']['trans_coord'])
    K = sc_object.uns['da_out']['k']
    K = int(K)
    P_hat = np.asarray(sc_object.uns['da_out']['P_hat'])
    mu_hat = np.asarray(sc_object.uns['da_out']['mu_hat'])
    
    centers = []
    labels = np.asarray(sc_object.uns['da_out']['class_order']).ravel()
    labels = labels.astype(int)-1
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
    
        
        state_reorder = np.array(range(K))
        state_reorder[0] = si
        state_reorder[-1] = sf
        state_reorder[sf+1:-1]=state_reorder[sf+1:-1]+1
        state_reorder[1:si]=state_reorder[1:si]-1
        tpt = msm.tpt(M, [si], [sf])
        Fsub = tpt.major_flux(fraction=flux_fraction)
        Fsubpercent = 100.0 * Fsub / tpt.total_flux
        
       
        plot_landscape(sc_object,alpha_land=alpha_land, show_colorbar = show_colorbar)
        nw.plot_network(Fsubpercent, state_scale=size_state,pos=centers, arrow_label_format="%3.1f",arrow_label_size = size_text,ax=plt.gca(), max_width=1000, max_height=1000)
        plt.scatter(*projection.T, s=size_point, linewidth=0, c=cluster_colors, alpha=alpha_point)
        plt.axis('off')        
                
                
                
                