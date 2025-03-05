# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 19:24:03 2024

@author: Wen
"""

get_ipython().run_line_magic('reset', '-sf')
style_path = r'C:\Wen\OneDrive\CodeHub\Python\style_paper.mplstyle'

import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
# import h5py
import pickle as pkl
import scipy.io as sio

from matplotlib.patches import Rectangle
from matplotlib import ticker

from utils import fun_pev
from utils import fun_loadData

#%% load SFC 
path = r'C:\Wen\Project\Osc_control\Data\monkey_combined\SFC_examlple_over_time'
file = '\\Ocean_SFC_over_time_example_mats2' 

data = sio.loadmat(path+file)

# sfc = data['SFC_example_pT1_0710_ch156_in_256'][:,3:]# example 1(backup)
sfc = data['SFC_example_pT1_0710_ch130_in_256'][:,3:] # example 1
# sfc = data['SFC_example_pT2_0710_ch64_in_256'][:,3:] # example 1
#%%
s_onset = [3,14]
with plt.style.context(style_path):
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"][3:]
    fig, axs = plt.subplots(1, 1, figsize=(1.3, 4), dpi=300) # , constrained_layout=True
    plt.subplots_adjust(wspace=0.05, hspace=0.2)
    
    cmap = plt.cm.inferno
    h1 = axs.pcolormesh(sfc, cmap=cmap, vmin=0.03, vmax=0.25)
    # h1 = axs.pcolormesh(sfc, cmap=cmap, vmin=0.03, vmax=0.15)
    axs.vlines(s_onset, 0, sfc.shape[0], color='w', linestyle='--')
    
    axs.set_xticks(s_onset, ['S1','S2'])
    axs.set_ylabel('Channel')
    cbar = fig.colorbar(h1, aspect=50)

#%% plot fw 排序
path = r'C:\Wen\Project\Osc_control\Data\Result_ocean\SPK\tdr_frontal_L2_croval\tdr_mem_delay1'
file = '\\0710_weight_mem_delay.mat' 

# fw1_example_neuidx = 95 
fw_example_neuidx = 108 # example 1
# fw_example_neuidx = 124 # example 1(backup)
# fw_example_neuidx = 53 # example 2

fw1 = sio.loadmat(path+file)['r1_weight']
fw2 = sio.loadmat(path+file)['r2_weight']
    
with plt.style.context(style_path):
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    fig, axs = plt.subplots(1, 1, figsize=(1.5, 1), dpi=300) # , constrained_layout=True

    data = {'values': fw1[:,0], 'neuidx':np.arange(len(fw1))}
    df = pd.DataFrame(data)
    # 按照'values'列从大到小排序
    sorted_df = df.sort_values(by='values', ascending=False).reset_index()
    idx = sorted_df.query('neuidx == @fw_example_neuidx').index
    
    axs.scatter(np.arange(len(fw1)), sorted_df['values'], color=colors[0],s=2)
    axs.scatter(idx, sorted_df.loc[idx, 'values'], color='k', s=40)
    
    axs.set_xticks([10, 60, 110])
    axs.set_xlabel('Unit')
    axs.set_ylabel('Feature weight')
    

with plt.style.context(style_path):
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    fig, axs = plt.subplots(1, 1, figsize=(1.5, 1), dpi=300) # , constrained_layout=True

    data = {'values': fw2[:,0], 'neuidx':np.arange(len(fw2))}
    df = pd.DataFrame(data)
    # 按照'values'列从大到小排序
    sorted_df = df.sort_values(by='values', ascending=False).reset_index()
    idx = sorted_df.query('neuidx == @fw_example_neuidx').index
    
    axs.scatter(np.arange(len(fw1)), sorted_df['values'], color=colors[1],s=2)
    axs.scatter(idx, sorted_df.loc[idx, 'values'], color='k', s=40)
    
    axs.set_xticks([10, 60, 110])
    axs.set_xlabel('Unit')
    axs.set_ylabel('Feature weight')



    