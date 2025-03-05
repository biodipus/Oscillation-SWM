# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 16:05:22 2023

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
# %%
fpath = r'C:\Wen\Project\Osc_control\Data\Data_ocean\SPK\FR_segment_win100step50'
savepath = r'C:\Wen\Project\Osc_control\Data\Result_ocean\SPK\PEV_seqL2'
days = ['0729']
seqL = 2
#%% load spk data
fr_all = [None] * len(days)
ch_info_all = [None] * len(days)

select_rule = 0
region = ['PFC','PM']#, 'PM', 'PAR', 'M1']
for i in range(len(days)):
    file = fpath + '\\'+days[i]+'_fr_L2_win100step50_selChs.pkl'
    info, neu_sit_info, fr_arr, mkts_merge, spk_arr, mktime = fun_loadData.load_fr_data(file, region, select_rule)
    # info, fr_arr = fun_loadData.creat_psudo_trial(info, fr_arr) # creat psudo trials
    fr_all[i] = fr_arr # neuron*trial*time
    ch_info_all[i] = neu_sit_info

if len(fr_all) > 1:
    fr_all = np.concatenate(fr_all, axis=0)
    ch_info_all = pd.concat(ch_info_all).reset_index()
else:    
    fr_all = fr_all[0]
    ch_info_all = ch_info_all[0]
#%% plot raster
figpath = r'C:\Wen\Project\Osc_control\Data\Result_ocean\SPK\raster_plots'
for neuID in range(59, 60): #range(fr_all.shape[0]): (9,10)
    print(neuID)
    with plt.style.context(style_path):
        colors = plt.rcParams["axes.prop_cycle"].by_key()["color"][3:]
        fig, axs = plt.subplots(2, 3, figsize=(2, 1.5), dpi=300, gridspec_kw={'width_ratios': [3,3,2]}) # , constrained_layout=True
        plt.subplots_adjust(wspace=0, hspace=0.2)
        k = 0
        for i in range(6):
            linewidth = 0.3
            idx = info.loc[:,'Target_2'] == i+1
            x = spk_arr[neuID, idx, 750:1500] # neuron*trial*time
            axs[0,0].get_xaxis().set_visible(False)
            axs[0,0].set_ylabel('Trial #')
            axs[0,0].add_patch(Rectangle((250, 0), 250, spk_arr.shape[1], edgecolor=None, facecolor="gray", alpha=0.05))
            for t in range(x.shape[0]):
                for j in range(x.shape[1]):
                    if(x[t,j] == 1):
                        x1 = [t+k, t+k+0.5]
                        x2 = [j, j]
                        axs[0,0].plot(x2,x1,color=colors[i], linewidth=linewidth)
            
            x = spk_arr[neuID, idx, 1500:2250]
            axs[0,1].get_xaxis().set_visible(False)
            axs[0,1].get_yaxis().set_visible(False)
            axs[0,1].spines['left'].set_visible(False)
            axs[0,1].add_patch(Rectangle((250, 0), 250, spk_arr.shape[1], edgecolor=None, facecolor="gray", alpha=0.05))
            for t in range(x.shape[0]):
                for j in range(x.shape[1]):
                    if(x[t,j] == 1):
                        x1 = [t+k, t+k+0.5]
                        x2 = [j, j]
                        axs[0,1].plot(x2,x1,color=colors[i], linewidth=linewidth)
                        
            x = spk_arr[neuID, idx, 2250:2750]
            axs[0,2].get_xaxis().set_visible(False)
            axs[0,2].get_yaxis().set_visible(False)
            axs[0,2].spines['left'].set_visible(False)
            for t in range(x.shape[0]):
                for j in range(x.shape[1]):
                    if(x[t,j] == 1):
                        x1 = [t+k, t+k+0.5]
                        x2 = [j, j]
                        axs[0,2].plot(x2,x1,color=colors[i], linewidth=linewidth)
            k = k+t    
        
        for i in range(6):
            idx = info.loc[:,'Target_1'] == i+1
            x = fr_arr[neuID, idx, 15:30].mean(0)
            axs[1,0].plot(x, color=colors[i])
            
            x = fr_arr[neuID, idx, 30:45].mean(0)
            axs[1,1].plot(x, color=colors[i])
        
            x = fr_arr[neuID, idx, 45:55].mean(0)
            axs[1,2].plot(x, color=colors[i])
        
        ymin = 0
        ymax = 48
        for i in range(3):
            axs[1,i].set_ylim([ymin, ymax])
            
        axs[1,0].add_patch(Rectangle((5, ymin), 4.3, ymax, edgecolor=None, facecolor="gray", alpha=0.15))
        axs[1,0].set_xticks([5])
        axs[1,0].set_xticklabels(['S1'])
        axs[1,0].set_ylabel('FR (Hz)')
        
        axs[1,1].get_yaxis().set_visible(False)
        axs[1,1].spines['left'].set_visible(False)
        axs[1,1].add_patch(Rectangle((5, ymin), 4.3, ymax, edgecolor=None, facecolor="gray", alpha=0.15))
        axs[1,1].set_xticks([5])
        axs[1,1].set_xticklabels(['S2'])
        
        axs[1,2].get_yaxis().set_visible(False)
        axs[1,2].spines['left'].set_visible(False)
        axs[1,2].set_xticks([3.75])
        axs[1,2].set_xticklabels(['Delay'])
        # plt.suptitle(days[0] + ' ' + str(neuID))
        plt.show()
        
        # fig.savefig(figpath + '\\S1_' + str(neuID) + '_' + days[0] + '.png', dpi='figure', bbox_inches='tight')

#%% plot example raw data with LFP and SPK
import h5py

fpath = r'D:\Osc_control\Data\ocean\LFP\Preprocessed\trial_cuts'
days = ['0713']

infoPath = r'C:\Wen\Projects\Osc_control\Data\Data_ocean\SPK\sorting' 

sfreq = 1000
info = [None]*len(days)
lfp_raw = [None]*len(days)
mk_time = [None]*len(days)
bad_trial_idx = [None]*len(days)
for i in range(len(days)):
    day = days[i]
    lfppath = os.path.join(fpath, 'ocean2021-' + day + '_CueFrameTrialData_LB1500_RB5500_denoise.mat')
    lfp_raw[i] = h5py.File(lfppath)['LFP_trial_data']
    mk_time[i] = np.array(h5py.File(lfppath)['trial_events']).T
mk_time = np.concatenate(mk_time)

lfp_raw = np.concatenate(lfp_raw, axis=0)
lfp_raw = np.transpose(lfp_raw, [0,2,1])

with plt.style.context(style_path):
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"][0:]
    fig, axs = plt.subplots(1, 1, figsize=(1, 0.6), dpi=300) # , constrained_layout=True
    axs.plot(np.arange(500), lfp_raw[10,26,2000:2500], linewidth=0.3, color='k')
    axs.get_xaxis().set_visible(False)
    axs.get_yaxis().set_visible(False)
    axs.spines['left'].set_visible(False)
    axs.spines['bottom'].set_visible(False)
    
    figpath = r'C:\Wen\Project\Replay\Result_ocean\plot_example\figure'
    # fig.savefig(figpath + '\\Raw_example.pdf', dpi=300, bbox_inches='tight')
#%% plot power example 
lfppath = r'C:\Wen\Project\Osc_control\Data\Result_ocean\plot_example(neuron_LFP)_data'
# file = '0710_WB103_TH_mats'
# file = '0710_WB57_TH_HG_mats'# gamma memory
# file = '0710_WB32_TH_HG_mats' # gamma entry
file = '0729_WB70_TH_HG_mats' # gamma s1
lfpdata = sio.loadmat(lfppath+'\\'+file+'.mat')

# fr_arr = np.concatenate([lfpdata['TH_hold_T1_bins'],lfpdata['TH_hold_T2_bins'],lfpdata['TH_hold_Rule_bins']], axis=1)
fr_arr = np.concatenate([lfpdata['HG_hold_T1_bins'],lfpdata['HG_hold_T2_bins'],lfpdata['HG_hold_Rule_bins']], axis=1)
fr_arr = np.concatenate([fr_arr[:,:32], fr_arr[:,-10:]], axis=1)
with plt.style.context(style_path):
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"][3:]
    fig, axs = plt.subplots(2, 3, figsize=(2, 1.5), dpi=800, gridspec_kw={'width_ratios': [3,3,2]}) # , constrained_layout=True
    plt.subplots_adjust(wspace=0.0, hspace=0.2)
    
    #### plot TF
    k = 0
    fr_arr_sort = np.zeros(fr_arr.shape)
    for i in range(6):
        idx = info.loc[:,'Target_2'] == i+1
        t = idx.sum()
        fr_arr_sort[k:k+t, :] = fr_arr[idx, :]
        # axs[0,0].plot([16.5,16.5],[k+10,k+t-10], color=colors[i], linewidth=2)
        # axs[0,1].plot([16.5,16.5],[k+10,k+t-10], color=colors[i], linewidth=2)
        axs[0,2].plot([10.5,10.5],[k+10,k+t-10], color=colors[i], linewidth=2)
        
        axs[0,0].plot([0, 16],[k+t, k+t], color='w', linewidth=0.3)
        axs[0,1].plot([0, 16],[k+t, k+t], color='w', linewidth=0.3)
        axs[0,2].plot([0, 10.5],[k+t, k+t], color='w', linewidth=0.3)
        k = k + t
        
    # vmin, vmax = 0.019, 0.023
    vmin, vmax = 0.008, 0.010
    # vmin, vmax = 0.01, 0.02
    cmap = plt.cm.gray
    axs[0,0].pcolormesh(fr_arr_sort[:, 0:16], cmap = cmap, vmin=vmin, vmax=vmax)
    axs[0,0].get_xaxis().set_visible(False)
    axs[0,0].set_ylabel('Trial')
    axs[0,0].add_patch(Rectangle((5.5, 0), 5, fr_arr.shape[0], edgecolor=None, facecolor="gray", alpha=0.25))
    
    axs[0,1].pcolormesh(fr_arr_sort[:, 16:32], cmap = cmap, vmin=vmin, vmax=vmax)
    axs[0,1].get_xaxis().set_visible(False)
    axs[0,1].get_yaxis().set_visible(False)
    axs[0,1].spines['left'].set_visible(False)
    axs[0,1].add_patch(Rectangle((5.5, 0), 5, fr_arr.shape[0], edgecolor=None, facecolor="gray", alpha=0.25))

    axs[0,2].pcolormesh(fr_arr_sort[:, 32:], cmap = cmap, vmin=vmin, vmax=vmax)
    axs[0,2].get_xaxis().set_visible(False)
    axs[0,2].get_yaxis().set_visible(False)
    axs[0,2].spines['left'].set_visible(False)
    
    #### plot average lines                
    for i in range(6):
        idx = info.loc[:,'Target_2'] == i+1
        x = fr_arr[idx, 0:16].mean(0)
        axs[1,0].plot(x, color=colors[i])

        x = fr_arr[idx, 16:32].mean(0)
        axs[1,1].plot(x, color=colors[i])
        
        x = fr_arr[idx, 32:].mean(0)
        axs[1,2].plot(x, color=colors[i])
    
    ylim = [0.0055, 0.009]
    # ylim = [0.0055, 0.0085]
    # ylim = [0.008, 0.022] # theta
    for i in range(3):
        axs[1,i].set_ylim(ylim)
        
    axs[1,0].add_patch(Rectangle((5, ylim[0]), 5, ylim[1], edgecolor=None, facecolor="gray", alpha=0.15))
    axs[1,0].set_xticks([7.5])
    axs[1,0].set_xticklabels(['S1'])
    axs[1,0].set_ylabel('Amplitude')
    # axs[1,0].yaxis.label.set_color('black')
    
    formatter = ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    formatter.set_powerlimits((-1,1)) 
    axs[1,0].yaxis.set_major_formatter(formatter)
    axs[1,0].get_yaxis().get_offset_text().set(size=6, position=(-0.5, 0.4))
    # axs[1,0].get_yaxis().get_offset_text().set_position((0, 0.5))
    axs[1,0].yaxis.offsetText.set_visible(False) # 隐藏 offset
    
    # axs[1,0].ticklabel_format(style='sci', axis='y')
    # axs[1,0].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    
    axs[1,1].get_yaxis().set_visible(False)
    axs[1,1].spines['left'].set_visible(False)
    axs[1,1].add_patch(Rectangle((5, ylim[0]), 5, ylim[1], edgecolor=None, facecolor="gray", alpha=0.15))
    axs[1,1].set_xticks([7.5])
    axs[1,1].set_xticklabels(['S2'])
    
    # axs[1,2].get_xaxis().set_visible(False)
    axs[1,2].get_yaxis().set_visible(False)
    axs[1,2].spines['left'].set_visible(False)
    axs[1,2].set_xticks([3.75])
    axs[1,2].set_xticklabels(['Delay'])
    
    figpath = r'C:\Wen\Project\Replay\Result_ocean\plot_example\figure'
    # fig.savefig(figpath + '\\Theta_S2_' + file + '.pdf', dpi='figure', bbox_inches='tight')


#%% plot fig1 TF example
file = r'C:\Wen\Project\Osc_control\Data\Result_ocean\plot_example(neuron_LFP)_data\Ocean_Fig1_TF_example(0713ch32).mat'
data = sio.loadmat(file)

file = r'C:\Wen\Project\Osc_control\Data\Result_ocean\plot_example(neuron_LFP)_data\Ocean_Fig1_TF_example_NoERPremoval.mat'
data2 = sio.loadmat(file) # ERP not removed

file = r'C:\Wen\Project\Osc_control\Data\Result_ocean\plot_example(neuron_LFP)_data\Ocean0713_Ch32_FreqVec.mat'
yaxs = sio.loadmat(file)['y_vec']

# y = data['img_cropped'] # ERP removed
y = data2['test_trim']  # ERP not removed
with plt.style.context(style_path):
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"][2:]
    fig, axs = plt.subplots(1, 1, figsize=(2.5, 1.5), dpi=800) # , constrained_layout=True
    plt.subplots_adjust(wspace=0.05, hspace=0.2)

    vmin, vmax = -1.8, 1.8
    cmap = plt.cm.jet
    f = axs.pcolormesh(y, cmap = cmap, vmin=vmin, vmax=vmax) # 
    axs.get_xaxis().set_visible(False)
    # axs.get_yaxis().set_visible(False)
    
    axs.vlines([data['S1on'], data['S2on']], 0, 200, colors='w')
    axs.set_yticks([30, 60, 107, 176], ['4', '10', '24', '120'])
    axs.set_ylabel('Frequency(Hz)')
    fig.colorbar(f, orientation='vertical')


    