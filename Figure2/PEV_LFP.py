# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 14:08:42 2022

@author: Wen
"""

get_ipython().run_line_magic('reset', '-sf')
style_path = r'C:\Wen\OneDrive\CodeHub\Python\style_paper.mplstyle'

import numpy as np
import pandas as pd
import os
from pathlib import Path

import matplotlib.pyplot as plt
# import h5py
import pickle as pkl
import scipy.io as sio
from scipy import stats

import statsmodels.api as sm
from statsmodels.formula.api import ols

from utils import fun_loadData
from utils import fun_pev
#%% functions
def remove_baseline(x):
    basemean = np.nanmean(x[:,:,:10], axis=2)
    for t in range(x.shape[2]):
        x[:,:,t] -= basemean 
    return x
#%%
seqL = 2
monkey = 'ocean'

if 'ocean' in monkey:
    # fpath = r'C:\Wen\Project\Osc_control\Data\Data_ocean\LFP\power13days'
    fpath = r'C:\Wen\Project\Osc_control\Data\Data_ocean\LFP\power_remove_erp'
    # filename = Path("C:\Wen\Project\Osc_control\Data\Data_ocean\LFP\power_remove_erp")
    alldays = ['0709','0710','0713','0714','0715','0720','0721','0722','0723','0726','0727','0729','0730'] #'0701','0702','0703','0705','0706','0707','0708',
    load_sel_chs = 0
    savepath = r'C:\Wen\Project\Osc_control\Data\Result_ocean\LFP\PEV_seqL2_frontal'
    
if 'grootsw' == monkey:
    sfreq = 1000
    seqL = 2
    fpath = r'C:\Wen\Project\Osc_control\Data\Data_grootsw\LFP_sw\power_remove_erp'
    alldays = ['04-15','05-07','05-09','05-13','05-14','05-15','05-16','05-20','05-21','05-22']
    load_sel_chs = 0
    savepath = r'C:\Wen\Project\Osc_control\Data\Result_grootsw\LFP_sw\PEV_seqL2_frontal'

#%% load LFP data
basecorr = False
select_rule = 0 # 1:repeat, 2:mirror, 0:all

region = ['PFC','PM','frontal'] #,'PM', 'PAR', 'M1']
osc_band = ['4_8', '80_120']
bandid = 0

aov_result_all = [None] * len(alldays)
for dayid in range(1,2): # len(alldays)
    days = [alldays[dayid]]
    # days = alldays
    fr_all = [None] * len(days)
    ch_info_all = [None] * len(days)
    for i in range(len(days)):
        if 'ocean' in monkey:
            file = fpath + '\\'+'frontal_tf200_'+osc_band[bandid]+'_morlet_removeMean_'+days[i]+'.pkl' # ocean
            # file = fpath + '\\'+'ch233_tf200_'+osc_band[bandid]+'_denoise_morlet_'+days[i]+'.pkl' # ocean
        if 'grootsw' in monkey:
            file = fpath + '\\'+'frontal_L2_tf200_'+osc_band[bandid]+'_morlet_'+days[i]+'.pkl' # 
            
        #### load LFP power neu_sit_info: ch*trial*time
        if load_sel_chs:
            neupath = r'C:\Wen\Project\Osc_control\Data\Data_lemon\SPK\FR_segment_win100step50_frontal'
            neufile = neupath + '\\2021-'+days[i]+'_fr_L2_win100step50_selChs.pkl'
            info, neu_sit_info, fr_arr, mkts_merge = fun_loadData.load_tf_data_selchs(file, region, select_rule, days[i], neufile) 
        else:
            info, neu_sit_info, fr_arr, mkts_merge = fun_loadData.load_tf_data(file, region, select_rule) 
            
        if len(days) > 1:
            # creat psudo trials; merge channels
            info, fr_arr = fun_loadData.creat_psudo_trial(info, fr_arr, seqL) 
            
        fr_all[i] = fr_arr # neuron*trial*time
        ch_info_all[i] = neu_sit_info
    
    if len(fr_all) > 1:
        idx = [idx for idx in range(len(fr_all)) if fr_all[idx].shape[1] == 30]
        fr_all = [fr_all[j] for j in idx]
        fr_all = np.concatenate(fr_all, axis=0)
        ch_info_all = [ch_info_all[j] for j in idx]
        ch_info_all = pd.concat(ch_info_all).reset_index(drop=True)
    else:    
        fr_all = fr_all[0]
        ch_info_all = ch_info_all[0]
        
    n_feature = fr_all.shape[0]
    n_sample = fr_all.shape[1]
    n_time = fr_all.shape[2]
    # %% PEV for each neuron
    x = fr_all.copy() # neuron*trial*time
    df = info.copy()
    
    pev, pval = fun_pev.pev(x, df, seqL)
    pkl.dump({'pev':pev, 'pval':pval,'ch_info_all':ch_info_all}, open(savepath + '\\removeERP20241014_pev_'+osc_band[bandid]+'_frontal_'+days[0]+'.pkl','wb'))
    #%% plot pev
    pev = pkl.load(open(savepath + '\\pev_'+osc_band[bandid]+'_frontal_'+days[0]+'.pkl', 'rb'))
    pval = pev['pval']
    pev = pev['pev']
    # mask = pval<0.01
    # pev = pev*mask
    title = osc_band[bandid] + ' pev (frontal)'
    sel_area = ch_info_all['area'].isin(region).values
    pev_region = pev[sel_area, 15:55, :]
    
    
    fun_pev.plot_pev(pev_region, seqL, mkts_merge, title)
    plt.show()
    # %% baseline permutation of each neuron
    x = fr_all.copy()
    df = info.copy()
    
    pevsf = fun_pev.pev_shuffle_lfp(x, df, seqL)
    pkl.dump({'pev':pevsf}, open(savepath + '\\removeERP20241014_pev_'+osc_band[bandid]+'_frontal_shuffle_'+days[0]+'.pkl','wb'))
    #%% permutation test 
    pev = pkl.load(open(savepath + '\\removeERP20241014_pev_'+osc_band[bandid]+'_frontal_'+days[0]+'.pkl','rb'))['pev']
    pevsf = pkl.load(open(savepath + '\\removeERP20241014_pev_'+osc_band[bandid]+'_frontal_shuffle_'+days[0]+'.pkl','rb'))['pev']
    pval = fun_pev.permu_test(pev, pevsf, seqL)    
    pkl.dump({'pval':pval, 'ch_info_all':ch_info_all}, open(savepath + '\\removeERP20241014_pev_'+osc_band[bandid]+'_frontal_permu_'+days[0]+'.pkl','wb'))
    #%% plot permutation mask
    do_plot=1
    if do_plot:
        pval = pkl.load(open(savepath + '\\pev_'+osc_band[bandid]+'_frontal_permu_'+days[0]+'.pkl', 'rb'))['pval']
        x = np.concatenate([pval[:,mkts_merge[1,0]:mkts_merge[1,2],:],pval[:,mkts_merge[0,3]-5:mkts_merge[0,3],:]], axis=1)
        
        # 画原始pev，盖上mask
        pev = np.concatenate([pev[:,15:45,:], pev[:,50:55,:]], axis=1)
        mask = x<0.05
        pev = pev * mask
        #####
        
        title = osc_band[bandid] + ' significant time (frontal)_' + days[0]
        figpath = r'C:\Wen\OneDrive\Project\Replay\Manuscript\Figures\Figure2'
        savename = [] #figpath + '\\' + title + '.pdf'
        # fun_pev.plot_pval(x, seqL, mkts_merge, title, savename)
        fun_pev.plot_pev_mask(pev, seqL, mkts_merge, title, savename) # 画原始pev，盖上mask
        plt.show()
    
