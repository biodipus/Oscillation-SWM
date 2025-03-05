# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 13:48:41 2022

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
from scipy import stats

import statsmodels.api as sm
from statsmodels.formula.api import ols

from utils import fun_loadData
from utils import fun_pev
# %%
monkey = 'ocean'
if 'ocean' in monkey:
    alldays = ['0709','0710','0713','0714','0715','0720','0721','0722','0723','0726','0727','0729','0730'] #'0701','0702','0703','0705','0706','0707','0708',
    # alldays = alldays[7:]
    fpath = r'C:\Wen\Project\Osc_control\Data\Data_ocean\SPK\FR_segment_win100step50'
    savepath = r'C:\Wen\Project\Osc_control\Data\Result_ocean\SPK\PEV_seqL2_frontal'
    seqL = 2
    # days= ['0710']

#%% load data
select_rule = 0
region = ['PFC','PM','frontal']#, 'PM', 'PAR', 'M1']
for dayid in range(1,2):
    days = [alldays[dayid]]
    print(days[0])
    fr_all = [None] * len(days)
    ch_info_all = [None] * len(days)
    for i in range(len(days)):
        file = fpath + '\\'+days[i]+'_fr_L2_win100step50_selChs.pkl'
        info, neu_sit_info, fr_arr, mkts_merge, spk_arr, mktime = fun_loadData.load_fr_data(file, region, select_rule)
        if len(days) > 1:
            info, fr_arr = fun_loadData.creat_psudo_trial(info, fr_arr) # creat psudo trials
        fr_all[i] = fr_arr # neuron*trial*time
        ch_info_all[i] = neu_sit_info
    
    if len(days) > 1:
        fr_all = np.concatenate(fr_all, axis=0)
        ch_info_all = pd.concat(ch_info_all).reset_index()
    else:    
        fr_all = fr_all[0]
        ch_info_all = ch_info_all[0]
    #%% PEV
    x = fr_all.copy() # neuron*trial*time
    df = info.copy()
    pev, pval = fun_pev.pev(x, df, seqL)
    
    pkl.dump({'pev':pev,'pval':pval,'ch_info_all':ch_info_all}, open(savepath + '\\pev_frontal_'+days[0]+'.pkl', 'wb'))
    #%% plot pev
    title = 'Spk pev (frontal) ' + days[0]
    sel_area = ch_info_all['area'].isin(region)
    pev_region = pev[sel_area, mkts_merge[1,0]:mkts_merge[1,3], :]
    
    fun_pev.plot_pev(pev_region, seqL, mkts_merge, title)
    plt.show()
    # %% baseline permutation of each neuron
    x = fr_all.copy()
    df = info.copy()
    pevsf = fun_pev.pev_shuffle_spk(x, df, seqL)
    
    pkl.dump({'pev': pevsf}, open(savepath + '\\pev_frontal_shuffle_'+days[0]+'.pkl', 'wb'))
    # %% permutation test of PEV
    pev = pkl.load(open(savepath + '\\pev_frontal_'+days[0]+'.pkl', 'rb'))['pev']
    pevsf = pkl.load(open(savepath + '\\pev_frontal_shuffle_'+days[0]+'.pkl', 'rb'))['pev']
    
    pval = fun_pev.permu_test(pev, pevsf, seqL)
    
    pkl.dump({'pval': pval}, open(savepath + '\\pev_frontal_permu_'+ days[0]+'.pkl', 'wb'))
    #%% plot permutation mask
    pval = pkl.load(open(savepath + '\\pev_frontal_permu_'+days[0]+'.pkl', 'rb'))['pval']
    x = np.concatenate([pval[:,mkts_merge[1,0]:mkts_merge[1,2],:],pval[:,mkts_merge[0,3]-5:mkts_merge[0,3],:]], axis=1)
    
    # 画原始pev，盖上mask
    pev = np.concatenate([pev[:,15:45,:], pev[:,50:55,:]], axis=1)
    mask = x<0.05
    pev = pev * mask  
    #####
    
    title = 'SPK significant time (frontal)_' + days[0]
    figpath = r'C:\Wen\OneDrive\Project\Replay\Manuscript\Figures\Figure2'
    savename = [] #figpath + '\\' + title + '.pdf'
    # fun_pev.plot_pval(x, seqL, mkts_merge, title, savename)
    fun_pev.plot_pev_mask(pev, seqL, mkts_merge, title, savename)# 画原始pev，盖上mask
    plt.show()
    
   