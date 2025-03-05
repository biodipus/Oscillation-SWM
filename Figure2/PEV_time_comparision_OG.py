# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 16:32:04 2023

@author: Wen
"""

get_ipython().run_line_magic('reset', '-sf')
style_path = r'C:\Wen\OneDrive\CodeHub\Python\style_paper.mplstyle'

import numpy as np
import pandas as pd
import os
# import h5py
import pickle as pkl
import scipy.io as sio


import matplotlib.pyplot as plt
import seaborn as sns

from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols

from utils import fun_loadData
from utils import fun_pev

#%% load ocean
days = ['0709','0710','0713','0714','0715','0720','0721','0722','0723','0726','0727','0729','0730'] #'0701','0702','0703','0705','0706','0707','0708',
# spkpath = r'C:\Wen\Project\Osc_control\Data\Result_ocean\SPK\PEV_seqL2'
# lfppath = r'C:\Wen\Project\Osc_control\Data\Result_ocean\LFP\PEV_seqL2'

pevpath = r'C:\Wen\Project\Osc_control\Data\Result_ocean\LFP\PEV_seqL2_frontal'

on_time_all = [None] * len(days)
for i in range(len(days)):
    pval_theta = pkl.load(open(pevpath + '\\removeERP20241014_pev_4_8_frontal_permu_'+days[i]+'.pkl', 'rb'))['pval']
    pval_gamma = pkl.load(open(pevpath + '\\pev_80_120_frontal_permu_'+days[i]+'.pkl', 'rb'))['pval']
    
    pval_theta = pval_theta[:,15:55,:]
    pval_gamma = pval_gamma[:,15:55,:]
    ####
    # sit_info = neu_sit_info.copy()
    # pval_theta = pval_theta[sit_info['chid in pfc'],:,:]
    # pval_gamma = pval_gamma[sit_info['chid in pfc'],:,:]

    time_distribution = np.zeros([5000, 4])+np.nan
    
    aa = pval_theta[:,:,0]<0.05
    aa = (np.where(aa)[1]-3)*0.05
    time_distribution[:len(aa), 0] = aa
    
    aa = pval_theta[:,:,1]<0.05
    aa = (np.where(aa)[1]-19)*0.05
    time_distribution[:len(aa), 1] = aa
    
    aa = pval_gamma[:,:,0]<0.05
    aa = (np.where(aa)[1]-4)*0.05
    time_distribution[:len(aa), 2] = aa
    
    aa = pval_gamma[:,:,1]<0.05
    aa = (np.where(aa)[1]-19)*0.05
    time_distribution[:len(aa), 3] = aa
    
    on_time_all[i] = time_distribution
#%% load grootsw lev
pevpath = r'C:\Wen\Project\Osc_control\Data\Result_grootsw\LFP_sw\PEV_seqL2_frontal'
fpath = r'C:\Wen\Project\Osc_control\Data\Data_grootsw\LFP_sw\power_remove_erp'
days = ['04-15','05-07','05-09','05-13','05-14','05-15','05-16','05-20','05-21','05-22']

# select_rule = 0 # 1:repeat, 2:mirror, 0:all
# region = ['PFC','PM','frontal'] #,'PM', 'PAR', 'M1']
# file = fpath + '\\'+'frontal_L2_tf200_4_8_morlet_'+days[i]+'.pkl' # 
# info, neu_sit_info, fr_arr, mkts_merge = fun_loadData.load_tf_data(file, region, select_rule) 
    
on_time_all = [None] * len(days)
for i in range(len(days)):
    pval_theta = pkl.load(open(pevpath + '\\removeERP20241014_pev_4_8_frontal_permu_'+days[i]+'.pkl', 'rb'))['pval']
    pval_gamma = pkl.load(open(pevpath + '\\removeERP20241014_pev_80_120_frontal_permu_'+days[i]+'.pkl', 'rb'))['pval']
    pval_theta = pval_theta[:,15:55,:]
    pval_gamma = pval_gamma[:,15:55,:]
    ####
    # sit_info = neu_sit_info.copy()
    # pval_theta = pval_theta[sit_info['chid in pfc'],:,:]
    # pval_gamma = pval_gamma[sit_info['chid in pfc'],:,:]

    time_distribution = np.zeros([5000, 4])+np.nan
    
    aa = pval_theta[:,:,0]<0.05
    aa = (np.where(aa)[1]-3)*0.05
    time_distribution[:len(aa), 0] = aa
    
    aa = pval_theta[:,:,1]<0.05
    aa = (np.where(aa)[1]-19)*0.05
    time_distribution[:len(aa), 1] = aa
    
    aa = pval_gamma[:,:,0]<0.05
    aa = (np.where(aa)[1]-4)*0.05
    time_distribution[:len(aa), 2] = aa
    
    aa = pval_gamma[:,:,1]<0.05
    aa = (np.where(aa)[1]-19)*0.05
    time_distribution[:len(aa), 3] = aa
    
    on_time_all[i] = time_distribution

#%% plot sig time distribution
from scipy.stats import ks_2samp

x = np.concatenate(on_time_all, axis=0)

# 进行双样本K-S检验
x1 = x[:,0]
x1 = x1[~np.isnan(x1)]
x2 = x[:,2]
x2 = x2[~np.isnan(x2)]
ks_2samp(x1, x2)

x1 = x[:,1]
x1 = x1[~np.isnan(x1)]
x2 = x[:,3]
x2 = x2[~np.isnan(x2)]
ks_2samp(x1, x2)

with plt.style.context(style_path):
    fig, axs = plt.subplots(1, 2, figsize=(3, 1.7), dpi=300, constrained_layout=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"][3:5][::-1]
    # axs[0].hist(x[:,[0,2]], bins=10, range=[0.0, 2.0], alpha=0.8)
    # axs[1].hist(x[:,[1,3]], bins=10, range=[0.0, 1.4], alpha=0.8)
    ax = sns.kdeplot(data=x[:,[0,2]], common_norm=False,fill=True,alpha=.5,linewidth=0.6,ax=axs[0], palette=colors)
    axs[0].set_xlabel('Time from S1 onset')
    
    sns.kdeplot(data=x[:,[1,3]], legend=False, common_norm=False,fill=True,alpha=.5,linewidth=0.6,ax=axs[1], palette=colors)
    axs[1].set_xlabel('Time from S2 onset')
    
    ax.legend(labels=['HG', 'Theta'], loc='upper right')

kurtosis_data = np.zeros([len(days), 4])
for i in range(len(days)):
    x = on_time_all[i]
    for j in range(4):
        # kurtosis_data[i,j] = kurtosis(x[:,j], fisher=False, nan_policy='omit')  # fisher=False表示计算常规峰度
        kurtosis_data[i,j] = np.nanstd(x[:,j])  # fisher=False表示计算常规峰度
stats.ttest_rel(kurtosis_data[:,0], kurtosis_data[:,2])
stats.ttest_rel(kurtosis_data[:,1], kurtosis_data[:,3])
#%%
xx = x[:,0]
xx = xx[~np.isnan(xx)]
# 计算均值和标准误
mean = np.mean(xx)
sem = stats.sem(xx)  # 标准误

# 计算95%置信区间
confidence = 0.95
h = sem * stats.t.ppf((1 + confidence) / 2, len(xx) - 1)

# 区间范围
ci_lower = mean - h
ci_upper = mean + h
#%%
with plt.style.context(style_path):
    fig, axs = plt.subplots(1, 2, figsize=(3, 1.7), dpi=300, constrained_layout=True)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"][3:5][::-1]
    
    y1 = kurtosis_data
    for i in range(y1.shape[0]):
        axs[0].plot([1.1,1.9], y1[i,[0,2]], linestyle='-', color='gray', alpha=0.7, zorder=0)
        axs[0].scatter([1], y1[i,[0]], color=colors[0], edgecolor=[], alpha=0.8, zorder=1)
        axs[0].scatter([2], y1[i,[2]], color=colors[1], edgecolor=[], alpha=0.8, zorder=1)
        
        axs[1].plot([1.1,1.9], y1[i,[1,3]], linestyle='-', color='gray', alpha=0.7, zorder=0)
        axs[1].scatter([1], y1[i,[1]], color=colors[0], edgecolor=[], alpha=0.8, zorder=1)
        axs[1].scatter([2], y1[i,[3]], color=colors[1], edgecolor=[], alpha=0.8, zorder=1)
        
    axs[0].set_xlim([0.5, 2.5])
    axs[0].set_ylim([0.1, 0.65])
    axs[0].set_xticks([1,2], ['Theta', 'HG'])
    axs[0].set_ylabel('Standard deviation')
    
    axs[1].set_xlim([0.5, 2.5])
    axs[1].set_ylim([0.0, 0.35])
    axs[1].set_xticks([1,2], ['Theta', 'HG'])
    axs[1].set_ylabel('Standard deviation')




        
