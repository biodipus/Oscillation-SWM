# -*- coding: utf-8 -*-
"""
Created on Fri Dec 29 15:17:35 2023

@author: Wen
"""
get_ipython().run_line_magic('reset', '-sf')
style_path = r'C:\Wen\OneDrive\CodeHub\Python\style_paper.mplstyle'

import sys 
import numpy as np
import pandas as pd
from sklearn import preprocessing
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import random
from scipy.sparse.linalg import eigsh as ssl_eigsh
from scipy import linalg
import scipy.io as sio
import pickle as pkl
import glob

import seaborn as sns

#%% spk memory ~ acc1+acc2
path = r'C:\Wen\Project\Osc_control\Data\monkey_combined\SPK_fw_SFC_corr\2-Seq_perm_stats_Switch_SFCxSPK-FW'

file = '\\GrootSwitch10Days_ANOVA-AllSig_SPK_Mem_FW_x_SFCacc_perm_stats_2Seq.mat'
# file = '\\Ocean13Days_ANOVA-AllSig_SPK_Mem_FW_x_SFCacc_perm_stats_2Seq.mat'


# Beta
# path = r'C:\Wen\Project\Osc_control\Data\monkey_combined\Beta_SFC_FW_corr' # Beta
# file = '\\Ocean13Days_ANOVA-AllSig_SPK_Mem_FW_x_BT-SFC_perm_stats_2Seq.mat' 
# file = '\\GrootSwitch10Days_ANOVA-AllSig_SPK_Mem_FW_x_BT-SFC_perm_stats_2Seq.mat' 


file = path + file
data = sio.loadmat(file)

with plt.style.context(style_path):
    fig = plt.figure(figsize = (3.8, 1.7), dpi = 300)
    ax = fig.add_subplot(1, 2, 1)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    
    df = pd.DataFrame({'Acc R1b1':data['beta_hold_S1'][0,:], 'Acc R1b2':data['beta_hold_S1'][1,:],
                       'Acc R2b1':data['beta_hold_S2'][0,:], 'Acc R2b2':data['beta_hold_S2'][1,:]})
    
    #### plot rank-1
    sns.violinplot(data=df.loc[:,['Acc R1b1','Acc R1b2']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.plot([0,1], [data['S1_R1_acc_beta'][0], data['S1_R2_acc_beta'][0]], color='gray')
    ax.scatter([0,1], [data['S1_R1_acc_beta'], data['S1_R2_acc_beta']], color = colors[0],
               marker='+', s=60, zorder=3)
    
    ax.set_xlim([-0.6, 1.6])
    ax.set_ylim([-0.5, 0.75])
    # ax.set_ylim([-2.3, 2.3])
    ax.set_xticklabels(['S1', 'S2'], rotation=0)
    ax.set_xlabel('SFC strength')
    ax.set_ylabel('Coefficient')
    # ax.set_title('Rank-1')
    
    #### plot rank-2
    ax = fig.add_subplot(1, 2, 2)
    sns.violinplot(data=df.loc[:,['Acc R2b1','Acc R2b2']], bw_adjust=.5, cut=1, linewidth=0, palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.plot([0,1], [data['S2_R1_acc_beta'][0], data['S2_R2_acc_beta'][0]], color='gray')
    ax.scatter([0,1], [data['S2_R1_acc_beta'], data['S2_R2_acc_beta']], color = colors[1],
               marker='+', s=60, zorder=3)
    
    ax.set_xlim([-0.6, 1.6])
    ax.set_ylim([-0.5, 0.75])
    ax.set_xticklabels(['S1', 'S2'], rotation=0)
    ax.set_xlabel('SFC strength')
    # ax.set_yticklabels([])
    # ax.set_title('Rank-2')
    
    plt.tight_layout()
#%% spk memory ~ acc1+acc2 firing rate
path = r'C:\Wen\Project\Osc_control\Data\monkey_combined\SPK_fw_SFC_corr\2-Seq_perm_Switch_FiringRateControls\SFC_x_FW_stats'
# file = '\\GrootSwitch10Days_ANOVA-AllSig_SPK_Mem_FW_x_FRctrl_perm_stats_2Seq'
file = '\\Ocean13Days_ANOVA-AllSig_SPK_Mem_FW_x_FRctrl_perm_stats_2Seq'

file = path + file
data = sio.loadmat(file)

with plt.style.context(style_path):
    fig = plt.figure(figsize = (3.8, 2), dpi = 300)
    ax = fig.add_subplot(1, 2, 1)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    
    df = pd.DataFrame({'Acc R1b1':data['beta_hold_S1'][0,:], 'Acc R1b2':data['beta_hold_S1'][1,:],
                       'Acc R2b1':data['beta_hold_S2'][0,:], 'Acc R2b2':data['beta_hold_S2'][1,:]})
    
    #### plot rank-1
    sns.violinplot(data=df.loc[:,['Acc R1b1','Acc R1b2']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.plot([0,1], [data['S1_R1_acc_beta'][0], data['S1_R2_acc_beta'][0]], color='gray')
    ax.scatter([0,1], [data['S1_R1_acc_beta'], data['S1_R2_acc_beta']], color = colors[0],
               marker='+', s=60, zorder=3)
    
    ax.set_xlim([-0.6, 1.6])
    ax.set_ylim([-1.5, 1.5])
    # ax.set_ylim([-2.3, 2.3])
    ax.set_xticklabels(['FR (S1)', 'FR (S2)'], rotation=30)
    ax.set_ylabel('Coefficient')
    # ax.set_title('Rank-1')
    
    #### plot rank-2
    ax = fig.add_subplot(1, 2, 2)
    sns.violinplot(data=df.loc[:,['Acc R2b1','Acc R2b2']], bw_adjust=.5, cut=1, linewidth=0, palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.plot([0,1], [data['S2_R1_acc_beta'][0], data['S2_R2_acc_beta'][0]], color='gray')
    ax.scatter([0,1], [data['S2_R1_acc_beta'], data['S2_R2_acc_beta']], color = colors[1],
               marker='+', s=60, zorder=3)
    
    ax.set_xlim([-0.6, 1.6])
    ax.set_ylim([-1.5, 1.5])
    # ax.set_ylim([-2.3, 2.3])
    ax.set_xticklabels(['FR (S1)', 'FR (S2)'], rotation=30)
    # ax.set_yticklabels([])
    # ax.set_title('Rank-2')
    
    plt.tight_layout()
#%% spk memory ~ acc1+acc2, neuron with equal FR S1S2
path = r'C:\Wen\Project\Osc_control\Data\monkey_combined\SPK_fw_SFC_corr\neuron_with_equalFR_S1S2'
file = '\\Ocean13Days_S1vsS2_FR_NoDiff_SPK_Mem_FW_x_SFCacc_perm_stats_2Seq'
# file = '\\GrootSwitch10Days_S1vsS2_FR_NoDiff_SPK_Mem_FW_x_SFCacc_perm_stats_2Seq'

file = path + file
data = sio.loadmat(file)

with plt.style.context(style_path):
    fig = plt.figure(figsize = (3, 1.7), dpi = 300)
    ax = fig.add_subplot(1, 2, 1)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    
    df = pd.DataFrame({'Acc R1b1':data['beta_hold_S1'][0,:], 'Acc R1b2':data['beta_hold_S1'][1,:],
                       'Acc R2b1':data['beta_hold_S2'][0,:], 'Acc R2b2':data['beta_hold_S2'][1,:]})
    
    #### plot rank-1
    sns.violinplot(data=df.loc[:,['Acc R1b1','Acc R1b2']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.plot([0,1], [data['S1_R1_acc_beta'][0], data['S1_R2_acc_beta'][0]], color='gray')
    ax.scatter([0,1], [data['S1_R1_acc_beta'], data['S1_R2_acc_beta']], color = colors[0],
               marker='+', s=60, zorder=3)
    
    ax.set_xlim([-0.6, 1.6])
    ax.set_ylim([-1., 1.])
    # ax.set_ylim([-2.3, 2.3])
    ax.set_xticklabels(['S1', 'S2'], rotation=0)
    ax.set_xlabel('SFC strength')
    ax.set_ylabel('Coefficient')
    # ax.set_title('Rank-1')
    
    #### plot rank-2
    ax = fig.add_subplot(1, 2, 2)
    sns.violinplot(data=df.loc[:,['Acc R2b1','Acc R2b2']], bw_adjust=.5, cut=1, linewidth=0, palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.plot([0,1], [data['S2_R1_acc_beta'][0], data['S2_R2_acc_beta'][0]], color='gray')
    ax.scatter([0,1], [data['S2_R1_acc_beta'], data['S2_R2_acc_beta']], color = colors[1],
               marker='+', s=60, zorder=3)
    
    ax.set_xlim([-0.6, 1.6])
    ax.set_ylim([-1., 1.])
    # ax.set_ylim([-2.3, 2.3])
    ax.set_xticklabels(['S1', 'S2'], rotation=0)
    ax.set_xlabel('SFC strength')
    
    plt.tight_layout()
#%% spk memory ~ acc1+acc2 theta power
path = r'C:\Wen\Project\Osc_control\Data\monkey_combined\SPK_fw_SFC_corr\2-Seq_perm_Switch_FiringRateControls\SFC_x_FW_stats'
# file = '\\Ocean13Days_ANOVA-AllSig_SPK_Mem_FW_x_TH-ERSP_perm_stats_2Seq'
file = '\\GrootSwitch10Days_ANOVA-AllSig_SPK_Mem_FW_x_TH-ERSP_perm_stats_2Seq'

file = path + file
data = sio.loadmat(file)

with plt.style.context(style_path):
    fig = plt.figure(figsize = (3, 2), dpi = 300)
    ax = fig.add_subplot(1, 2, 1)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    
    df = pd.DataFrame({'Acc R1b1':data['beta_hold_S1'][0,:], 'Acc R1b2':data['beta_hold_S1'][1,:],
                       'Acc R2b1':data['beta_hold_S2'][0,:], 'Acc R2b2':data['beta_hold_S2'][1,:]})
    
    #### plot rank-1
    sns.violinplot(data=df.loc[:,['Acc R1b1','Acc R1b2']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.plot([0,1], [data['S1_R1_acc_beta'][0], data['S1_R2_acc_beta'][0]], color='gray')
    ax.scatter([0,1], [data['S1_R1_acc_beta'], data['S1_R2_acc_beta']], color = colors[0],
               marker='+', s=60, zorder=3)
    
    ax.set_xlim([-0.6, 1.6])
    ax.set_ylim([-1.5, 1.5])
    # ax.set_ylim([-2.3, 2.3])
    ax.set_xticklabels(['Power (S1)', 'Power (S2)'], rotation=30)
    ax.set_ylabel('Coefficient')
    ax.set_title('Rank-1')
    
    #### plot rank-2
    ax = fig.add_subplot(1, 2, 2)
    sns.violinplot(data=df.loc[:,['Acc R2b1','Acc R2b2']], bw_adjust=.5, cut=1, linewidth=0, palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.plot([0,1], [data['S2_R1_acc_beta'][0], data['S2_R2_acc_beta'][0]], color='gray')
    ax.scatter([0,1], [data['S2_R1_acc_beta'], data['S2_R2_acc_beta']], color = colors[1],
               marker='+', s=60, zorder=3)
    
    ax.set_xlim([-0.6, 1.6])
    ax.set_ylim([-1.5, 1.5])
    # ax.set_ylim([-2.3, 2.3])
    ax.set_xticklabels(['FR (S1)', 'FR (S2)'], rotation=30)
    ax.set_yticklabels([])
    ax.set_title('Rank-2')
    
    plt.tight_layout()
#%% spk memory ~ mod1 + mod2
path = r'C:\Wen\Project\Osc_control\Data\Result_ocean\SFC_fw_corr'
file = '\\Ocean13Days_SPK_Mem_FW_x_SFCacc_perm_stats_2Seq_ModulatorCtrl'

file = path + file
data = sio.loadmat(file)

with plt.style.context(style_path):
    fig = plt.figure(figsize = (2, 1.6), dpi = 300)
    ax = fig.add_subplot(1, 2, 1)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    
    df = pd.DataFrame({'Acc R1b1':data['beta_hold_S1'][0,:], 'Acc R1b2':data['beta_hold_S1'][1,:],
                       'Acc R2b1':data['beta_hold_S2'][0,:], 'Acc R2b2':data['beta_hold_S2'][1,:]})
    
    #### plot rank-1
    sns.violinplot(data=df.loc[:,['Acc R1b1','Acc R1b2']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.plot([0,1], [data['S1_R1_acc_beta'][0], data['S1_R2_acc_beta'][0]], color='gray')
    ax.scatter([0,1], [data['S1_R1_acc_beta'], data['S1_R2_acc_beta']], color = colors[0],
               marker='_', s=60, zorder=3)
    
    ax.set_xlim([-0.6, 1.6])
    # ax.set_ylim([-2.9, 2.9])
    ax.set_ylim([-0.3, 0.3])
    ax.set_xticklabels(['Mod(R1)', 'Mod(R2)'], rotation=30)
    ax.set_ylabel('Coefficient')
    ax.set_title('Rank-1')
    
    #### plot rank-2
    ax = fig.add_subplot(1, 2, 2)
    sns.violinplot(data=df.loc[:,['SFC R2b1','SFC R2b2']], bw_adjust=.5, cut=1, linewidth=0, palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.plot([0,1], [data['S2_R1_acc_beta'][0], data['S2_R2_acc_beta'][0]], color='gray')
    ax.scatter([0,1], [data['S2_R1_acc_beta'], data['S2_R2_acc_beta']], color = colors[1],
               marker='_', s=60, zorder=3)
    
    ax.set_xlim([-0.6, 1.6])
    # ax.set_ylim([-2.9, 2.9])
    ax.set_ylim([-0.3, 0.3])
    ax.set_xticklabels(['Mod(R1)', 'Mod(R2)'], rotation=30)
    ax.set_yticklabels([])
    ax.set_title('Rank-2')
    
    plt.tight_layout()
#%% spk memory ~ acc1+acc2+acc3, L3
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.formula.api import ols

path = r'C:\Wen\Project\Osc_control\Data\Result_ocean\SFC_fw_corr'

file = '\\Ocean3Seq_FW_concat_NewFW2024_srm_ZC_Frontal_321_adj.mat'
fw = sio.loadmat(path+file)

file = '\\Ocean3Seq_SFCvec_concat_BaseComp_NoSubtract_AddEnt_TrueIdxVerify_Frontal_NoMask.mat'
sfc = sio.loadmat(path + file)

file = '\\Ocean3Seq_SPK_ANOVA_concat_Frontal'
anova = sio.loadmat(path + file)
idx = ~np.isnan(anova['ANOVA_vec_R3'])

df = pd.DataFrame({'fw1': fw['SPK_S1_D_Comb_ansc'][:,0], 'fw2': fw['SPK_S2_D_Comb_ansc'][:,0], 'fw3': fw['SPK_S3_D_Comb_ansc'][:,0],
                   'sfc1': sfc['SFC_T1c_CombZ'][:,0],'sfc2': sfc['SFC_T2c_CombZ'][:,0],'sfc3': sfc['SFC_T3c_CombZ'][:,0]})
df = df.loc[idx, :]

formula = "fw3 ~  sfc1 + sfc2 + sfc3"
mod = smf.ols(formula, data=df)
res = mod.fit()
print(res.summary())

plt.scatter(df['sfc1'], df['sfc2'])
plt.scatter(df['sfc2'], df['sfc3'])
#%% spk memory, L3
path = r'C:\Wen\Project\Osc_control\Data\monkey_combined\SPK_fw_SFC_corr\SPK_Mem_FW_x_SFCacc_perm_stats_3Seq'
file = '\\Ocean_SPK_Mem_FW_x_SFCacc_perm_stats_3Seq_3XLM'

file = path + file
data = sio.loadmat(file)

with plt.style.context(style_path):
    fig = plt.figure(figsize = (3.6, 1.6), dpi = 300)
    ax = fig.add_subplot(1, 3, 1)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    
    df = pd.DataFrame({'Acc R1b1':data['beta_hold_S1'][0,:], 'Acc R1b2':data['beta_hold_S1'][1,:],'Acc R1b3':data['beta_hold_S1'][2,:],
                       'Acc R2b1':data['beta_hold_S2'][0,:], 'Acc R2b2':data['beta_hold_S2'][1,:],'Acc R2b3':data['beta_hold_S2'][2,:],
                       'Acc R3b1':data['beta_hold_S3'][0,:], 'Acc R3b2':data['beta_hold_S3'][1,:],'Acc R3b3':data['beta_hold_S3'][2,:]})
    
    #### plot rank-1
    sns.violinplot(data=df.loc[:,['Acc R1b1','Acc R1b2','Acc R1b3']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.plot([0,1,2], [data['S1_R1_acc_beta'][0], data['S1_R2_acc_beta'][0], data['S1_R3_acc_beta'][0]], color='gray')
    ax.scatter([0,1,2], [data['S1_R1_acc_beta'], data['S1_R2_acc_beta'], data['S1_R3_acc_beta']], color = colors[0],
               marker='+', s=60, zorder=3)
    
    # ax.set_xlim([-0.3, 1.6])
    # ax.set_ylim([-0.01, 0.015])
    ax.set_ylim([-0.5, 0.5])
    ax.set_xticklabels(['SFC(S1)', 'SFC(S2)', 'SFC(S3)'], rotation=30)
    ax.set_ylabel('Coefficient')
    ax.set_title('Memory S1', color=colors[0])
    
    #### plot rank-2
    ax = fig.add_subplot(1, 3, 2)
    sns.violinplot(data=df.loc[:,['Acc R2b1','Acc R2b2','Acc R2b3']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.plot([0,1,2], [data['S2_R1_acc_beta'][0], data['S2_R2_acc_beta'][0], data['S2_R3_acc_beta'][0]], color='gray')
    ax.scatter([0,1,2], [data['S2_R1_acc_beta'], data['S2_R2_acc_beta'], data['S2_R3_acc_beta']], color = colors[1],
               marker='+', s=60, zorder=3)
    
    # ax.set_xlim([-0.3, 1.6])
    # ax.set_ylim([-0.01, 0.015])
    ax.set_ylim([-0.5, 0.5])
    ax.set_xticklabels(['SFC(S1)', 'SFC(S2)', 'SFC(S3)'], rotation=30)
    ax.set_yticklabels([])
    ax.set_ylabel(' ')
    ax.set_title('Memory S2', color=colors[1])
    
    #### plot rank-3
    ax = fig.add_subplot(1, 3, 3)
    sns.violinplot(data=df.loc[:,['Acc R3b1','Acc R3b2','Acc R3b3']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.plot([0,1,2], [data['S3_R1_acc_beta'][0], data['S3_R2_acc_beta'][0], data['S3_R3_acc_beta'][0]], color='gray')
    ax.scatter([0,1,2], [data['S3_R1_acc_beta'], data['S3_R2_acc_beta'], data['S3_R3_acc_beta']], color = colors[2],
               marker='+', s=60, zorder=3)
    
    # ax.set_xlim([-0.3, 1.6])
    # ax.set_ylim([-0.01, 0.015])
    ax.set_ylim([-0.5, 0.5])
    ax.set_xticklabels(['SFC(S1)', 'SFC(S2)', 'SFC(S3)'], rotation=30)
    ax.set_yticklabels([])
    ax.set_ylabel(' ')
    ax.set_title('Memory S3', color=colors[2])
    
    # plt.tight_layout()

#%% spk memory, L3, pair
path = r'C:\Wen\Project\Osc_control\Data\monkey_combined\SPK_fw_SFC_corr\SPK_Mem_FW_x_SFCacc_perm_stats_3Seq'
file = '\\Ocean_SPK_Mem_FW_x_SFCacc_perm_stats_3Seq_2XLM'

file = path + file
data = sio.loadmat(file)

with plt.style.context(style_path):
    fig = plt.figure(figsize = (3.6, 1.6), dpi = 300)
    ax = fig.add_subplot(1, 3, 1)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    
    df = pd.DataFrame({'Acc R2b1':data['beta_hold_1'][0,:], 'Acc R2b2':data['beta_hold_1'][1,:],
                       'Acc R3b1':data['beta_hold_3'][0,:], 'Acc R3b3':data['beta_hold_3'][1,:],
                       'Acc R3b2':data['beta_hold_2'][0,:], 'Acc R3b3_2':data['beta_hold_2'][1,:]})
    
    #### plot rank-2
    sns.violinplot(data=df.loc[:,['Acc R2b1','Acc R2b2']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.plot([0,1], [data['S2_R1_acc_beta'][0], data['S2_R2_acc_beta'][0]], color='gray')
    ax.scatter([0,1], [data['S2_R1_acc_beta'], data['S2_R2_acc_beta']], color = colors[1],
               marker='+', s=60, zorder=3)
    
    ax.set_xlim([-0.5, 1.8])
    ax.set_ylim([-.5, .5])
    ax.set_xticklabels(['Acc(R1)', 'Acc(R2)'], rotation=30)
    ax.set_ylabel('Coefficient')
    ax.set_title('Rank-2 (1-2 pair)')
    
    #### plot rank-3, acc1-3
    ax = fig.add_subplot(1, 3, 2)
    sns.violinplot(data=df.loc[:,['Acc R3b1','Acc R3b3']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.plot([0,1], [data['S3_R1_acc_beta'][0], data['S3_R3_acc_beta2'][0]], color='gray')
    ax.scatter([0,1], [data['S3_R1_acc_beta'], data['S3_R3_acc_beta2']], color = colors[2],
               marker='+', s=60, zorder=3)
    
    ax.set_xlim([-0.5, 1.8])
    ax.set_ylim([-.5, .5])
    ax.set_xticklabels(['Acc(R1)', 'Acc(R3)'], rotation=30)
    ax.set_yticklabels([])
    ax.set_ylabel(' ')
    ax.set_title('Rank-3 (1-3 pair)')
    
    #### plot rank-3, acc2-3
    ax = fig.add_subplot(1, 3, 3)
    sns.violinplot(data=df.loc[:,['Acc R3b2','Acc R3b3_2']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.plot([0,1], [data['S3_R2_acc_beta'][0], data['S3_R3_acc_beta1'][0]], color='gray')
    ax.scatter([0,1], [data['S3_R2_acc_beta'], data['S3_R3_acc_beta1']], color = colors[2],
               marker='+', s=60, zorder=3)
    
    ax.set_xlim([-0.5, 1.8])
    ax.set_ylim([-.5, .5])
    ax.set_xticklabels(['Acc(R2)', 'Acc(R3)'], rotation=30)
    ax.set_yticklabels([])
    ax.set_ylabel(' ')
    ax.set_title('Rank-3 (2-3 pair)')
    
    # plt.tight_layout()
#%% spk memory, stepwise
path = r'C:\Wen\Project\Osc_control\Data\monkey_combined\FW_SFC_stepwise'
file = '\\Groot_SPK_Mem_FW_x_SFCacc_perm_stats_3Seq_StepLM'

file = path + file
data = sio.loadmat(file)

with plt.style.context(style_path):
    fig = plt.figure(figsize = (3.6, 1.6), dpi = 300)
    ax = fig.add_subplot(1, 3, 1)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    
    df = pd.DataFrame({'Acc R1b1':data['beta_hold_S1'][0,:], 'Acc R1b2':0,'Acc R1b3':data['beta_hold_S1'][1,:],
                       'Acc R2b1':data['beta_hold_S2'][0,:], 'Acc R2b2':data['beta_hold_S2'][1,:],'Acc R2b3':0,
                       'Acc R3b1':0, 'Acc R3b2':0,'Acc R3b3':data['beta_hold_S3'][0,:]})
    
    #### plot rank-1
    sns.violinplot(data=df.loc[:,['Acc R1b1','Acc R1b2','Acc R1b3']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.scatter([0,2], [data['S1_R1_acc_beta'], data['S1_R3_acc_beta']], color = colors[0],
               marker='+', s=60, zorder=3)
    # ax.scatter([1], [0], color = colors[0], linestyle='--',
    #            marker='+', s=60, zorder=3)
    
    # ax.set_xlim([-0.3, 1.6])
    # ax.set_ylim([-0.01, 0.015])
    ax.set_ylim([-0.5, 0.5])
    ax.set_xticklabels(['SFC(S1)', 'SFC(S2)', 'SFC(S3)'], rotation=30)
    ax.set_ylabel('Coefficient')
    ax.set_title('Memory S1', color=colors[0])
    
    #### plot rank-2
    ax = fig.add_subplot(1, 3, 2)
    sns.violinplot(data=df.loc[:,['Acc R2b1','Acc R2b2','Acc R2b3']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.scatter([0,1], [data['S2_R1_acc_beta'], data['S2_R2_acc_beta']], color = colors[1],
               marker='+', s=60, zorder=3)
    
    # ax.set_xlim([-0.3, 1.6])
    # ax.set_ylim([-0.01, 0.015])
    ax.set_ylim([-0.5, 0.5])
    ax.set_xticklabels(['SFC(S1)', 'SFC(S2)', 'SFC(S3)'], rotation=30)
    ax.set_yticklabels([])
    ax.set_ylabel(' ')
    ax.set_title('Memory S2', color=colors[1])
    
    #### plot rank-3
    ax = fig.add_subplot(1, 3, 3)
    sns.violinplot(data=df.loc[:,['Acc R3b1','Acc R3b2','Acc R3b3']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.scatter([2], [data['S3_R3_acc_beta']], color = colors[2],
               marker='+', s=60, zorder=3)
    
    # ax.set_xlim([-0.3, 1.6])
    # ax.set_ylim([-0.01, 0.015])
    ax.set_ylim([-0.5, 0.5])
    ax.set_xticklabels(['SFC(S1)', 'SFC(S2)', 'SFC(S3)'], rotation=30)
    ax.set_yticklabels([])
    ax.set_ylabel(' ')
    ax.set_title('Memory S3', color=colors[2])  
#%% spk memory, stepwise, L3 error trial, ocean
path = r'C:\Wen\Project\Osc_control\Data\monkey_combined\SFC_fw_error_L3'
file = '\\Ocean_SPK_Mem_FW_x_SFCacc_perm_stats_3Seq_StepLM_Err'

file = path + file
data = sio.loadmat(file)

with plt.style.context(style_path):
    fig = plt.figure(figsize = (3.6, 1.6), dpi = 300)
    ax = fig.add_subplot(1, 3, 1)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    
    #### groot
    df = pd.DataFrame({'Acc R1b1':data['beta_hold_S1'][0,:], 'Acc R1b2':0,'Acc R1b3':data['beta_hold_S1'][1,:],
                       'Acc R2b1':data['beta_hold_S2'][0,:], 'Acc R2b2':0,'Acc R2b3':data['beta_hold_S2'][1,:],
                       'Acc R3b1':data['beta_hold_S3'][0,:], 'Acc R3b2':data['beta_hold_S3'][1,:],'Acc R3b3':0})
    
    #### plot rank-1
    sns.violinplot(data=df.loc[:,['Acc R1b1','Acc R1b2','Acc R1b3']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.scatter([0,2], [data['S1_R1_acc_beta'],data['S1_R3_acc_beta']], color = colors[0],
               marker='+', s=60, zorder=3)
    # ax.scatter([1], [0], color = colors[0], linestyle='--',
    #            marker='+', s=60, zorder=3)
    
    # ax.set_xlim([-0.3, 1.6])
    # ax.set_ylim([-0.01, 0.015])
    ax.set_ylim([-0.5, 0.5])
    ax.set_xticklabels(['S1', 'S2', 'S3'], rotation=0)
    ax.set_xlabel('SFC strength')
    ax.set_ylabel('Coefficient')
    ax.set_title('Memory S1', color=colors[0])
    
    #### plot rank-2
    ax = fig.add_subplot(1, 3, 2)
    sns.violinplot(data=df.loc[:,['Acc R2b1','Acc R2b2','Acc R2b3']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.scatter([0,2], [data['S2_R1_acc_beta'], data['S2_R3_acc_beta']], color = colors[1],
               marker='+', s=60, zorder=3)
    
    # ax.set_xlim([-0.3, 1.6])
    # ax.set_ylim([-0.01, 0.015])
    ax.set_ylim([-0.5, 0.5])
    ax.set_xticklabels(['S1', 'S2', 'S3'], rotation=0)
    ax.set_yticklabels([])
    ax.set_ylabel(' ')
    ax.set_xlabel('SFC strength')
    ax.set_title('Memory S2', color=colors[1])
    
    #### plot rank-3
    ax = fig.add_subplot(1, 3, 3)
    sns.violinplot(data=df.loc[:,['Acc R3b1','Acc R3b2','Acc R3b3']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.scatter([0,1], [data['S3_R1_acc_beta'], data['S3_R2_acc_beta']], color = colors[2],
               marker='+', s=60, zorder=3)
    
    # ax.set_xlim([-0.3, 1.6])
    # ax.set_ylim([-0.01, 0.015])
    ax.set_ylim([-0.5, 0.5])
    ax.set_xticklabels(['S1', 'S2', 'S3'], rotation=0)
    ax.set_yticklabels([])
    ax.set_ylabel(' ')
    ax.set_xlabel('SFC strength')
    ax.set_title('Memory S3', color=colors[2])   
#%% spk memory, stepwise, L3 error trial, groot
path = r'C:\Wen\Project\Osc_control\Data\monkey_combined\SFC_fw_error_L3'
file = '\\Groot_SPK_Mem_FW_x_SFCacc_perm_stats_3Seq_StepLM_Err'

file = path + file
data = sio.loadmat(file)

with plt.style.context(style_path):
    fig = plt.figure(figsize = (3.6, 1.6), dpi = 300)
    ax = fig.add_subplot(1, 3, 1)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    
    #### groot
    df = pd.DataFrame({'Acc R1b1':data['beta_hold_S1'][0,:], 'Acc R1b2':0,'Acc R1b3':0,
                       'Acc R2b1':data['beta_hold_S2'][0,:], 'Acc R2b2':data['beta_hold_S2'][1,:],'Acc R2b3':0,
                       'Acc R3b1':0, 'Acc R3b2':0,'Acc R3b3':data['beta_hold_S3'][0,:]})
    
    #### plot rank-1
    sns.violinplot(data=df.loc[:,['Acc R1b1','Acc R1b2','Acc R1b3']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.scatter([0], [data['S1_R1_acc_beta']], color = colors[0],
               marker='+', s=60, zorder=3)
    # ax.scatter([1], [0], color = colors[0], linestyle='--',
    #            marker='+', s=60, zorder=3)
    
    # ax.set_xlim([-0.3, 1.6])
    # ax.set_ylim([-0.01, 0.015])
    ax.set_ylim([-0.5, 0.5])
    ax.set_xticklabels(['S1', 'S2', 'S3'], rotation=0)
    ax.set_ylabel('Coefficient')
    ax.set_xlabel('SFC strength')
    ax.set_title('Memory S1', color=colors[0])
    
    #### plot rank-2
    ax = fig.add_subplot(1, 3, 2)
    sns.violinplot(data=df.loc[:,['Acc R2b1','Acc R2b2','Acc R2b3']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.scatter([0,1], [data['S2_R1_acc_beta'], data['S2_R2_acc_beta']], color = colors[1],
               marker='+', s=60, zorder=3)
    
    # ax.set_xlim([-0.3, 1.6])
    # ax.set_ylim([-0.01, 0.015])
    ax.set_ylim([-0.5, 0.5])
    ax.set_xticklabels(['S1', 'S2', 'S3'], rotation=0)
    ax.set_yticklabels([])
    ax.set_ylabel(' ')
    ax.set_xlabel('SFC strength')
    ax.set_title('Memory S2', color=colors[1])
    
    #### plot rank-3
    ax = fig.add_subplot(1, 3, 3)
    sns.violinplot(data=df.loc[:,['Acc R3b1','Acc R3b2','Acc R3b3']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.scatter([2], [data['S3_R3_acc_beta']], color = colors[2],
               marker='+', s=60, zorder=3)
    
    # ax.set_xlim([-0.3, 1.6])
    # ax.set_ylim([-0.01, 0.015])
    ax.set_ylim([-0.5, 0.5])
    ax.set_xticklabels(['S1', 'S2', 'S3'], rotation=0)
    ax.set_yticklabels([])
    ax.set_ylabel(' ')
    ax.set_xlabel('SFC strength')
    ax.set_title('Memory S3', color=colors[2])   
#%% spk memory, linear moldel, L3 error trial swap(1-2-3)
path = r'C:\Wen\Project\Osc_control\Data\monkey_combined\SFC_fw_error_L3'
file = '\\Ocean_SPK_Mem_FW_x_SFCacc_perm_stats_3Seq_StepLM_Err.mat'
# file = '\\Groot_SPK_Mem_FW_x_SFCacc_perm_stats_3Seq_StepLM_Err.mat'

file = path + file
data = sio.loadmat(file)

with plt.style.context(style_path):
    fig = plt.figure(figsize = (4.8, 1.7), dpi = 300)
    ax = fig.add_subplot(1, 3, 1)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    
    df = pd.DataFrame({'Acc R1b1':data['beta_hold_S1'][0,:], 'Acc R1b2':data['beta_hold_S1'][1,:],'Acc R1b3':data['beta_hold_S1'][2,:],
                       'Acc R2b1':data['beta_hold_S2'][0,:], 'Acc R2b2':data['beta_hold_S2'][1,:],'Acc R2b3':data['beta_hold_S2'][2,:],
                       'Acc R3b1':data['beta_hold_S3'][0,:], 'Acc R3b2':data['beta_hold_S3'][1,:],'Acc R3b3':data['beta_hold_S3'][2,:]})
    
    #### plot rank-1
    sns.violinplot(data=df.loc[:,['Acc R1b1','Acc R1b2','Acc R1b3']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.plot([0,1,2], [data['S1_R1_acc_beta'][0], data['S1_R2_acc_beta'][0], data['S1_R3_acc_beta'][0]], color='gray')
    ax.scatter([0,1,2], [data['S1_R1_acc_beta'], data['S1_R2_acc_beta'], data['S1_R3_acc_beta']], color = colors[0],
               marker='+', s=60, zorder=3)
    
    ax.set_xlim([-0.6, 2.6])
    ax.set_ylim([-0.5, 0.75])
    # ax.set_ylim([-2.3, 2.3])
    ax.set_xticklabels(['S1', 'S2', 'S3'], rotation=0)
    ax.set_xlabel('SFC strength')
    ax.set_ylabel('Coefficient')
    # ax.set_title('Rank-1')
    
    #### plot rank-2
    ax = fig.add_subplot(1, 3, 2)
    sns.violinplot(data=df.loc[:,['Acc R2b1','Acc R2b2','Acc R2b3']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.plot([0,1,2], [data['S2_R1_acc_beta'][0], data['S2_R2_acc_beta'][0], data['S2_R3_acc_beta'][0]], color='gray')
    ax.scatter([0,1,2], [data['S2_R1_acc_beta'], data['S2_R2_acc_beta'], data['S2_R3_acc_beta']], color = colors[1],
               marker='+', s=60, zorder=3)
    
    ax.set_xlim([-0.6, 2.6])
    ax.set_ylim([-0.5, 0.75])
    # ax.set_ylim([-2.3, 2.3])
    ax.set_xticklabels(['S1', 'S2', 'S3'], rotation=0)
    ax.set_xlabel('SFC strength')
    ax.set_ylabel('Coefficient')
    
    #### plot rank-3
    ax = fig.add_subplot(1, 3, 3)
    sns.violinplot(data=df.loc[:,['Acc R3b1','Acc R3b2','Acc R3b3']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.plot([0,1,2], [data['S3_R1_acc_beta'][0], data['S3_R2_acc_beta'][0], data['S3_R3_acc_beta'][0]], color='gray')
    ax.scatter([0,1,2], [data['S3_R1_acc_beta'], data['S3_R2_acc_beta'], data['S3_R3_acc_beta']], color = colors[2],
               marker='+', s=60, zorder=3)
    
    ax.set_xlim([-0.6, 2.6])
    ax.set_ylim([-0.5, 0.75])
    # ax.set_ylim([-2.3, 2.3])
    ax.set_xticklabels(['S1', 'S2', 'S3'], rotation=0)
    ax.set_xlabel('SFC strength')
    ax.set_ylabel('Coefficient')
    
    plt.tight_layout()
#%% spk memory, linear moldel, L3 error trial, swap(2-3)
path = r'C:\Wen\Project\Osc_control\Data\monkey_combined\SFC_fw_error_L3'
# file = '\\Groot_SPK_Mem_FW_x_SFCacc_perm_stats_3Seq_2pred_OrderErr.mat'
file = '\\Ocean_SPK_Mem_FW_x_SFCacc_perm_stats_3Seq_2pred_OrderErr.mat'

file = path + file
data = sio.loadmat(file)

with plt.style.context(style_path):
    fig = plt.figure(figsize = (3.2, 1.7), dpi = 300)
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    
    df = pd.DataFrame({'Acc R2b2':data['beta_hold_S2'][0,:], 'Acc R2b3':data['beta_hold_S2'][1,:],
                       'Acc R3b2':data['beta_hold_S2'][0,:], 'Acc R3b3':data['beta_hold_S2'][1,:]})
    
    #### plot rank-1

    
    #### plot rank-2
    ax = fig.add_subplot(1, 2, 1)
    sns.violinplot(data=df.loc[:,['Acc R2b2','Acc R2b3']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.plot([0,1], [ data['S2_R2_acc_beta'][0], data['S2_R3_acc_beta'][0]], color='gray')
    ax.scatter([0,1], [data['S2_R2_acc_beta'], data['S2_R3_acc_beta']], color = colors[1],
               marker='+', s=60, zorder=3)
    
    ax.set_xlim([-0.6, 1.6])
    ax.set_ylim([-0.5, 0.75])
    # ax.set_ylim([-2.3, 2.3])
    ax.set_xticklabels(['S2', 'S3'], rotation=0)
    ax.set_xlabel('SFC strength')
    ax.set_ylabel('Coefficient')
    
    #### plot rank-3
    ax = fig.add_subplot(1, 2, 2)
    sns.violinplot(data=df.loc[:,['Acc R3b2','Acc R3b3']], bw_adjust=.5, cut=0.5, linewidth=0, 
                   palette=['gray'])
    plt.setp(ax.collections, alpha=.5)
    
    ax.plot([0,1], [data['S3_R2_acc_beta'][0], data['S3_R3_acc_beta'][0]], color='gray')
    ax.scatter([0,1], [data['S3_R2_acc_beta'], data['S3_R3_acc_beta']], color = colors[2],
               marker='+', s=60, zorder=3)
    
    ax.set_xlim([-0.6, 1.6])
    ax.set_ylim([-0.5, 0.75])
    # ax.set_ylim([-2.3, 2.3])
    ax.set_xticklabels(['S2', 'S3'], rotation=0)
    ax.set_xlabel('SFC strength')
    ax.set_ylabel('Coefficient')
    
    plt.tight_layout()


