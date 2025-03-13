# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 19:29:19 2023

@author: Wen
"""

style_path = r'C:\Wen\OneDrive\CodeHub\Python\style_paper.mplstyle'

import numpy as np
import pandas as pd
import os

from scipy import stats
import statsmodels.api as sm
from statsmodels.formula.api import ols

import matplotlib.pyplot as plt
from matplotlib import ticker
#%%
def eta_squared(aov):
    aov['eta_sq'] = 'NaN'
    aov['eta_sq'] = aov[:-1]['sum_sq']/sum(aov['sum_sq'])
    return aov
def omega_squared(aov):
    mse = aov['sum_sq'][-1]/aov['df'][-1]
    aov['omega_sq'] = 'NaN'
    aov['omega_sq'] = (aov[:-1]['sum_sq']-(aov[:-1]['df']*mse))/(sum(aov['sum_sq'])+mse)
    return aov

def plot_activity(x, seqL, mkts_merge, title, savename):
    with plt.style.context(style_path):
        fig, axs = plt.subplots(1, 1, figsize=(2, 1.5), dpi=300) # , constrained_layout=True
        vmin, vmax = 0, 0.002  # set min and max ERDS values in plo
        # vmin, vmax = 15, 60 # set min and max ERDS values in plo
        f = axs.pcolormesh(x, vmin=vmin, vmax=vmax) #vmin=vmin, vmax=vmax, ,  cmap=plt.get_cmap('RdBu_r')
        axs.vlines(mkts_merge[0,1:3]-mkts_merge[1,0], 0, x.shape[0], color='w', linestyle = '--')
        
        axs.set_xticks(mkts_merge[0, 1:3]-mkts_merge[1, 0], ['S1', 'S2'])
        
        scientific = 0
        if scientific == 1:
            formatter = ticker.ScalarFormatter(useMathText=True)
            formatter.set_scientific(True) 
            formatter.set_powerlimits((-1, 1)) 
            axs.yaxis.set_major_formatter(formatter)
        axs.yaxis.set_label_text('Channel')
        
        cbar = fig.colorbar(f)
        scientific = 1
        if scientific == 1:
            cbar.formatter.set_scientific(True)
            cbar.formatter.set_powerlimits((-1, 1))
            cbar.formatter.set_useMathText(True)
        
def plot_pev(pev, seqL, mkts_merge, title):
    with plt.style.context(style_path):
        fig, axs = plt.subplots(1, seqL, figsize=(3, 1.5), dpi=300, constrained_layout=True)
        for i in range(seqL):
            f = axs[i].pcolormesh(pev[:,:,i], vmin=-0, vmax=0.14,
                              cmap=plt.get_cmap('viridis'))
            axs[i].vlines(mkts_merge[0,1:3]-mkts_merge[1,0], 0, pev.shape[0], color='w', linestyle='--')
            axs[i].vlines(mkts_merge[1,1:3]-mkts_merge[1,0], 0, pev.shape[0], color='w', linestyle = '-', linewidth=2)
            axs[i].set_title('Rank' + str(i+1))

        axs[0].set_ylabel('Channel')
        # plt.subplots_adjust(left=0, bottom=0, right=1, top=0.85, hspace=0.1, wspace=0.15)
        plt.suptitle(title)
        fig.colorbar(f)
        # plt.margins(8, 0.5)
       
def plot_pval(pval, seqL, mkts_merge, title, savename):
    x = pval < 0.05
    sort = 0
    with plt.style.context(style_path):
        fig, axs = plt.subplots(1, seqL, figsize=(3.5, 1.5), dpi=300, constrained_layout=True)
        for i in range(seqL):
            xx = x[:,:,i]
            
            # Sort the rows based on the column indices
            if sort == 1:
                idx = np.argmax(xx, axis=1)
                idx[idx == 0] = xx.shape[1]
                idx = np.argsort(idx)[::-1]
                xx = xx[idx,:]
            
            f = axs[i].pcolormesh(xx, vmin=-0, vmax=1,
                              cmap=plt.get_cmap('viridis'))
            axs[i].vlines(mkts_merge[0,1:3]-mkts_merge[1,0], 0, x.shape[0], color='w', linestyle='--')
            axs[i].set_xticks(mkts_merge[0,1:3]-mkts_merge[1,0], ['S1', 'S2'])
            axs[i].set_title(' ')
            
            formatter = ticker.ScalarFormatter(useMathText=True)
            axs[i].yaxis.set_major_formatter(formatter)

        axs[0].set_ylabel('Channel')
        
def plot_pev_mask(pval, seqL, mkts_merge, title, savename):
    x = pval
    sort = 0
    with plt.style.context(style_path):
        fig, axs = plt.subplots(1, seqL, figsize=(3.5, 1.5), dpi=300, constrained_layout=True)
        for i in range(seqL):
            xx = x[:,:,i]
            
            # Sort the rows based on the column indices
            if sort == 1:
                idx = np.argmax(xx, axis=1)
                idx[idx == 0] = xx.shape[1]
                idx = np.argsort(idx)[::-1]
                xx = xx[idx,:]
            
            f = axs[i].pcolormesh(xx, vmin=-0, vmax=0.14,
                              cmap=plt.get_cmap('viridis'))
            axs[i].vlines(mkts_merge[0, 1:3]-mkts_merge[1,0], 0, x.shape[0], color='w', linestyle='--')
            axs[i].set_xticks(mkts_merge[0, 1:3]-mkts_merge[1,0], ['S1', 'S2'])
            axs[i].set_title(' ')
            
            formatter = ticker.ScalarFormatter(useMathText=True)
            axs[i].yaxis.set_major_formatter(formatter)

        axs[0].set_ylabel('Channel')
        fig.colorbar(f, orientation='vertical')
        
def plot_pev_single_region(x, ch_sit, mkts_merge):
    x = pev.copy()
    areas = ['PFC','PM','PAR']
    for area in areas:
        sel_area = np.where(ch_sit['area']==area)[0]
        xx = x[sel_area,:]
        with plt.style.context(style_path):
            fig = plt.figure(figsize=(4, 4), dpi=300)
            ax = fig.add_subplot(111)
            f = ax.pcolormesh(xx, vmin=-0, vmax=0.12, cmap=plt.get_cmap('viridis')) #vmin=vmin, vmax=vmax,
            ax.vlines(mkts_merge[0], 0, xx.shape[0], color='w', linestyle = '--')
            fig.colorbar(f)
            ax.set_ylabel('Channel')
            
def pev(x, df, seqL):
    pev = np.zeros((x.shape[0], x.shape[2], 2))
    pval = np.zeros((x.shape[0], x.shape[2], 2))
    formular = 'fr ~ C(Target_1) + C(Target_2)'
    if seqL==3:
        pev = np.zeros((x.shape[0], x.shape[2], 3))
        pval = np.zeros((x.shape[0], x.shape[2], 3))
        formular = 'fr ~ C(Target_1) + C(Target_2) + C(Target_3)'

    for neuid in range(pev.shape[0]):
        if np.mod(neuid, 10) == 0:
            print(neuid)
        for time in range(pev.shape[1]):
            df['fr'] = x[neuid, :, time]
            # ols model
            model = ols(formular, df, missing='drop').fit()
            if len(list(model.conf_int().index)) > 2:
                aov_table = sm.stats.anova_lm(model, typ=2)
                # eta_squared(aov_table)
                aov_table = omega_squared(aov_table)
                for i in range(seqL):
                    pev[neuid, time, i] = aov_table['omega_sq'][i]  
                    pval[neuid, time, i] = aov_table['PR(>F)'][i]
    return pev, pval

def pev_shuffle_spk(x, df, seqL):
    pevsf = np.zeros((x.shape[0], 1000, 2))
    formular = 'fr ~ C(Target_1) + C(Target_2)'
    if seqL==3:
        pevsf = np.zeros((x.shape[0], 1000, 3))
        formular = 'fr ~ C(Target_1) + C(Target_2) + C(Target_3)'
        
    for neuid in range(pevsf.shape[0]):
        if np.mod(neuid, 10) == 0:
            print(neuid)
        df['fr'] = x[neuid, :, 15:20].mean(1)  # baseline FR
        for per in range(1000):
            peridx = np.random.choice(np.arange(df.shape[0]), 200, replace=True)
            dftmp = df.loc[peridx, :]
            
            model = ols(formular, dftmp, missing='drop').fit()
            aov_table = sm.stats.anova_lm(model, typ=2)
            aov_table = omega_squared(aov_table)
            
            for i in range(seqL):
                pevsf[neuid, per, i] = aov_table['omega_sq'][i]
    return pevsf

def pev_shuffle_lfp(x, df, seqL):
    pevsf = np.zeros((x.shape[0], 1000, 2))
    formular = 'fr ~ C(Target_1) + C(Target_2)'
    if seqL==3:
        pevsf = np.zeros((x.shape[0], 1000, 3))
        formular = 'fr ~ C(Target_1) + C(Target_2) + C(Target_3)'
        
    for neuid in range(pevsf.shape[0]):
        if np.mod(neuid, 10) == 0:
            print(neuid)
        df['fr'] = x[neuid, :, 15:20].mean(1)  # baseline FR
        for per in range(1000):
            peridx = np.random.choice(np.arange(df.shape[0]), 200, replace=True)
            dftmp = df.loc[peridx, :]
            
            model = ols(formular, dftmp, missing='drop').fit()
            aov_table = sm.stats.anova_lm(model, typ=2)
            aov_table = omega_squared(aov_table)
            
            for i in range(seqL):
                pevsf[neuid, per, i] = aov_table['omega_sq'][i]
    return pevsf

def find_pev_sig_time(x, p_alpha):
    window = 3
    pvalmask = (x < p_alpha).astype(np.int32)
    pev_sig = pd.DataFrame(pvalmask).rolling(window=window, axis=1).sum()
    r, c = np.where(pev_sig==3)
    pev_sig = np.zeros(pev_sig.shape)
    for i in range(len(r)):
        pev_sig[r[i], c[i]-window+1:c[i]+1] = 1
    return pev_sig

def permu_test(pev, pevsf, seqL):
    pval = np.zeros(pev.shape)
    for neu in range(pev.shape[0]):
        x = pev[neu, :, 0]
        for t in range(len(x)):
            pval[neu, t, 0] = sum(pevsf[neu, :, 0] > x[t])/1000

    for neu in range(pev.shape[0]):
        x = pev[neu, :, 1]
        for t in range(len(x)):
            pval[neu, t, 1] = sum(pevsf[neu, :, 1] > x[t])/1000
            
    if seqL == 3:
        for neu in range(pev.shape[0]):
            x = pev[neu, :, 2]
            for t in range(len(x)):
                pval[neu, t, 2] = sum(pevsf[neu, :, 2] > x[t])/1000
    
    return pval

def find_sig_to_baseline(x, baseline):
    pval = np.zeros([2, x.shape[1]])
    for t in range(x.shape[1]):
        pval[:,t] = stats.ttest_rel(x[:,t], baseline)
    return pval











    