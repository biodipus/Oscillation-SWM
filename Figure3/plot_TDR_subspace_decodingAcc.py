# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 15:42:21 2024

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
import itertools

from matplotlib.patches import Rectangle
from matplotlib import ticker

import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.formula.api import ols
import statsmodels.stats.multitest as smm
from scipy.stats import wilcoxon, mannwhitneyu
# from glmcc import GLMCC  # model
# from glmcc import spiketime_relative
# import Spike_contrast as SC

from utils import fun_pev
from utils import fun_loadData

from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn import svm
#%% functions
def corr_spikes_sorted(spike1, spike2, tbin, tau_max, resolution):
    tau_max_i = int(tau_max / resolution)
    tbin_i = int(tbin / resolution)
    cross = np.zeros(int(2 * tau_max_i / tbin_i + 1), 'd')
    j0 = 0
    for spki in spike1:
        j = j0
        while j < len(spike2) and spike2[j] - spki < -tau_max_i - tbin_i / 2.0:
            j += 1
        j0 = j
        while j < len(spike2) and spike2[j] - spki < tau_max_i + tbin_i / 2.0:
            cross[int(
                (spike2[j] - spki + tau_max_i + 0.5 * tbin_i) / tbin_i)] += 1.0
            j += 1
    return cross

def merge_vectors_to_matrix(vectors):
    max_length = max(len(vec) for vec in vectors)
    matrix = np.zeros((len(vectors), max_length))
    for i, vec in enumerate(vectors):
        matrix[i, :len(vec)] = vec
    return matrix

def prob_targroup(x, y):
    prob = np.zeros([x.shape[0], 2, 100])+np.nan
    for r in range(100):
        knn_classifier = KNeighborsClassifier(n_neighbors=6)
        rkf = StratifiedKFold(n_splits=2, shuffle=True, random_state=None)
        k=0
        for i, (train_idx, test_idx) in enumerate(rkf.split(x[:,:], y)):
        # for train_idx, test_idx in rkf.split(x):
            X_train, X_test = x[train_idx], x[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]
            knn_classifier.fit(X_train, y_train)
        
            probabilities = knn_classifier.predict_proba(X_test)
            aa = np.array([probabilities[i, int(y_test[i]-1)] for i in range(len(probabilities))])
            prob[test_idx, k, r] = aa
            k += 1
    prob = np.nanmean(prob,1)
    prob = np.nanmean(prob,1)
    return prob

def prob_targroup_svm(xfix, x, y, shuffle):
    from sklearn.metrics import accuracy_score
    # prob = np.zeros([x.shape[0], 3, 1000])+np.nan
    # print(['shuffle:', shuffle])
    train_r = 0
    test_r = 0
    
    #### normalize X
    scaler = StandardScaler()
    xfix = scaler.fit_transform(xfix).T
    scaler = StandardScaler()
    x = scaler.fit_transform(x).T
    # xfix = xfix.T
    # x = x.T
    
    y = y.reshape(-1,1)
    
    C = 1
    gamma = 1
    n_fold = 2 #y.shape[0]
    n_repeats = 5
    # print([n_fold, n_repeats])

    clf = svm.SVC(kernel='linear', C=C, gamma=gamma, probability=True)
    scores =  np.full((n_repeats, n_fold), np.nan)
    
    for r in range(n_repeats):
        # if np.mod(r,10) == 0:
        #     print(r)
        if shuffle == 1:
            np.random.shuffle(y[:,0])
            
        # rkf = StratifiedKFold(n_splits=n_fold, shuffle=True, random_state=None)  # for within rank
        rkf = KFold(n_splits=n_fold, shuffle=True, random_state=None)  # for cross rank
        for i, (train, test) in enumerate(rkf.split(x[:,:].T, y[:,train_r])):
            y_train = y[train, train_r]
            x_train = xfix[:,train].T
                    
            clf.fit(x_train, y_train)
            
            x_test = x[:,test].T
            y_test = y[test, test_r]
            
            # 预测类别
            # probabilities = clf.predict_proba(x_test)
            # predictions = clf.predict(x_test)
            # aa = np.array([probabilities[i, int(y_test[i]-1)] for i in range(probabilities.shape[0])])
            # prob[test, i, r] = aa
            
            scores[r,i] = clf.score(x_test, y_test)

    # prob = np.nanmean(prob,1)
    # prob = np.nanmean(prob,1)
    scores = scores.mean(1)
    # scores = scores.mean(0).mean(0)
    return scores

def find_sig(A, l):
    # 用于存储连续出现三个及以上0的索引
    result_indices = []
    count = 0
    zero_indices = []
    
    for index in range(len(A)):
        if A[index] == 1:
            count += 1
            zero_indices.append(index)
        else:
            # 如果有连续的 0，检查是否达到三个或以上
            if count >= l:
                result_indices.extend(zero_indices[-count:])  # 记录这部分索引
            count = 0
            zero_indices = []
    
    # 检查最后一段
    if count >= 3:
        result_indices.extend(zero_indices[-count:])
    return result_indices
# %% 
monkey = 'grootsw'
#### ocean
if 'ocean' == monkey:
    days = ['0709','0710','0713','0714','0715','0720','0721','0722','0723','0726','0727','0729','0730']
    # days = ['0710']
    fpath = r'C:\Wen\Project\Osc_control\Data\Data_ocean'
    seqL = 2
#### groot switch
if 'grootsw' == monkey:
    days = ['04-15','05-07','05-09','05-13','05-14','05-15','05-16','05-20','05-21','05-22'] # '03-23',,'03-27'
    # days = ['03-11'] 
    fpath = r'C:\Wen\Project\Osc_control\Data\Data_grootsw'
    seqL = 2
#%% load matlab decoding result
# osc_band = ['4_8','80_120']
# bandid = 1
select_rule = 0
region = ['PFC','PM','frontal'] #, 'PM', 'PAR', 'M1']

infomative_all = [None] * len(days)
infomative_all_shuffle = [None] * len(days)
info_all = [None] * len(days)
for d in range(len(days)):
    print(d)
    acc_all = [None] * 12
    acc_all_sf = [None] * 12
    if 'ocean' == monkey:
        # SPK
        fpath = r'C:\Wen\Project\Osc_control\Data\Data_ocean\SPK\FR_segment_win100step50'
        file = fpath + '\\'+days[d]+'_fr_L2_win100step50_selChs.pkl'
        info, neu_sit_info, fr_arr, mkts_merge, spk_arr, mktime = fun_loadData.load_fr_data(file, region, select_rule)
        
        #### spk
        tdrpath = r'C:\Wen\Project\Osc_control\Data\Result_ocean\SPK\tdr_frontal_L2_croval\tdr_sens'
        tdrfile = tdrpath + '\\' + days[d]+'_weight_enc250_CV.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3).mean(2)
        acc_all[0] = acc[:,[0]]
        acc_all[1] = acc[:,[1]]
        tdrfile = tdrpath + '\\' + days[d]+'_weight_enc250_CV_shuffle.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3)
        acc_all_sf[0] = acc[:,[0]]
        acc_all_sf[1] = acc[:,[1]]
        
        tdrpath = r'C:\Wen\Project\Osc_control\Data\Result_ocean\SPK\tdr_frontal_L2_croval\tdr_mem_delay1'
        tdrfile = tdrpath + '\\' + days[d]+'_weight_mem_delay_CV.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3).mean(2)
        acc_all[2] = acc[:,[0]]
        acc_all[3] = acc[:,[1]]
        tdrfile = tdrpath + '\\' + days[d]+'_weight_mem_delay_CV_shuffle.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3)
        acc_all_sf[2] = acc[:,[0]]
        acc_all_sf[3] = acc[:,[1]]
        
        #### theta
        tdrpath = r'C:\Wen\Project\Osc_control\Data\Result_ocean\LFP\tdr_frontal_L2_croval\tdr_sens'
        tdrfile = tdrpath + '\\' + days[d]+'_weight_enc250_4_8_CV.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3).mean(2)
        acc_all[4] = acc[:,[0]]
        acc_all[5] = acc[:,[1]]
        tdrfile = tdrpath + '\\' + days[d]+'_weight_enc250_4_8_CV_shuffle.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3)
        acc_all_sf[4] = acc[:,[0]]
        acc_all_sf[5] = acc[:,[1]]
        
        tdrpath = r'C:\Wen\Project\Osc_control\Data\Result_ocean\LFP\tdr_frontal_L2_croval\tdr_mem_delay1'
        tdrfile = tdrpath + '\\' + days[d]+'_weight_mem_delay_4_8_CV.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3).mean(2)
        acc_all[6] = acc[:,[0]]
        acc_all[7] = acc[:,[1]]
        tdrfile = tdrpath + '\\' + days[d]+'_weight_mem_delay_4_8_CV_shuffle.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3)
        acc_all_sf[6] = acc[:,[0]]
        acc_all_sf[7] = acc[:,[1]]
        
        #### gamma
        tdrpath = r'C:\Wen\Project\Osc_control\Data\Result_ocean\LFP\tdr_frontal_L2_croval\tdr_sens'
        tdrfile = tdrpath + '\\' + days[d]+'_weight_enc250_80_120_CV.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3).mean(2)
        acc_all[8] = acc[:,[0]]
        acc_all[9] = acc[:,[1]]
        tdrfile = tdrpath + '\\' + days[d]+'_weight_enc250_80_120_CV_shuffle.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3)
        acc_all_sf[8] = acc[:,[0]]
        acc_all_sf[9] = acc[:,[1]]
        
        tdrpath = r'C:\Wen\Project\Osc_control\Data\Result_ocean\LFP\tdr_frontal_L2_croval\tdr_mem_delay1'
        tdrfile = tdrpath + '\\' + days[d]+'_weight_mem_delay_80_120_CV.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3).mean(2)
        acc_all[10] = acc[:,[0]]
        acc_all[11] = acc[:,[1]]
        tdrfile = tdrpath + '\\' + days[d]+'_weight_mem_delay_80_120_CV_shuffle.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3)
        acc_all_sf[10] = acc[:,[0]]
        acc_all_sf[11] = acc[:,[1]]
        
        
        acc_all = np.concatenate(acc_all, axis=1)
        acc_all = acc_all.T
        # acc_all = np.concatenate((acc_all[15:45,:], acc_all[50:75,:]), axis=0).T
        infomative_all[d] = acc_all
        
        acc_all_sf = np.concatenate(acc_all_sf, axis=1) # time*type*shuffle_rep
        infomative_all_shuffle[d] = acc_all_sf
    
    if 'grootsw' == monkey:
        fpath = r'C:\Wen\Project\Osc_control\Data\Data_grootsw\SPK_sw\FR_segment_win100step50_frontal'
        file = fpath + '\\'+days[d]+'_fr_L2_win100step50_selChs.pkl'
        info, neu_sit_info, fr_arr, mkts_merge, spk_arr, mktime = fun_loadData.load_fr_data(file, region, select_rule)

        #### spk
        tdrpath = r'C:\Wen\Project\Osc_control\Data\Result_grootsw\SPK_sw\tdr_frontal_L2_croval\tdr_sens'
        tdrfile = tdrpath + '\\' + days[d]+'_weight_enc250_CV.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3).mean(2)
        acc_all[0] = acc[:,[0]]
        acc_all[1] = acc[:,[1]]
        tdrfile = tdrpath + '\\' + days[d]+'_weight_enc250_CV_shuffle.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3)
        acc_all_sf[0] = acc[:,[0]]
        acc_all_sf[1] = acc[:,[1]]
        
        tdrpath = r'C:\Wen\Project\Osc_control\Data\Result_grootsw\SPK_sw\tdr_frontal_L2_croval\tdr_mem_delay1'
        tdrfile = tdrpath + '\\' + days[d]+'_weight_mem_delay_CV.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3).mean(2)
        acc_all[2] = acc[:,[0]]
        acc_all[3] = acc[:,[1]]
        tdrfile = tdrpath + '\\' + days[d]+'_weight_mem_delay_CV_shuffle.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3)
        acc_all_sf[2] = acc[:,[0]]
        acc_all_sf[3] = acc[:,[1]]
        
        #### theta
        tdrpath = r'C:\Wen\Project\Osc_control\Data\Result_grootsw\LFP_sw\tdr_frontal_L2_croval\tdr_sens'
        tdrfile = tdrpath + '\\' + days[d]+'_weight_enc250_4_8_CV.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3).mean(2)
        acc_all[4] = acc[:,[0]]
        acc_all[5] = acc[:,[1]]
        tdrfile = tdrpath + '\\' + days[d]+'_weight_enc250_4_8_CV_shuffle.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3)
        acc_all_sf[4] = acc[:,[0]]
        acc_all_sf[5] = acc[:,[1]]
        
        tdrpath = r'C:\Wen\Project\Osc_control\Data\Result_grootsw\LFP_sw\tdr_frontal_L2_croval\tdr_mem_delay1'
        tdrfile = tdrpath + '\\' + days[d]+'_weight_mem_delay_4_8_CV.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3).mean(2)
        acc_all[6] = acc[:,[0]]
        acc_all[7] = acc[:,[1]]
        tdrfile = tdrpath + '\\' + days[d]+'_weight_mem_delay_4_8_CV_shuffle.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3)
        acc_all_sf[6] = acc[:,[0]]
        acc_all_sf[7] = acc[:,[1]]
        
        #### gamma
        tdrpath = r'C:\Wen\Project\Osc_control\Data\Result_grootsw\LFP_sw\tdr_frontal_L2_croval\tdr_sens'
        tdrfile = tdrpath + '\\' + days[d]+'_weight_enc250_80_120_CV.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3).mean(2)
        acc_all[8] = acc[:,[0]]
        acc_all[9] = acc[:,[1]]
        tdrfile = tdrpath + '\\' + days[d]+'_weight_enc250_80_120_CV_shuffle.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3)
        acc_all_sf[8] = acc[:,[0]]
        acc_all_sf[9] = acc[:,[1]]
        
        tdrpath = r'C:\Wen\Project\Osc_control\Data\Result_grootsw\LFP_sw\tdr_frontal_L2_croval\tdr_mem_delay1'
        tdrfile = tdrpath + '\\' + days[d]+'_weight_mem_delay_80_120_CV.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3).mean(2)
        acc_all[10] = acc[:,[0]]
        acc_all[11] = acc[:,[1]]
        tdrfile = tdrpath + '\\' + days[d]+'_weight_mem_delay_80_120_CV_shuffle.mat'
        acc = sio.loadmat(tdrfile)['validationAccuracy'].mean(3)
        acc_all_sf[10] = acc[:,[0]]
        acc_all_sf[11] = acc[:,[1]]
       
        #### grootsw, LFP的数据取了整个trial
        for ii in range(4,12):
            acc_all[ii] = acc_all[ii][:85,:]
            acc_all_sf[ii] = acc_all_sf[ii][:85,:]
            
        acc_all = np.concatenate(acc_all, axis=1)
        acc_all = acc_all.T
        # acc_all = np.concatenate((acc_all[:,:45], acc_all[:,55:]), axis=1)
        infomative_all[d] = acc_all
        
        acc_all_sf = np.concatenate(acc_all_sf, axis=1) # time*type*shuffle_rep
        # acc_all_sf = np.concatenate((acc_all_sf[:45,:,:], acc_all_sf[55:,:,:]), axis=0)
        infomative_all_shuffle[d] = acc_all_sf
        
    if 'ocean' == monkey:
        infomative_ocean = infomative_all
        infomative_ocean_shuffle = infomative_all_shuffle
    if 'grootsw' == monkey:
        infomative_grootsw = infomative_all
        infomative_grootsw_shuffle = infomative_all_shuffle
#%% calculate mean alpha of all days
time = acc_all.shape[1]
dayNum = len(days)
shuffle_acc = [None] * 12
top_5_percent = np.zeros([12,time,dayNum])
first_sig_time = np.zeros([12,time,dayNum])+np.nan
for i in range(12): #
    # shuffle_acc[i] = np.concatenate([infomative_all_shuffle[d][:,i:i+1,:] for d in range(dayNum)], axis=2)
    shuffle_acc[i] = np.concatenate([infomative_all_shuffle[d][:,i:i+1,:] for d in range(dayNum)], axis=1)
    # shuffle_acc[i] = shuffle_acc[i].mean(1)
    for t in range(time):
        for d in range(dayNum):
            top_5_percent[i,t,d] = np.percentile(shuffle_acc[i][t,d,:], 99)
            
            mean = np.mean(shuffle_acc[i][t,d,:])
            std_dev = np.std(shuffle_acc[i][t,d,:], ddof=1)  # ddof=1 表示样本标准差
            # 计算99%的置信区间
            confidence_level = 0.99
            z_score = stats.norm.ppf((1 + confidence_level) / 2)  # 99% 置信区间对应的z分数
            # 下限和上限
            lower_bound = mean - z_score * std_dev
            upper_bound = mean + z_score * std_dev
            first_sig_time[i,t,d] = infomative_all[d][i,t] > upper_bound
            
top_5_percent = top_5_percent.mean(2)
# top_5_percent = top_5_percent[:,:,[0,1,2,4,5,6,7,8,9,10]].mean(2)
#%% entry/mem subspace sigtime
plotTime = np.arange(15,55)
titles = ['SPK']*4 + ['Theta']*4 + ['Gamma']*4
dayNum = len(infomative_all)

x = [None] * 12
pval_all = np.zeros([len(plotTime), 12])
for i in range(12): # spk, theta, gamma 3*4(E1,E2,M1,M2)=12
    x[i] = np.concatenate([infomative_all[d][i:i+1,:] for d in range(dayNum)], axis=0)
    baseline = x[i][:,:5].mean(1) #reshape([-1])
    shuffle_acc_tmp = top_5_percent[i,plotTime]
    x[i] = x[i][:, plotTime]
    for t in range(len(plotTime)):
        tmp = x[i][:,t]
        pval_all[t, i] = tmp.mean() < shuffle_acc_tmp[t]
        

pvals_corrected = pval_all#*np.sqrt(time)
sig_time1 = pvals_corrected<0.05

sig_time = np.ones_like(sig_time1)
for i in range(sig_time.shape[1]):
    result_indices = find_sig(sig_time1[:,i],3)
    sig_time[result_indices,i] = 0
sig_time = ~sig_time
#%% temp subspace sig time
infomative_all = infomative_ocean.copy()
# infomative_all = infomative_ocean + infomative_groot
titles = ['Transient']*3
dayNum = len(infomative_all)
time = 30 #infomative_all[0].shape[1]

x = [None] * 12
pval_all = np.zeros([time, 12])
for i in range(12): # spk, theta, gamma 3*4(E1,E2,M1,M2)=12
    x[i] = np.concatenate([infomative_all[d][12+i:12+i+1,:] for d in range(dayNum)], axis=0)
    # baseline = x[i][:,15:20].reshape([-1]) #
    # shuffle_acc_tmp = shuffle_acc[i][45:75,:]
    shuffle_acc_tmp = top_5_percent[12+i,45:75]
    x[i] = x[i][:,45:75]
    for t in range(time):
        tmp = x[i][:,t]
        # tval, pval_all[t, i] = stats.ttest_rel(baseline, tmp, alternative='less')
        # tval, pval_all[t, i] = stats.ttest_ind(baseline, tmp, alternative='less')
        # tval, pval_all[t, i] = stats.ttest_ind(tmp, 1/6, alternative='greater')
        # tval, pval_all[t, i] = stats.ttest_ind(baseline, tmp.mean(), alternative='less')
        # pval_all[t, i] = sum(shuffle_acc_tmp[t,:]>tmp.mean())/shuffle_acc_tmp.shape[1]
        pval_all[t, i] = tmp.mean() < shuffle_acc_tmp[t]
        
    # pvals_corrected = pval_all[:, i]*np.sqrt(time)
    # sig_time = pvals_corrected<0.001
    
#### correct pvals
pvals_corrected = pval_all.copy()
# for i in range(12):
    # reject, pvals_corrected[:,i], _, _ = smm.multipletests(pval_all[:,i], alpha=0.05, method='fdr_bh')

pvals_corrected = pval_all#*np.sqrt(time)
sig_time = pvals_corrected<0.01
#%% plot entry/mem
with plt.style.context(style_path):
    for i in range(3):
        fig, axs = plt.subplots(1, 1, figsize=(3.6, 0.8), dpi=300) # , constrained_layout=True
        colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        axs.plot(x[i*4].mean(0), color = colors[0], linewidth=1.5)
        axs.plot(x[i*4+1].mean(0), color = colors[1], linewidth=1.5)
        # axs[0].set_title(str(i+1))
        for d in range(dayNum):
            axs.plot(x[i*4][d,:], color=colors[0], alpha=0.3, linewidth=0.5)
            axs.plot(x[i*4+1][d,:], color=colors[1], alpha=0.3, linewidth=0.5)
        
        axs.set_xticks([4, 19, 29], ['S1', 'S2', 'Delay'])
        # ymax = x[i*4].mean(0).max() + 0.15
        # ymin = x[i*4].mean(0).min() - 0.05
        ymin = 0.1
        ymax = 0.41
        # if i>1:
        #     ymax = 0.58
        axs.set_ylim([ymin, ymax])
        axs.set_ylabel('Accuracy')
        # axs.set_title(titles[i])
        
        axs.add_patch(Rectangle((4, 0), 5, ymax, edgecolor=None, facecolor="gray", alpha=0.15))
        axs.add_patch(Rectangle((19, 0), 5, ymax, edgecolor=None, facecolor="gray", alpha=0.15))
        # axs.plot([29, 29], [ymin-0.001, ymax], color='w', linewidth=2, clip_on=False, zorder=100)
        
        #### plot sig line
        p_se1 = np.where(sig_time[:,i*4]==1)[0]
        p_se2 = np.where(sig_time[:,i*4+1]==1)[0]
        if len(p_se1) > 0:
            for s in range(len(p_se1)):
                axs.plot([p_se1[s], p_se1[s]+1], [ymax-0.02, ymax-0.02], color=colors[0])
        if len(p_se2) > 0:
            for s in range(len(p_se2)):
                axs.plot([p_se2[s], p_se2[s]+1], [ymax-0.04, ymax-0.04], color=colors[1])
            
        #### mem
        fig, axs = plt.subplots(1, 1, figsize=(3.6, 0.8), dpi=300) # , constrained_layout=True
        colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        axs.plot(x[i*4+2].mean(0), color = colors[0], linewidth=1.5)
        axs.plot(x[i*4+3].mean(0), color = colors[1], linewidth=1.5)
        # axs[0].set_title(str(i+1))
        for d in range(dayNum):
            axs.plot(x[i*4+2][d,:], color=colors[0], alpha=0.3, linewidth=0.5)
            axs.plot(x[i*4+3][d,:], color=colors[1], alpha=0.3, linewidth=0.5)
        
        axs.set_xticks([4, 19, 29], ['S1', 'S2', 'Delay'])
        # ymax = x[i*4+2].mean(0).max() + 0.15
        # ymin = x[i*4+2].mean(0).min() - 0.05
        ymin = 0.1
        ymax = 0.41
        # if i>1:
            # ymax = 0.58
        axs.set_ylim([ymin, ymax])
        axs.set_ylabel('Accuracy')
        # axs.set_title(titles[i])
        
        axs.add_patch(Rectangle((4, 0), 5, ymax, edgecolor=None, facecolor="gray", alpha=0.15))
        axs.add_patch(Rectangle((19, 0), 5, ymax, edgecolor=None, facecolor="gray", alpha=0.15))
        # axs.plot([29, 29], [ymin-0.001, ymax], color='w', linewidth=2, clip_on=False, zorder=100)
        #### plot sig line
        p_se1 = np.where(sig_time[:,i*4+2]==1)[0]
        p_se2 = np.where(sig_time[:,i*4+3]==1)[0]
        if len(p_se1) > 0:
            for s in range(len(p_se1)):
                axs.plot([p_se1[s], p_se1[s]+1], [ymax-0.02, ymax-0.02], color=colors[0])
        if len(p_se2) > 0:
            for s in range(len(p_se2)):
                axs.plot([p_se2[s], p_se2[s]+1], [ymax-0.04, ymax-0.04], color=colors[1])
#%% theta entry vs. gamma entry
# x = np.stack(infomative_ocean,axis=2)[:,15:55,:]
x = np.stack(infomative_grootsw,axis=2)[:,15:55,:]
spk_ent1 = x[0,7:9,:].mean(0)
spk_ent2 = x[1,22:24,:].mean(0)
spk_ent = np.concatenate((spk_ent1, spk_ent2))
spk_ent = (spk_ent1 + spk_ent2)/2

theta_ent1 = x[4,7:9,:].mean(0)
theta_ent2 = x[5,22:24,:].mean(0)
theta_ent = np.concatenate((theta_ent1, theta_ent2))
theta_ent = (theta_ent1 + theta_ent2)/2

gamma_ent1 = x[8,7:9,:].mean(0)
gamma_ent2 = x[9,22:24,:].mean(0)
gamma_ent = np.concatenate((gamma_ent1, gamma_ent2))
gamma_ent = (gamma_ent1 + gamma_ent2)/2

plt.scatter(theta_ent2, spk_ent2)
tval, pval_all = stats.pearsonr(spk_ent, theta_ent)

#### memory
theta_mem1 = x[6,-3:,:].mean(0)
theta_mem2 = x[7,-3:,:].mean(0)
# theta_mem = np.concatenate((theta_mem1, theta_mem2))
theta_mem = (theta_mem1 + theta_mem2)/2

gamma_mem1 = x[10, -3:, :].mean(0)
gamma_mem2 = x[11, -3:, :].mean(0)
# gamma_mem = np.concatenate((gamma_mem1, gamma_mem2))
gamma_mem = (gamma_mem1 + gamma_mem2)/2

tval, pval_all = stats.ttest_rel(theta_ent, gamma_ent)
tval, pval_all = stats.ttest_rel(gamma_ent, gamma_mem)

with plt.style.context(style_path):
    fig, axs = plt.subplots(1, 1, figsize=(1.2, 1.5), dpi=300) # , constrained_layout=True
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for i in range(len(theta_ent)):
        axs.plot([1.1, 1.9], [theta_ent[i], gamma_ent[i]], color='gray', linestyle='-', linewidth=1, alpha=0.3)
    axs.scatter(np.full(len(theta_ent), 1), theta_ent, color='k', edgecolor=[], alpha=0.8, zorder=1)
    axs.scatter(np.full(len(gamma_ent), 2), gamma_ent, color='k', edgecolor=[], alpha=0.8, zorder=1)
    axs.set_xlim([0.5,2.5])
    axs.set_ylim([0.12,0.6])
    axs.set_xticks([1,2], ['Theta', 'Gamma'])
    axs.set_ylabel('Accuracy')
    axs.set_title('Entry')

with plt.style.context(style_path):    
    fig, axs = plt.subplots(1, 1, figsize=(1.2, 1.5), dpi=300) # , constrained_layout=True
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for i in range(len(theta_ent)):
        axs.plot([1.1, 1.9], [theta_mem[i], gamma_mem[i]], color='gray', linestyle='-', linewidth=1, alpha=0.3)
    axs.scatter(np.full(len(theta_mem), 1), theta_mem, color='k', edgecolor=[], alpha=0.8, zorder=1)
    axs.scatter(np.full(len(gamma_mem), 2), gamma_mem, color='k', edgecolor=[], alpha=0.8, zorder=1)
    axs.set_xlim([0.5,2.5])
    axs.set_ylim([0.05,0.4])
    axs.set_xticks([1,2], ['Theta', 'Gamma'])
    axs.set_title('Memory')

#%% 计算memory subspace显著的最早时刻
def find_first_one_indices(matrix):
    # 将输入转换为numpy数组以便于操作
    matrix = np.array(matrix)
    # 创建一个列表来存储每一列第一个1的索引
    indices = []
    for col in range(matrix.shape[1]):
        # 使用argmax找出每一列第一个等于1的索引
        index = np.argmax(matrix[:, col] == 1)
        # 如果一列中没有1，则将index设为None
        if matrix[:, col][index] != 1:
            indices.append(None)
        else:
            indices.append(index)
    return indices

def find_last_one_indices(matrix):
    # 将输入转换为numpy数组以便于操作
    matrix = np.array(matrix)
    # 创建一个列表来存储每一列第一个1的索引
    indices = []
    for col in range(matrix.shape[1]):
        # 使用argmax找出每一列第一个等于1的索引
        index = np.where(matrix[:, col] == 1)[0][-1]
        # 如果一列中没有1，则将index设为None
        if matrix[:, col][index] != 1:
            indices.append(None)
        else:
            indices.append(index)
    return indices

all_sig_on = [None]*6 # [spk-e1, spk-e2, spk-m1, spk-m2, theta-e1, theta-e2]
all_sig_end = [None]*2 # [theta-e1, theta-e2]
#### spk entry
x = first_sig_time[:,15:55,:]
first_sig_time_m1 = x[0,:,:]
first_sig_time_m2 = x[1,:,:]

sig_time = np.ones_like(first_sig_time_m1)-1
for i in range(sig_time.shape[1]):
    result_indices = find_sig(first_sig_time_m1[:,i], 5)
    sig_time[result_indices,i] = 1

sig_time = find_first_one_indices(sig_time)
sig_time = np.array(sig_time)
all_sig_on[0] = (sig_time-3)*0.05

sig_time = np.ones_like(first_sig_time_m2)-1
for i in range(sig_time.shape[1]):
    result_indices = find_sig(first_sig_time_m2[:,i], 5)
    sig_time[result_indices,i] = 1

sig_time = find_first_one_indices(sig_time)
sig_time = np.array(sig_time)
all_sig_on[1] = (sig_time-18)*0.05

#### spk memory
x = first_sig_time[:,15:55,:]
first_sig_time_m1 = x[2,:,:]
first_sig_time_m2 = x[3,:,:]

sig_time = np.ones_like(first_sig_time_m1)-1
for i in range(sig_time.shape[1]):
    result_indices = find_sig(first_sig_time_m1[:,i], 5)
    sig_time[result_indices,i] = 1

sig_time = find_first_one_indices(sig_time)
sig_time = np.array(sig_time)
all_sig_on[2] = (sig_time-3)*0.05

sig_time = np.ones_like(first_sig_time_m2)-1
for i in range(sig_time.shape[1]):
    result_indices = find_sig(first_sig_time_m2[:,i], 5)
    sig_time[result_indices,i] = 1

sig_time = find_first_one_indices(sig_time)
sig_time = np.array(sig_time)
all_sig_on[3] = (sig_time-18)*0.05

#### theta
first_sig_time_m1 = x[4,:,:]
first_sig_time_m2 = x[5,:,:]

sig_time = np.ones_like(first_sig_time_m1)-1
for i in range(sig_time.shape[1]):
    result_indices = find_sig(first_sig_time_m1[:,i], 5)
    sig_time[result_indices,i] = 1
sig_time = find_first_one_indices(sig_time)
sig_time = np.array(sig_time)
all_sig_on[4] = (sig_time-3)*0.05

sig_time = np.ones_like(first_sig_time_m2)-1
for i in range(sig_time.shape[1]):
    result_indices = find_sig(first_sig_time_m2[:,i], 5)
    sig_time[result_indices,i] = 1
sig_time = find_first_one_indices(sig_time)
sig_time = np.array(sig_time)
all_sig_on[5] = (sig_time-18)*0.05

#### theta last sig
first_sig_time_m1 = x[4,:,:]
first_sig_time_m2 = x[5,:,:]

sig_time = np.ones_like(first_sig_time_m1)-1
for i in range(sig_time.shape[1]):
    result_indices = find_sig(first_sig_time_m1[:,i], 5)
    sig_time[result_indices,i] = 1
sig_time = find_last_one_indices(sig_time)
sig_time = np.array(sig_time)
all_sig_end[0] = (sig_time-3)*0.05

sig_time = np.ones_like(first_sig_time_m2)-1
for i in range(sig_time.shape[1]):
    result_indices = find_sig(first_sig_time_m2[:,i], 5)
    sig_time[result_indices,i] = 1
sig_time = find_last_one_indices(sig_time)
sig_time = np.array(sig_time)
all_sig_end[1] = (sig_time-18)*0.05

#%% 合并O和G的结果
if 'ocean' == monkey:
    all_sig_on_O = all_sig_on
    all_sig_end_O = all_sig_end
if 'grootsw' == monkey:
    all_sig_on_G = all_sig_on
    all_sig_end_G = all_sig_end


#%%
# [spk-e1, spk-e2, spk-m1, spk-m2, theta-e1, theta-e2]
sig_time1 = np.hstack([all_sig_on_O[0], all_sig_on_G[0]])
sig_time2 = np.hstack([all_sig_on_O[1], all_sig_on_G[1]])
sig_time3 = np.hstack([all_sig_on_O[2], all_sig_on_G[2]])
sig_time4 = np.hstack([all_sig_on_O[3], all_sig_on_G[3]])
sig_time5 = np.hstack([all_sig_on_O[4], all_sig_on_G[4]])
sig_time6 = np.hstack([all_sig_on_O[5], all_sig_on_G[5]])

sig_time_end1 = np.hstack([all_sig_end_O[0], all_sig_end_G[0]])
sig_time_end2 = np.hstack([all_sig_end_O[1], all_sig_end_G[1]])


theta_dur1 = [np.arange(sig_time5[i],sig_time_end1[i]+0.01,0.05) for i in range(23)]
theta_dur1 = np.concatenate(theta_dur1)

theta_dur2 = [np.arange(sig_time6[i],sig_time_end2[i]+0.01,0.05) for i in range(23)]
theta_dur2 = np.concatenate(theta_dur2)
#%% plot memory sig time distribution and theta overlap
import seaborn as sns
from scipy.stats import gaussian_kde

with plt.style.context(style_path):
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    fig, axs = plt.subplots(1,1, figsize=(1.4, .9), dpi=300) # , constrained_layout=True
    
    df = pd.DataFrame({'spk_e1':sig_time1, 'spk_m1':sig_time3})
    
    kde = gaussian_kde(theta_dur1, bw_method=0.3)
    x_values = np.linspace(min(theta_dur1), max(theta_dur1), 1000)
    kde_values = kde(x_values)
    kde_values = 3*kde_values/kde_values.max()
    axs.plot(x_values, kde_values, color=colors[4], linewidth=1)
    axs.fill_between(x_values, kde_values, color=colors[4], alpha=0.5)
    
    y = np.zeros([len(df),1])+1
    y_jittered = y + np.random.normal(0, 0.2, size=y.shape)
    x_jittered = np.random.normal(0, 0.007, size=y.shape)
    axs.scatter(df['spk_e1'].values+x_jittered[:,0], y_jittered, facecolor='k',s=20)
    
    y = np.zeros([len(df),1])+2
    y_jittered = y + np.random.normal(0, 0.2, size=y.shape)
    x_jittered = np.random.normal(0, 0.007, size=y.shape)
    axs.scatter(df['spk_m1'].values+x_jittered[:,0], y_jittered,s=20)
    
    rectangle = Rectangle((0, 0), 0.25, 3, linewidth=2, edgecolor='none', facecolor='gray', alpha=0.1)
    axs.add_patch(rectangle)
    axs.set_xlim([-0.1, 0.6])
    axs.set_xticks([0, 0.25, 0.5])
    axs.set_ylim([-0, 3.8])
    axs.set_xlabel('Time from S1 onset')
    axs.set_ylabel('')
    plt.gca().get_yaxis().set_visible(False)
    # plt.gca().get_xaxis().set_visible(False)
    plt.gca().spines['left'].set_visible(False)       # Hides the left spine (y-axis line)
    
    
    fig, axs = plt.subplots(1,1, figsize=(1.4, .9), dpi=300) # , constrained_layout=True
    df = pd.DataFrame({'spk_e2':sig_time2, 'spk_m2':sig_time4})
    
    kde = gaussian_kde(theta_dur2, bw_method=0.3)
    x_values = np.linspace(min(theta_dur2), max(theta_dur2), 1000)
    kde_values = kde(x_values)
    kde_values = 3*kde_values/kde_values.max()
    axs.plot(x_values, kde_values, color=colors[4], linewidth=1)
    axs.fill_between(x_values, kde_values, color=colors[4], alpha=0.5)
    
    y = np.zeros([len(df),1])+1
    y_jittered = y + np.random.normal(0, 0.2, size=y.shape)
    x_jittered = np.random.normal(0, 0.007, size=y.shape)
    axs.scatter(df['spk_e2'].values+x_jittered[:,0], y_jittered, facecolor='k',s=20)
    
    y = np.zeros([len(df),1])+2
    y_jittered = y + np.random.normal(0, 0.2, size=y.shape)
    x_jittered = np.random.normal(0, 0.007, size=y.shape)
    axs.scatter(df['spk_m2'].values+x_jittered[:,0], y_jittered, facecolor=colors[1],s=20)
    
    rectangle = Rectangle((0, 0), 0.25, 3, linewidth=2, edgecolor='none', facecolor='gray', alpha=0.1)
    axs.add_patch(rectangle)
    axs.set_xlim([-0.1, 0.6])
    axs.set_xticks([0, 0.25, 0.5])
    axs.set_ylim([-0, 3.8])
    axs.set_xlabel('Time from S2 onset')
    axs.set_ylabel('')
    plt.gca().get_yaxis().set_visible(False)
    # plt.gca().get_xaxis().set_visible(False)
    plt.gca().spines['left'].set_visible(False)       # Hides the left spine (y-axis line)
#%% plot memory sig time distribution and theta overlap
import seaborn as sns
from scipy.stats import gaussian_kde

with plt.style.context(style_path):
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    fig, axs = plt.subplots(1,1, figsize=(2.5, 1.5), dpi=300) # , constrained_layout=True
    
    df = pd.DataFrame({'theta_e1':sig_time5, 'spk_m1':sig_time3, 'theta_e2':sig_time6, 'spk_m2':sig_time4})
    
    y = np.zeros([len(df),1])+1
    y_jittered = y + np.random.normal(0, 0.2, size=y.shape)
    x_jittered = np.random.normal(0, 0.007, size=y.shape)
    axs.scatter(df['theta_e1'].values+x_jittered[:,0], y_jittered, facecolor=colors[4],s=30)
    
    y = np.zeros([len(df),1])+2
    y_jittered = y + np.random.normal(0, 0.2, size=y.shape)
    x_jittered = np.random.normal(0, 0.007, size=y.shape)
    axs.scatter(df['spk_m1'].values+x_jittered[:,0], y_jittered,s=30)
    
    y = np.zeros([len(df),1])+3
    y_jittered = y + np.random.normal(0, 0.2, size=y.shape)
    x_jittered = np.random.normal(0, 0.007, size=y.shape)
    axs.scatter(df['theta_e2'].values+x_jittered[:,0], y_jittered, facecolor=colors[4],s=30)
    
    y = np.zeros([len(df),1])+4
    y_jittered = y + np.random.normal(0, 0.2, size=y.shape)
    x_jittered = np.random.normal(0, 0.007, size=y.shape)
    axs.scatter(df['spk_m2'].values+x_jittered[:,0], y_jittered,s=30)
    
    # rectangle = Rectangle((0, 0), 0.25, 3, linewidth=2, edgecolor='none', facecolor='gray', alpha=0.1)
    # axs.add_patch(rectangle)
    axs.set_xlim([-0.1, 0.6])
    axs.set_xticks([0, 0.25, 0.5])
    axs.set_yticks([1,2,3,4],['Theta S1', 'Memory S1', 'Theta S2', 'Memory S2'])
    axs.set_xlabel('Time from S1/S2 onset')

stats.ttest_rel(df['theta_e1'], df['spk_m1'])
stats.ttest_rel(df['theta_e2'], df['spk_m2'])

    
    
    
    
    