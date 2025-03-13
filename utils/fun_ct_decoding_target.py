# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 17:24:48 2023

@author: Wen
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn import svm
from sklearn.model_selection import StratifiedKFold, RepeatedStratifiedKFold, RepeatedKFold, KFold,LeaveOneOut
from sklearn import preprocessing
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.metrics import confusion_matrix, balanced_accuracy_score
from sklearn.inspection import permutation_importance
from sklearn.impute import SimpleImputer

style_path = r'C:\Wen\OneDrive\CodeHub\Python\style_paper.mplstyle'
#%%
def fill_missing_val(tf): ## fill missing values
    imp = SimpleImputer(missing_values=np.nan, strategy='mean')
    for i in range(tf.shape[2]):
        x = tf[:,:,i]
        tf[:,:,i] = imp.fit_transform(x)
    return tf
    
def norm_stackTime(x):
    n_sample = x.shape[0]
    n_feature = x.shape[1]
    n_time = x.shape[2]
    
    ### preprocessing ##
    scaler = StandardScaler()
    x = np.transpose(x, [1, 2, 0])
    x_norm = x.reshape(n_feature, n_time*n_sample, order='F')
    x_norm = scaler.fit_transform(x_norm.T).T
    x_norm = x_norm.reshape(n_feature, n_time, n_sample, order='F')
    x = np.transpose(x_norm, [2, 1, 0])
    return x

def ct_decoding(x, y, train_r, test_r, shuffle, n_repeats):
    n_feature = x.shape[0]
    n_sample = x.shape[1]
    n_time = x.shape[2]
    
    scaler = StandardScaler()
    print(x.shape, ['shuffle:', shuffle])
    
    shape_3d = x.shape
    xx = x.reshape(shape_3d[0], -1)
    xx = scaler.fit_transform(xx.T).T
    x = xx.reshape(shape_3d)

    #### svc params ####
    C = 1
    gamma = 1
    n_fold = 3 
    n_repeats = n_repeats
    if shuffle==1:
        n_repeats = n_repeats
    print([n_fold, n_repeats])

    clf = svm.SVC(kernel='linear', C=C, gamma=gamma)
    scores =  np.full((n_repeats, n_fold, n_time, n_time), np.nan)
    fea_weight = np.full((n_repeats, n_fold, x.shape[0], n_time), np.nan)
    for r in range(n_repeats):
        if np.mod(r,10) == 0:
            print(r)
        if shuffle == 1:
            np.random.shuffle(y[:,0])
            np.random.shuffle(y[:,1])
            
        rkf = StratifiedKFold(n_splits=n_fold, shuffle=True, random_state=None)  # for within rank
        for i, (train, test) in enumerate(rkf.split(x[:,:,0].T, y[:,train_r])):
            # print(train)
            y_train = y[train, train_r]
            for t in range(n_time): 
                x_train = x[:,train,t].T
                #### cross rank
                # if t < mkts_merge[0][1]:
                #     y_train = y[train, 0]
                # else:
                #     y_train = y[train, 1]
                    
                clf.fit(x_train, y_train)
                fea_weight[r,i,:,t] = clf.coef_.mean(0)
                for tt in range(n_time):
                    x_test = x[:,test,tt].T
                    y_test = y[test, test_r]
                    #### cross rank
                    # if tt < mkts_merge[0][1]:
                    #     y_test = y[test, 0]
                    # else:
                    #     y_test = y[test, 1]
                        
                    scores[r,i,t,tt] = clf.score(x_test, y_test)
    return scores, fea_weight

def ct_decoding_cr_rank(x, y, train_r, test_r, shuffle, n_repeats):
    n_feature = x.shape[0]
    n_sample = x.shape[1]
    n_time = x.shape[2]
    
    scaler = StandardScaler()
    print(x.shape, ['shuffle:', shuffle])
    shape_3d = x.shape
    xx = x.reshape(shape_3d[0], -1)
    xx = scaler.fit_transform(xx.T).T
    x = xx.reshape(shape_3d)

    #### svc params ####
    C = 1
    gamma = 1
    n_fold = 3 
    n_repeats = n_repeats
    print([n_fold, n_repeats])

    clf = svm.SVC(kernel='linear', C=C, gamma=gamma)
    scores =  np.full((n_repeats, n_fold, n_time, n_time), np.nan)
    fea_weight = np.full((n_repeats, n_fold, 15, x.shape[0], n_time), np.nan)
    for r in range(n_repeats):
        if np.mod(r,10) == 0:
            print(r)
        if shuffle == 1:
            np.random.shuffle(y[:,0])
            np.random.shuffle(y[:,1])
            
        rkf = StratifiedKFold(n_splits=n_fold, shuffle=True, random_state=None)  # for within rank
        for i, (train, test) in enumerate(rkf.split(x[:,:,0].T, y[:,train_r])):
            # print(train)
            y_train = y[train, train_r]
            for t in range(n_time): 
                x_train = x[:,train,t].T
                #### cross rank
                # if t < mkts_merge[0][1]:
                #     y_train = y[train, 0]
                # else:
                #     y_train = y[train, 1]
                    
                clf.fit(x_train, y_train)
                # fea_weight[r,i,:,:,t] = clf.coef_
                for tt in range(n_time):
                    x_test = x[:,test,tt].T
                    y_test = y[test, test_r]
                    #### cross rank
                    # if tt < mkts_merge[0][1]:
                    #     y_test = y[test, 0]
                    # else:
                    #     y_test = y[test, 1]
                        
                    scores[r,i,t,tt] = clf.score(x_test, y_test)
    return scores, fea_weight


def plot_acc_mat(scores_mean, rank, mkts_merge, figName):
    with plt.style.context(style_path):
        fig = plt.figure(figsize = (1.8, 1.6), dpi = 300)
        ax1 = fig.add_subplot(1, 1, 1)
        vmin = 0.15
        vmax = 0.55
        h1 = ax1.pcolormesh(scores_mean, cmap = plt.cm.RdBu_r, vmin=vmin, vmax=vmax)
        ax1.vlines(mkts_merge[0,:], 0, mkts_merge[1][-1], color='w', linestyle = '-')
        ax1.hlines(mkts_merge[0], 0, mkts_merge[1][-1], color='w', linestyle = '-')

        ax1.set_xlabel('Test time')
        ax1.set_ylabel('Train time')
        fig.colorbar(h1)
        plt.show()


