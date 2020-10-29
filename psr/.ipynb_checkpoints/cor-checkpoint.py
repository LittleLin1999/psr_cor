# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 15:15:14 2020

@author: Little Lin
"""

"""
去除时域干扰：
1. imshow : 能否清晰显示子脉冲图像
2. corr ： 是否存在干扰线
3. np.argmax(pulse) : 找到整个序列中最值点的位置（注意这个结果给的是flatten的index）
4. plot或者subplot

干扰周期的weight设为0
"""
import numpy as np
import matplotlib.pyplot as plt

#%% correlation function
'''
普通的互相关
'''
#phase difference
def serial_corr(wave1,wave2, lag=1):
    n = len(wave1) 
    wave1=wave1[: n]
    wave2=wave2[: n]
    y1 = wave1[lag:]
    y2 = wave2[:n-lag]
    corr = np.corrcoef(y1, y2, ddof=0)[0, 1]
    return corr

#auto-correlation
def autocorr(wave):
    lags = range(len(wave))
    corrs = [serial_corr(wave,wave,lag) for lag in lags]
    return lags, corrs

#cross-correlation
def crosscorr(wave1,wave2):
    lags = range(len(wave1))
    corrs = [serial_corr(wave1,wave2,lag) for lag in lags]
    return lags, corrs

def cor_map(ll,lags):
    cr = np.zeros((len(ll),len(ll)))
    n = np.zeros((len(ll),len(ll)))
    for j in range(len(ll)):
        for k in range(len(ll)):
            coe=serial_corr(ll[j], ll[k], lags)
            cr[j][k]=coe
            if coe<0:
                n[j][k]=abs(coe)
    return cr,n



#%%
'''
进阶版互相关
'''
def weighed_corr(w1,w2,l1,l2,lag):
    
    n = len(l1)
    
    s1 = l1[lag:]
    s2 = l2[:n-lag]
    
    ww1 = w1[lag:]
    ww2 = w2[:n-lag]
    w = ww1*ww2
    
    y1 = w*s1
    y2 = w*s2
    
    avg1 = np.sum(y1)/np.sum(w)
    avg2 = np.sum(y2)/np.sum(w)
    
    cov = np.sum((y1-avg1)*(y2-avg2))
    
    var1 = np.sum((y1-avg1)**2)
    var2 = np.sum((y2-avg2)**2)
    
    coeff = cov/np.sqrt(var1*var2)
    
    return coeff

    
def weighed_corr_map(w,l,lags):

    ll = l.T
    cr = np.zeros((len(ll),len(ll)))
    n = np.zeros((len(ll),len(ll)))
    for j in range(len(ll)):
        for k in range(len(ll)):
            coe = weighed_corr(w,w,ll[j], ll[k], lags)
            cr[j][k] = coe
            if coe<0:
                n[j][k]=abs(coe)
    return cr,n





