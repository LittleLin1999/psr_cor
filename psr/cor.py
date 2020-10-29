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
import scipy.signal as sig
import mpmath as mm
import scipy.special as sc
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
def sigma(n,rho):
    
    inner0 = n/2*mm.log(1-rho**2)
    
    inner1 = mm.log(mm.hyp3f2(3/2,n/2,n/2,1/2,n/2+1,rho**2))
  
    inner2 = -mm.log(n)
    inner = mm.mp.e**(inner0+inner1+inner2)
    
    outer0 = n/2*np.log(1-rho**2)+np.log(abs(rho))
    outer1 = 2*sc.gammaln(n/2+1/2)
    outer2 = mm.log(mm.hyp2f1(n/2+1/2,n/2+1/2,n/2+1,rho**2))
    outer3 = -sc.gammaln(n/2)
    outer4 = -sc.gammaln(n/2+1)
    outer = outer0+outer1+outer2+outer3+outer4
    outer = 2*outer
    
    return mm.sqrt(inner - mm.mp.e**(outer))

def weighed_corr(w1,w2,l1,l2,lag):
    
    n = len(l1)
    

    s1 = l1[lag:]
    s2 = l2[:n-lag]
    
    ww1 = w1[lag:]
    ww2 = w2[:n-lag]
    
    w = ww1*ww2
    
    y1 = w*s1
    y2 = w*s2
    
    avg1 = np.nansum(y1)/np.nansum(w)
    avg2 = np.nansum(y2)/np.nansum(w)

    cov = np.nansum((y1-avg1)*(y2-avg2))
    
    var1 = np.nansum((y1-avg1)**2)
    var2 = np.nansum((y2-avg2)**2)
    
   
    coeff = cov/np.sqrt(var1*var2)
    
    '''
    N = len(w[w!=0])
    #print(N)
    sg = np.nan
    if coeff < 0.6:
        sg = float(sigma(N,coeff))
    '''

    return coeff

    
def weighed_corr_map(w,l,lags,sigma):

    ll = l.T
    
    cr = np.zeros((len(ll),len(ll)))
    #delta = np.zeros((len(ll),len(ll)))
    for j in range(len(ll)):
        
        for k in range(len(ll)):
            
            coe = weighed_corr(w,w,ll[j], ll[k], lags)
            
            if abs(coe) >= 3*sigma:
                cr[j][k] = coe
            else:
                cr[j][k] = np.nan
            #delta[j][k]=delta_c
    return cr#,delta



def weighed_corr_map_freq(w1,w2,l1,l2,lags):
    
    l1 = l1.T
    l2 = l2.T
    cr = np.zeros((len(l2),len(l1)))

    for j in range(len(l2)):
        for k in range(len(l1)):
            coe = weighed_corr(w1,w2,l1[k], l2[j], lags)
            cr[j][k] = coe
  
    return cr


#互相关
def weighed_corr_noi(w1, w2, l1, l2, lag, noise1, noise2, flagnoi): #!!!两个noise对应两个相位处的噪声序列,flag决定要不要扣噪声
    n = len(l1)

    s1 = l1[lag:]
    s2 = l2[:n - lag]

    ww1 = w1[lag:]
    ww2 = w2[:n - lag]
    w = ww1 * ww2

    y1 = w * s1
    y2 = w * s2

    avg1 = np.nansum(y1) / np.nansum(w)
    avg2 = np.nansum(y2) / np.nansum(w)

    cov = np.nansum((y1 - avg1) * (y2 - avg2))

    var1 = np.nansum((y1 - avg1) ** 2)
    var2 = np.nansum((y2 - avg2) ** 2)

    if flagnoi == 1: #!!!
        noise1 = noise1[lag:] * w #!!!
        noise2 = noise2[:n - lag] * w #!!!
        var01 = var1 - np.sum(noise1) * (1 - 1 / np.sum(w)) #!!!
        var02 = var2 - np.sum(noise2) * (1 - 1 / np.sum(w)) #!!!
        if var01 > 0 and var02 > 0: #!!! 万一出现根号里面小于0,就强行不扣噪声
            #print('woop！')
            var1 = var01 #!!!
            var2 = var02 #!!!

    coeff = cov / np.sqrt(var1 * var2)

    N = len(w[w!=0])
    P1 = (y1-avg1)**2
    P2 = (y2-avg1)**2
    par1 = np.sum(P1*P2)
    par2 = var1*var2
    delta_c = (1-coeff**2)/N/np.sqrt(par2)*np.sqrt(par1)
    return coeff,delta_c

def weighed_corr_map_noi(w1, w2, l1, l2, lag, noise1, noise2, noil, noir, flagnoi): #!!! noil和noir是扣噪声区域的左右边界,单位是bin,并且以脉冲起始位置left作为0
    ll1 = l1.T
    ll2 = l2.T
    #cr,n分别记录全部相关图/负相关部分的相关图
    cr = np.zeros((len(ll1), len(ll2)))
    delta = np.zeros((len(ll1), len(ll2)))
    for j in range(len(ll1)):
        for k in range(len(ll2)):
            if lag == 0 and j == k: #!!!
                flagnoi = 0 #!!!
            elif j < noil or j > noir or k < noil or k > noir: #!!!
                flagnoi = 0 #!!!
            else: #!!!
                flagnoi = 1 #!!!
            coe,delta_c = weighed_corr_noi(w1, w2, ll1[j], ll2[k], lag, noise1, noise2, flagnoi) #!!!
            if abs(coe) >1:
                coe = np.nan
            cr[j][k] = coe
            delta[j][k] = delta_c
    return cr, delta



    
    