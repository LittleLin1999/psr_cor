# -*- coding: utf-8 -*-
"""
Analysis for subpulse drifting

@author: Littl
"""
import numpy as np

#%%
#LRFS
def LRFS(peak):
    
    LRFS = []
    phase = []
    
    for i in range(peak.shape[1]):
        fft = np.fft.rfft(peak[:,i]-peak[:,i].mean(),axis=-1)
        amp = abs(fft)
        i_fp = np.argmax(amp)
        fp = fft[i_fp]
        p = np.rad2deg(np.arctan2(np.imag(fp),np.real(fp)))  #phase
        LRFS.append(amp) #amplitude
        phase.append(p)
    LRFS = np.array(LRFS).T
    phase = np.array(phase)
    
    freq = np.fft.rfftfreq(peak.shape[0])
    return LRFS,phase,freq
'''    
data = np.load("C:\\Users\\Littl\\Documents\\GitHub\\pulsar\\psr1644\\sub_high_I.npy")
LRFS(data[:,770:820]) 
'''   
#Time variation of LRFS

def sliding_LRFS(peak,lenth):
    S = []
    for i in np.arange(0,peak.shape[0]-lenth,step=10):
        block = peak[i:i+lenth,:]  
        #LRFS of the chosen block
        fft = abs(np.fft.rfft(block-block.mean(0)[None,:],axis=0))
        S.append(fft.mean(1))
    S = np.array(S).T
    freq=np.fft.rfftfreq(lenth)
    
    return S,freq