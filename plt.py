# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 21:09:43 2020

@author: Littl
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 10:10:53 2020

@author: Littl
"""

import sys
sys.path.append('../')
import psr.cor as cor
import psr.drifting as drifting
import psr.plot as plot

import matplotlib.pyplot as plt
import numpy as np

star_name = 'J1401-6357'
peak_end = [50,150]
pulse_end = [90,130]
lags = [1]
mannual = []

import os
filelist = os.listdir(r'C:\Users\Littl\Desktop\psrchive\psrchive1401')[:-1][:]


#%% 
for filename in filelist:
    
    data = np.load(filename,allow_pickle=True)
    base_data = ((data[:,:,:peak_end[0]].sum(-1) + data[:,:,peak_end[1]:].sum(-1))/
                    (data[:,:,:peak_end[0]].shape[-1]+data[:,:,peak_end[1]:].shape[-1]))
    data = data - base_data[:,:,None]
    left=peak_end[0]
    right=peak_end[1]

    
    w = np.ones((data.shape[0],))
    mean = data[:,0].mean(1)
    distb = np.where(mean<=0)
    w[distb] = 0
    w[mannual] = 0

    freq = filename.split('.')[0].split('_')[-1]
  #  plot.plot_all(w,data,peak_end,pulse_end,lags,star_name,freq)


    for i in lags:
        #plot_cor_noise(True,w,w,data,data,peak_end,pulse_end,star_name,freq,i,bin_num=1024)
        plot.plot_cor(w,data[:,0],'I',peak_end,i,star_name,freq,bin_num=1024)
        plot.plot_cor_noise(True,w,w,data[:,0],data[:,0],'I',peak_end,pulse_end,star_name,freq,1,bin_num=1024)
       # Q = data[:,1]
       # V = data[:,2]
        #L = np.sqrt(Q**2+V**2)
        #base_L = (L[:,:peak_end[0]].sum(1)+L[:,peak_end[1]:].sum(1))/(L[:,:peak_end[0]].shape[1]+L[:,peak_end[1]:].shape[1])
      #  L = L- base_L[:,None]
        #plot.plot_cor(w,L,'L',peak_end,i,star_name,freq,bin_num=1024)
        #plot.plot_cor(w,data[:,-1],'V',peak_end,i,star_name,freq,bin_num=1024)

    
#%%

#plot.plot_cor_freq(w1,w2,data1,data2,[760,830],star_name,freq1,freq2,lag=0,bin_num=1024)
#plot.plot_cor_freq(w1,w3,data1,data3,[760,830],star_name,freq1,freq3,lag=0,bin_num=1024)
#plot.plot_cor_freq(w2,w3,data2,data3,[760,830],star_name,freq2,freq3,lag=0,bin_num=1024)    