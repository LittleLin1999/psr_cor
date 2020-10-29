# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 14:29:10 2020

@author: Littl
"""

import os
import sys
sys.path.append('../')
import psr.cor as cor
import psr.drifting as drifting
import psr.plot as plot

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

folder_name = ['J1825-0935']
pulse_end = { 'J1825-0935':[[165,235],[610,750]]}

mannual =  { 'J1825-0935':[2330,2331,2333,2335,2447,2448,3284]}


def group_consecutive(a):
    return np.split(a, np.where(np.diff(a) != 1)[0] + 1)


for psr_name in folder_name:
    
    
    file_list = os.listdir(r'.\\'+psr_name+'_v2')
    if psr_name == 'J0837-4135':
        file_list.remove('J0837-4135_704-1343MHz.npy')
        
    for filename in file_list:
        print(filename)
        
        data = np.load('.\\%s_v2\\%s'%(psr_name,filename))
        
        peak_end = pulse_end[psr_name]
    
        pulse_off = np.concatenate([data[:,:,:peak_end[0][0]],
                            data[:,:,peak_end[0][1]:peak_end[1][0]],
                            data[:,:,peak_end[1][1]:]],
                            axis=-1)
        base_data = np.nanmean(pulse_off,axis=-1)
        data = data - base_data[:,:,None]
       
        w = np.ones((data.shape[0],))
        mean = data[:,0,:].mean(-1)
        distb = np.where(mean<=0)
        w[distb] = 0
        w[mannual[psr_name]] = 0
        
        
        freq = filename.split('.')[0].split('_')[-1]
        freq_chn = []
        freq_chn.append(int(freq.split('-')[0]))
        freq_chn.append(int(freq.split('-')[-1].split('M')[0]))
        for c,channel in enumerate(freq_chn):
            if channel%2 == 1:
                freq_chn[c] += 1
        freq = '%d-%d MHz'%(freq_chn[0],freq_chn[1])
        
        phase = np.linspace(0,360,1025)[:1024]
        
        
    
        I_prof = data[:,0,:].mean(0)#(w[:,None]*data_all[:,0,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
    
        top = max(I_prof)
 
        
        Q = data[:,1]
        V = data[:,2]
        L = np.sqrt(Q**2+V**2)
        L_base = np.nanmean(np.concatenate([L[:,:peak_end[0][0]],L[:,peak_end[1][1]:],
                         L[:,peak_end[0][1]:peak_end[1][0]]],axis=-1))
        L = L- L_base
        L = L/top
        
        I = data[:,0]/top
        abs_V = abs(data[:,-1]/top)
        
        
        p_end = peak_end

        intensity = {'I':I,'L':L,'V':abs_V}
        
        for key,d in intensity.items():
            
            pulse_off11 = d[:,p_end[0][0]-50:p_end[0][0]]
            pulse_off12 = d[:,p_end[0][1]:p_end[0][1]+50]
            pulse_off21 = d[:,p_end[1][0]-50:p_end[1][0]]
            pulse_off22 = d[:,p_end[1][1]:p_end[1][1]+50]
            
            pulse_off = np.concatenate((pulse_off11,pulse_off12,pulse_off21,pulse_off22),
                                       axis=1)
            
            cr_off = cor.weighed_corr_map(w, pulse_off, 0,-10000)
            cr_off = np.concatenate([cr_off[:int(cr_off.shape[0]/3),-int(cr_off.shape[1]/3):],
                                        cr_off[-int(cr_off.shape[0]/3):,:int(cr_off.shape[1]/3)]])
            cr_off = np.std(cr_off)
        
            data_peak = np.concatenate([d[:,p_end[0][0]:p_end[0][1]],
                                        d[:,p_end[1][0]:p_end[1][1]]],axis=1)
            cr = cor.weighed_corr_map(w, data_peak,0,abs(cr_off))
            
             
            norm=mpl.colors.Normalize(vmin=-1,vmax=1) #将1设为最红，-1为最蓝，中间的0为白
            cmap = mpl.cm.seismic
            cmap.set_bad(color='#404040',alpha=0.15)
            corr_map = plt.imshow(cr,norm=norm, cmap=cmap, aspect='equal', origin='lower')
                                  
            plt.show()
            
            tril = np.tril(cr, k=0)
            diag_part = []
            flag = np.full(cr.shape,0.)
            for j in np.arange(len(cr)):
                sequence = np.where(tril[j]>0)[0]
                inner = group_consecutive(sequence)[-1]
                if len(inner)>3:
                    diag_part.append(list(inner))
                    flag[j][inner[0]:inner[-1]] = 1
                
            diag_part = np.array(diag_part)
            
            
            wnum = np.nanmean(flag,axis=1)
            wnum[np.where(wnum<=0)] = np.nan
            wnum[np.where(wnum>0)] = 1
            wnum = np.nansum(wnum)
            
            flag[np.where(flag==0)] = np.nan
            wsum = np.nansum(flag)
            
            
            length = wsum/wnum
            
            phase_length = length/abs(peak_end[0][0]-peak_end[0][-1])*abs(phase[peak_end[0][0]]-phase[peak_end[0][-1]])
            
            print(wsum,wnum,'length=%f'%phase_length)
            corr_map = plt.imshow(flag,norm=norm, cmap=cmap, aspect='equal', origin='lower')
                
            plt.show()
            
            info.append([psr_name,freq,key,phase_length])