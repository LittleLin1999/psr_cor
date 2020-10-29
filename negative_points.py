# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 16:22:34 2020

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

import pandas as pd

folder_name = ['J0837-4135', 'J0738-4042',
               'J1327-6222','J1644-4559',
                           'J1752-2806']


pulse_end = { 'J0738-4042':[590,800],'J0837-4135':[655,730], 'J0942-5552':[330,480],
             'J1327-6222':[65,130],'J1644-4559':[750,840],'J1401-6357':[80,135],
                                  'J1752-2806':[710,765],
            'J1901-0906':[710,880],'J1456-6843':[520,720],
            'J1651-4246':[400,750],'J1825-0935':[170,750]
        }

mannual =  { 
        'J0738-4042':[],'J0837-4135':[2339],
        'J0942-5552':[204,206,800,801,802,803,1201,1550,1551,1628,1629,1634,1635,1663,4789],
        'J1327-6222':[1866,2327,2328],'J1401-6357':[],'J1644-4559':[],
        'J1651-4246':[],'J1752-2806':[],'J1825-0935':[2330,2331,2333,2335,2447,2448,3284],
        'J1901-0906':[],'J1456-6843':[]
        
        }

info = []

for psr_name in folder_name:
    
    
    file_list = os.listdir(r'.\\'+psr_name+'_v2')
    if psr_name == 'J0837-4135':
        file_list.remove('J0837-4135_704-1343MHz.npy')
        
    for filename in file_list:
        print(filename)
        
        data = np.load('.\\%s_v2\\%s'%(psr_name,filename))
        
        peak_end = pulse_end[psr_name]
        left=peak_end[0]
        right=peak_end[1]        

            
        base_data = ((((data[:,:,:peak_end[0]])).sum(-1) + ((data[:,:,peak_end[1]:])).sum(-1))/
                        (data[:,:,:peak_end[0]].shape[-1]+data[:,:,peak_end[1]:].shape[-1]))
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
    
        
    
        I_prof = data[:,0,peak_end[0]:peak_end[1]].mean(0)#(w[:,None]*data_all[:,0,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
        Q_peak = data[:,1,peak_end[0]:peak_end[1]].mean(0)#(w[:,None]*data[:,1,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
        U_peak = data[:,2,peak_end[0]:peak_end[1]].mean(0)#(w[:,None]*data[:,2,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
        L_prof = np.sqrt(Q_peak**2+U_peak**2)
        V_prof = (w[:,None]*data[:,3,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
        
        top = max(I_prof)
        I_prof = I_prof/top
        L_prof = L_prof/top
        V_prof = V_prof/top
        
        
        Q = data[:,1]
        V = data[:,2]
        L = np.sqrt(Q**2+V**2)
        base_L = (L[:,:peak_end[0]].sum(1)+L[:,peak_end[1]:].sum(1))/(L[:,:peak_end[0]].shape[1]+L[:,peak_end[1]:].shape[1])
        L = L- base_L[:,None]
        L = L/top
        
        I = data[:,0]/top
        abs_V = abs(data[:,-1]/top)
        
        
        phase_peak = phase[peak_end[0]:peak_end[1]]

        intensity = {'I':I,'L':L,'V':abs_V}
        
        for key,d in intensity.items():
            
            pulse_off1 = d[:,peak_end[0]-50:peak_end[0]]
            pulse_off2 = d[:,peak_end[1]:peak_end[1]+50]
            pulse_off = np.concatenate((pulse_off1,pulse_off2),axis=1)
            
            cr_off = cor.weighed_corr_map(w, pulse_off, 0,-10000)
            cr_off = np.concatenate([cr_off[:int(cr_off.shape[0]/3),-int(cr_off.shape[1]/3):],
                                        cr_off[-int(cr_off.shape[0]/3):,:int(cr_off.shape[1]/3)]])
            cr_off = np.std(cr_off)
        
            data_peak = d[:,peak_end[0]:peak_end[1]]
            cr = cor.weighed_corr_map(w, data_peak,0,abs(cr_off))
            
            norm=mpl.colors.Normalize(vmin=-1,vmax=1) #将1设为最红，-1为最蓝，中间的0为白
            cmap = mpl.cm.seismic
            cmap.set_bad(color='#404040',alpha=0.15)
            corr_map = plt.imshow(cr,norm=norm, cmap=cmap, aspect='equal', origin='lower',
                                  extent=(phase_peak[0], phase_peak[-1], phase_peak[0], phase_peak[-1]))
            plt.show()
            
            vmin = np.nanmin(cr)
            
            pos_min = np.where(cr==vmin)
            pos_min = pos_min[0][0],pos_min[1][0]
            phase_min = phase_peak[pos_min[0]],phase_peak[pos_min[1]]
            dist = abs(phase_min[0]-phase_min[1])
            print(phase_min,vmin)
            
            info.append([psr_name,freq,key,phase_min,vmin,dist])
            
xml = pd.DataFrame(np.array(info),columns=('PSR','FREQ','INTENSITY','PHASE','VMIN','DIST'))
xml.to_csv('./revised_v2/negative_points.csv')
            