# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 16:10:25 2020

@author: Littl
"""

import os
import sys
sys.path.append('../')
import psr.cor as cor
import psr.drifting as drifting
import psr.plot as plot

import matplotlib.pyplot as plt
import numpy as np

folder_name = ['J0837-4135', 'J0942-5552',
               'J1327-6222','J1401-6357','J1644-4559',
                           'J1752-2806',
               'J1901-0906','J1456-6843',
               'J1651-4246','J0738-4042']#,
folder_name = ['J0738-4042']
pulse_end = { 'J0738-4042':[590,750],'J0837-4135':[655,735], 'J0942-5552':[330,480],
             'J1327-6222':[65,130],'J1644-4559':[750,840],'J1401-6357':[80,135],
                                  'J1752-2806':[710,765],
            'J1901-0906':[710,880],'J1456-6843':[520,720],
            'J1651-4246':[400,750],'J1825-0935':[170,750]
        }

mannual =  { 
        'J0738-4042':[2038, 2318, 2496, 2865, 2947, 3154, 3514, 5804, 5893, 7736, 8964,
        8965, 8966, 8967, 9300, 9301, 9304, 9306, 9307],
        'J0837-4135':list(range(2400,2500)),
        'J0942-5552':[204,206,800,801,802,803,1201,1550,1551,1628,1629,1634,1635,1663,4789],
        'J1327-6222':[1866,2327,2328],'J1401-6357':[],'J1644-4559':[],
        'J1651-4246':[],'J1752-2806':[],'J1825-0935':[2330,2331,2333,2335,2447,2448,3284],
        'J1901-0906':[],'J1456-6843':[]
        
        }


lags = [0,1,2,3,4,5]
for psr_name in folder_name:
    
    file_list = os.listdir(r'.\\'+psr_name+'_v2')
    
    os.chdir('./revised_v3/cor')
    
    for filename in file_list:
        #filename = 'J1456-6843_704-1727MHz.npy'
        data = np.load('..\\..\\%s_v2\\%s'%(psr_name,filename))
        
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
        #data[distb] = np.nan
        #data[mannual[psr_name]] = np.nan
        
        freq = filename.split('.')[0].split('_')[-1]
        freq_chn = []
        freq_chn.append(int(freq.split('-')[0]))
        freq_chn.append(int(freq.split('-')[-1].split('M')[0]))
        for c,channel in enumerate(freq_chn):
            if channel%2 == 1:
                freq_chn[c] += 1
        if freq_chn[0] > freq_chn[1]:
            freq_chn.reverse()
        freq = '%d - %d MHz'%(freq_chn[0],freq_chn[1])
        
        plot.plot_all(w,data,peak_end,pulse_end,lags,psr_name,freq)

        
    os.chdir('../../')
    
