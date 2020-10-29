# -*- coding: utf-8 -*-
"""
Created on Mon Jun 15 11:52:38 2020

@author: Littl
"""

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
#folder_name = ['J1825-0935']
pulse_end = { 'J0738-4042':[350,650],'J0837-4135':[650,730], 'J0942-5552':[330,480],
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


lags = [0,1,2,3,4,5]

    
file_list = os.listdir('./0837')
psr_name = 'J0837-4135'

for filename in file_list[::-1]:
    #filename = 'J1456-6843_704-1727MHz.npy'
    data = np.load('./0837/%s'%(filename))
    
    peak_end = pulse_end[psr_name]
    left=peak_end[0]
    right=peak_end[1]        

        
    #base_data = np.concatenate((data[:,:,:peak_end[0]],data[:,:,peak_end[1]:]),axis=-1)
    #base_data = np.nanmean(base_data,axis=-1)
    #data = data - base_data[:,:,None]
   
    w = np.ones((data.shape[0],))

    
    
    freq = filename.split('.')[0].split('_')[-1]

    
    plot.plot_all(w,data,peak_end,pulse_end,lags,psr_name,freq)

