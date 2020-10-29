# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 14:33:58 2020

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

folder_name = ['J0837-4135', 'J0942-5552',
               'J1327-6222','J1401-6357','J1644-4559',
                           'J1752-2806',
               'J1456-6843',
               'J1651-4246','J0738-4042']#,
#folder_name = ['J1327-6222']
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
for psr_name in folder_name:
    print(psr_name)
    file_list = os.listdir(r'../%s_v2'%psr_name)
    file_list = os.listdir(r'..\\'+psr_name+'_v2')
    
    #file_list[0], file_list[1], file_list[-1] = file_list[-1], file_list[0], file_list[1]
    file_list = file_list[::-1]
    
    phase = np.linspace(0,360,1024)
    peak_end = pulse_end[psr_name]
    phase_peak = phase[peak_end[0]:peak_end[1]]
    
   
    
    fig = plt.figure(figsize=(6,4*len(file_list)))

    grid=plt.GridSpec(len(file_list),1,fig,wspace=0.,hspace=0.)
    ax = []
    for j in range(len(file_list)):
        axis = plt.subplot(grid[j,0])
        ax.append(axis)
        
    for index,filename in enumerate(file_list):
        channel = filename.split('_')[-1].split('.')[0]
        channel = channel.split('-')[0],channel.split('-')[1].split('M')[0]
        channel = list(channel)
        for c,chn in enumerate(channel):
            chn = int(chn)
            if chn%2==1:
                chn+=1
            channel[c] = chn
        if channel[0] < 1000:
            sub = 0
        if channel[0] > 1000 and channel[0] < 2000:
            sub = 1
        if channel[0] > 2500:
            sub = 2
        data = np.load('..\%s_v2\\%s'%(psr_name,filename))
        pulse_off = np.concatenate([data[:,0,:peak_end[0]],
                            data[:,0,peak_end[1]:]],
                            axis=-1)
        base_data = np.nanmean(pulse_off,axis=-1)[:,None]
        data = data - base_data[:,:,None]
       
        w = np.ones((data.shape[0],))
        mean = data[:,0,:].mean(-1)
        distb = np.where(mean<=0)
        w[distb] = 0
        w[mannual[psr_name]] = 0
        
        data[np.where(w==0)] = np.nan
        
        u = np.nanmean(data[:,1,:],axis=0)#(w[:,None]*data[:,1,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
        q = np.nanmean(data[:,2,:],axis=0)#(w[:,None]*data[:,2,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
        
        L = np.sqrt(u**2+q**2)
        L_base = np.nanmean(np.concatenate([L[peak_end[1]:],L[:peak_end[0]]],
   
                                 axis=-1))
        L = L-L_base
        stdL = np.nanstd(np.concatenate([L[peak_end[1]:],L[:peak_end[0]]],
                                 axis=-1))
        snrL = np.abs(L/stdL)
        
        pu_err = np.nanstd(np.concatenate([u[peak_end[1]:],u[:peak_end[0]]],
                                 axis=-1))
        pq_err = np.nanstd(np.concatenate([q[peak_end[1]:],q[:peak_end[0]]],
                                 axis=-1))

        pu = u[peak_end[0]:peak_end[1]]
        pq = q[peak_end[0]:peak_end[1]]
        
        ang=np.arctan2(pu, pq)/np.pi*180./2.
        ang[ang<0] += 180
        
        ang_err2 = ((pu*pq_err)**2 + (pq*pu_err)**2)/(pq**2 + pu**2)**2
        ang_err = np.sqrt(ang_err2)/np.pi*180/2
        
        idx = np.where(snrL[peak_end[0]:peak_end[1]]>3)
    
        ax[sub].errorbar(phase_peak[idx],ang[idx],ang_err[idx],c='k',fmt='.',linewidth=0.2,ms=1.5)
        ax[sub].set_ylim(max(min(ang)-5,-180),min(max(ang)+5,180))
        ax[sub].locator_params(axis='y',nbins=4)
        ax[sub].tick_params(labelbottom=False)
        
        ax[sub].set_xlim(phase_peak[0],phase_peak[-1])
        
        if psr_name in ['J0942-5552','J0837-4135']:
            ax[sub].text(0.02,0.02,'%d - %d MHz'%(channel[0],channel[1]),
                          verticalalignment='bottom', 
                          horizontalalignment='left',transform=ax[sub].transAxes)
        else:
            ax[sub].text(0.02,0.97,'%d - %d MHz'%(channel[0],channel[1]),
                          verticalalignment='top', 
                          horizontalalignment='left',transform=ax[sub].transAxes)
        ax[-1].tick_params(labelbottom=True)
        
        plt.gca().xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
        #ax[-1][side].set_xlabel('phase/$^\circ$')
            
        
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.ylabel('PA/$^\circ$')     
    plt.xlabel('phase/$^\circ$')
    plt.savefig('%s_PA.pdf'%psr_name)
    plt.show()