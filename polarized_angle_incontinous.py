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

folder_name = ['J1825-0935']
text_pos = {'J1825-0935':[0.05,0.95,'top','left']}
pulse_end = { 'J1825-0935':[[165,235],[610,750]]
        }

mannual =  { 'J1825-0935':[2330,2331,2333,2335,2447,2448,3284]
        }
for psr_name in folder_name:
    print(psr_name)
    file_list = os.listdir(r'./%s_v2'%psr_name)
    file_list = os.listdir(r'.\\'+psr_name+'_v2')
    
    #file_list[0], file_list[1], file_list[-1] = file_list[-1], file_list[0], file_list[1]
    file_list = file_list[::-1]
    
    
    peak_end = pulse_end[psr_name]
    
    phase = np.linspace(0,360,1024)
    
    fig = plt.figure(figsize=(6,4*len(file_list)))

    grid=plt.GridSpec(len(file_list),3,fig,wspace=0.1,hspace=0.)
    ax = []
    for j in range(len(file_list)):
        axis = plt.subplot(grid[j,0]),plt.subplot(grid[j,1:])
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
            if psr_name == 'J1644-4559':
                peak_end = [750,1023]
                phase_peak = phase[peak_end[0]:peak_end[1]]
            sub = 0
        if channel[0] > 1000 and channel[0] < 2000:
            sub = 1
        if channel[0] > 2500:
            sub = 2
        data = np.load('.\%s_v2\\%s'%(psr_name,filename))
        pulse_off = np.concatenate([data[:,0,:peak_end[0][0]],
                            data[:,0,peak_end[0][1]:peak_end[1][0]],
                            data[:,0,peak_end[1][1]:]],
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
        L_base = np.nanmean(np.concatenate([L[peak_end[1][1]:],L[:peak_end[0][0]],
                                            L[peak_end[0][1]:peak_end[1][0]]],
                                 axis=-1))
        L = L-L_base
        stdL = np.nanstd(np.concatenate([L[peak_end[1][1]:],L[:peak_end[0][0]],
                                            L[peak_end[0][1]:peak_end[1][0]]],
                                 axis=-1))
        snrL = np.abs(L/stdL)
        
        pu_err = np.nanstd(np.concatenate([u[peak_end[1][1]:],u[:peak_end[0][0]],
                                            u[peak_end[0][1]:peak_end[1][0]]],
                                 axis=-1))
        pq_err = np.nanstd(np.concatenate([q[peak_end[1][1]:],q[:peak_end[0][0]],
                                            q[peak_end[0][1]:peak_end[1][0]]],
                                 axis=-1))
       
        #pu_err = np.nanstd(data[:,1,peak_end[0]:peak_end[1]],axis=0)
        #pq_err = np.nanstd(data[:,2,peak_end[0]:peak_end[1]],axis=0)
        
  
        
        for side,end in enumerate(peak_end):
            
            phase_peak = phase[end[0]:end[1]]
            
            pu = u[end[0]:end[1]]
            pq = q[end[0]:end[1]]
            
            ang=np.arctan2(pu, pq)/np.pi*180./2.
            ang[ang<0] += 180
            
            ang_err2 = ((pu*pq_err)**2 + (pq*pu_err)**2)/(pq**2 + pu**2)**2
            ang_err = np.sqrt(ang_err2)/np.pi*180/2
            
            idx = np.where(snrL[end[0]:end[1]]>3)
            
            ax[sub][side].errorbar(phase_peak[idx],ang[idx],ang_err[idx],c='k',fmt='.',linewidth=0.2,ms=1.5)
            ax[sub][side].set_ylim(max(min(ang)-5,-180),min(max(ang)+5,180))
            ax[sub][side].locator_params(axis='y',nbins=4)
            ax[sub][side].tick_params(labelbottom=False)
            #ax[sub][side].set_aspect(0.5)
            ax[-1][side].tick_params(labelbottom=True)
            
            ax[sub][side].set_xlim(phase_peak[0],phase_peak[-1])
            ax[sub][side].xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
            #ax[-1][side].set_xlabel('phase/$^\circ$')
            
        
        
        d = .01
        axis = ax[sub][0]
        
        axis.spines['right'].set_visible(False)
        axis.yaxis.tick_left()
        pos = text_pos[psr_name]
        ax[sub][0].text(pos[0],pos[1],'%d - %d MHz'%(channel[0],channel[1]),
                          verticalalignment=pos[2], 
                          horizontalalignment=pos[3],transform=ax[sub][0].transAxes)
        kwargs = dict(transform=axis.transAxes, color='k', clip_on=False)
        axis.plot((1-d,1+d),(-d,+d), **kwargs) # top-left diagonal
        axis.plot((1-d,1+d),(1-d,1+d), **kwargs) # bottom-left diagonal
        axis = ax[sub][1]
        axis.spines['left'].set_visible(False)
        
        #axis.yaxis.tick_right()
        axis.set_yticks([])
        
        kwargs.update(transform=axis.transAxes) # switch to the right axes
        axis.plot((-0.5*d,0.5*d),(-d,+d), **kwargs) # top-right diagonal
        axis.plot((-0.5*d,0.5*d),(1-d,1+d), **kwargs) # bottom-right diagonal
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.ylabel('PA/$^\circ$')     
    plt.xlabel('phase/$^\circ$')
    #plt.savefig('%s_PA.pdf'%psr_name)
    plt.show()