# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 13:26:18 2020

@author: Littl
"""
import os

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


import sys
sys.path.append('../')

import psr.cor as cor
import psr.drifting as drifting

abc = ['a','b','c']

folder_name = ['J1644-4559', 
               'J1327-6222','J1401-6357','J0837-4135',
                           'J1752-2806', 'J0738-4042']#,

folder_name = ['J1401-6357']

pulse_end = { 'J0738-4042':[590,800],'J0837-4135':[655,730], 'J0942-5552':[330,480],
             'J1327-6222':[65,130],'J1644-4559':[750,840],'J1401-6357':[80,135],
                                  'J1752-2806':[710,765],
            'J1901-0906':[710,880],'J1456-6843':[520,720],
            'J1651-4246':[375,750],'J1825-0935':[170,750]
        }

mannual =  { 
        'J0738-4042':[],'J0837-4135':[2339],
        'J0942-5552':[204,206,800,801,802,803,1201,1550,1551,1628,1629,1634,1635,1663,4789],
        'J1327-6222':[1866,2327,2328],'J1401-6357':[],'J1644-4559':[],
        'J1651-4246':[],'J1752-2806':[],'J1825-0935':[2330,2331,2333,2335,2447,2448,3284],
        'J1901-0906':[],'J1456-6843':[]
        
        }




plt.rcParams['font.sans-serif'] = ['Times New Roman'] #设置字体为罗马字体
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['text.usetex'] = False
plt.rcParams['mathtext.rm'] = 'Times New Roman'
plt.rcParams['mathtext.it'] = 'Times New Roman:italic'
plt.rcParams['mathtext.bf'] = 'Times New Roman:bold'
plt.rcParams['xtick.direction'] = 'in' #坐标轴刻度向内
plt.rcParams['ytick.direction'] = 'in'
#plt.rcParams['savefig.dpi'] = 500 #保存图片的分辨率
plt.rcParams['font.size']=30

plt.rcParams['xtick.bottom'] = True
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.left'] = True
plt.rcParams['ytick.right'] = True

plt.rcParams['xtick.labelbottom'] = False
plt.rcParams['xtick.labeltop'] = False
plt.rcParams['ytick.labelleft'] = False
plt.rcParams['ytick.labelright'] = False




def align_xaxis(axes):
       """Align zeros of the two axes, zooming them out by same ratio"""
       axes = np.array(axes)

       extrema = np.array([ax.get_xlim() for ax in axes])

       # reset for divide by zero issues
       for i in range(len(extrema)):
           if np.isclose(extrema[i, 0], 0.0):
               extrema[i, 0] = -1
           if np.isclose(extrema[i, 1], 0.0):
               extrema[i, 1] = 1

       # upper and lower limits
       lowers = extrema[:, 0].min()
       uppers = extrema[:, 1].max()

       extrema[:,0] = lowers
       extrema[:,1] = uppers
       # bump by 10% for a margin
       extrema[:, 0] *= 1.1
       extrema[:, 1] *= 1.1
           
       # set axes limits
       [axes[i].set_xlim(*extrema[i]) for i in range(len(extrema))]

for psr_name in folder_name:
    
    file_list = os.listdir(r'./%s_v2'%psr_name)
    
    file_list = file_list[::-1]
    
    for filename in file_list:
        
        # load data
        
        channel = filename.split('_')[-1].split('.')[0]
        channel = channel.split('-')[0],channel.split('-')[1].split('M')[0]
        channel = list(int(i) for i in channel)
        
        peak_end = pulse_end[psr_name]
        
        phase = np.linspace(0,360,1024)
        phase_peak = phase[peak_end[0]:peak_end[1]]
    
        data = np.load('.\%s_v2\\%s'%(psr_name,filename))
        base_data = ((((data[:,:,:peak_end[0]])).sum(-1) + ((data[:,:,peak_end[1]:])).sum(-1))/
                        (data[:,:,:peak_end[0]].shape[-1]+data[:,:,peak_end[1]:].shape[-1]))
        data = data - base_data[:,:,None]
       
        w = np.ones((data.shape[0],))
        mean = data[:,0,:].mean(-1)
        distb = np.where(mean<=0)
        s1 = data[:,0,:peak_end[0]].mean(-1)
        s2 = data[:,0,-20:].mean(-1)
        for q in range(1):
            w[np.argmax(s1)-2:np.argmax(s1)+2] = 0
            w[np.argmax(s2)-2:np.argmax(s2)+2] = 0
            s1[np.argmax(s1)] = 0
            s2[np.argmax(s2)] = 0
        w[distb] = 0
        w[mannual[psr_name]] = 0
    
        I_prof = np.nanmean(data[:,0,peak_end[0]:peak_end[1]],axis=0)
        Q_peak = np.nanmean(data[:,1,peak_end[0]:peak_end[1]],axis=0)
        U_peak = np.nanmean(data[:,2,peak_end[0]:peak_end[1]],axis=0)
        L_prof = np.sqrt(Q_peak**2+U_peak**2)
        V_prof = np.nanmean(data[:,3,peak_end[0]:peak_end[1]],axis=0)
        
        top = max(I_prof)
        I_prof = I_prof/top
        L_prof = L_prof/top
        V_prof = V_prof/top
        
        
        Q = data[:,1]
        U = data[:,2]
        L = np.sqrt(Q**2+U**2)
        base_L = np.concatenate((L[:,:peak_end[0]],L[:,peak_end[1]:]),axis=1)
        base_L = np.nanmean(base_L,axis=1)
        L = L- base_L[:,None]
        L = L/top
        
        I = data[:,0]/top
        abs_V = abs(data[:,-1]/top)
    
        data[np.where(w==0)] = np.nan
        
        u = np.nanmean(data[:,1,:],axis=0)#(w[:,None]*data[:,1,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
        q = np.nanmean(data[:,2,:],axis=0)#(w[:,None]*data[:,2,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
         
        
        stdL = np.nanstd(np.concatenate([L.mean(0)[peak_end[1]:],L.mean(0)[:peak_end[0]]],
                                 axis=-1))
        snrL = np.abs(L.mean(0)/stdL)
        
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

        
        # plot
        
        fig = plt.figure(figsize=(6,16))
    
        grid=plt.GridSpec(8,3,fig,wspace=0,hspace=0)
        
        # PA
        ax=plt.subplot(grid[0,0:2])     
        ax.errorbar(phase_peak[idx],ang[idx],ang_err[idx],c='k',fmt='.',linewidth=0.2,ms=1.5)
        ax.set_ylim(max(min(ang)-5,-180),min(max(ang)+5,180))    
        ax.set_yticks([30,60,90,120,150,180])
        ax.locator_params(axis='y',nbins=3)           
        ax.set_xlim(phase_peak[0],phase_peak[-1])
        ax.set_title("%s %d-%d MHz"%(psr_name,channel[0],channel[1]),fontsize=33,loc='left')
        ax.set_ylabel('PA/$^\circ$')
        plt.tick_params(labelleft=True)
         
        
        profile = {'I':I_prof,'L':L_prof,'V':V_prof}
        ls = {'I':['black','solid'],'L':['red','dashed'],'V':['blue','dotted']}
        columns = {'I':I,'L':L,'V':abs_V}
        
        # all profile
        ax=plt.subplot(grid[-1,0:2])     
        plt.sca(ax) #选择一个子图，进行后续设置
        for i,key in enumerate(profile):
            plt.plot(phase_peak,profile[key],label=key,c=ls[key][0],linewidth=1,linestyle=ls[key][1])     
        plt.locator_params(axis='x',nbins=3)
        plt.xlim(phase_peak[0],phase_peak[-1])
        plt.yticks([0,0.5,1])
        #plt.ylabel('normalised flux',fontsize=25)
        #plt.xlabel(r'phase/$^\circ$',fontsize=25)
        plt.legend(loc='lower left', bbox_to_anchor=(1.05,-0.2),borderaxespad = 0.0,ncol=1,frameon=False)
        #plt.tick_params(labelbottom=True)
        
        # profile respectively      
        axis = []
        for j,key in enumerate(profile):
            ax=plt.subplot(grid[1+j*2:3+j*2,-1])
            plt.sca(ax)
            plt.locator_params(axis='y',nbins=3)
            plt.locator_params(axis='x',nbins=3)
            plt.plot(profile[key]/max(abs(profile[key])),phase_peak,label=key,c=ls[key][0],
                         linestyle=ls[key][1],linewidth=1)
            plt.ylim(phase_peak[0],phase_peak[-1])
            ax.axes.xaxis.set_ticklabels([])
            plt.locator_params(axis='y',nbins=3)
            ax.yaxis.set_label_position("right")
            ax.set_ylabel(key,fontsize=35)
            axis.append(ax)
            
            d = columns[key]
            pulse_off1 = d[:,peak_end[0]-50:peak_end[0]]
            pulse_off2 = d[:,peak_end[1]:peak_end[1]+50]
            pulse_off = np.concatenate((pulse_off1,pulse_off2),axis=1)
            
            cr_off = cor.weighed_corr_map(w, pulse_off, 0,-10000)
    
            cr_off = np.concatenate([cr_off[:int(cr_off.shape[0]/3),-int(cr_off.shape[1]/3):],
                                    cr_off[-int(cr_off.shape[0]/3):,:int(cr_off.shape[1]/3)]])
  
            cr_off = np.std(cr_off)
            
            data_peak = d[:,peak_end[0]:peak_end[1]]
            cr = cor.weighed_corr_map(w, data_peak, 1,abs(cr_off))
            
            
            ax = plt.subplot(grid[1+j*2:3+j*2,0:2])
            plt.sca(ax)
            
            norm=mpl.colors.Normalize(vmin=-1,vmax=1) #将1设为最红，-1为最蓝，中间的0为白
            cmap = mpl.cm.seismic
            cmap.set_bad(color='#404040',alpha=0.15)
            corr_map = plt.imshow(cr,norm=norm, cmap=cmap, aspect='auto', origin='lower',
                                  extent=(phase_peak[0], phase_peak[-1], phase_peak[0], phase_peak[-1]))

            ax.axes.yaxis.set_ticklabels([])
            ax.tick_params(labeltop=False)
            plt.locator_params(axis='x',nbins=3)
            plt.locator_params(axis='y',nbins=3)
            
    
        align_xaxis(axis)   
    
        plt.savefig('.\\revised_v2\\cor_PA\\%s_%d-%d_corxPA_lag=1.pdf'%(psr_name,channel[0],channel[1]),bbox_inches='tight',pad_inches = 0)
        plt.show()