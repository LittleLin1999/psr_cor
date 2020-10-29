# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 09:02:36 2020

@author: Littl
"""
import os

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

abc = ['a','b','c']

folder_name = ['J0837-4135', 'J0942-5552',
               'J1327-6222','J1401-6357','J1644-4559',
                           'J1752-2806',
               'J1456-6843',
               'J1651-4246','J0738-4042']#,
folder_name = ['J0738-4042']

pulse_end = { 'J0738-4042':[420,625],'J0837-4135':[650,735], 'J0942-5552':[330,480],
             'J1327-6222':[65,130],'J1644-4559':[750,840],'J1401-6357':[80,135],
                                  'J1752-2806':[710,765],
            'J1901-0906':[710,880],'J1456-6843':[520,720],
            'J1651-4246':[375,750],'J1825-0935':[170,750]
        }
#'J0738-4042':[590,800]

mannual =  { 
        'J0738-4042':[2038, 2318, 2496, 2865, 2947, 3154, 3514, 5804, 5893, 7736, 8964,
        8965, 8966, 8967, 9300, 9301, 9304, 9306, 9307],
        'J0837-4135':list(range(2400,2500)),
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



flux = ['I','L','V']

def align_yaxis(axes):
       """Align zeros of the two axes, zooming them out by same ratio"""
       axes = np.array(axes)

       extrema = np.array([ax.get_ylim() for ax in axes])

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
       [axes[i].set_ylim(*extrema[i]) for i in range(len(extrema))]

for psr_name in folder_name:
    
    file_list = os.listdir(r'./%s_v3'%psr_name)
    
    
    #file_list[0], file_list[1], file_list[-1] = file_list[-1], file_list[0], file_list[1]
    file_list = file_list[::-1]
    
    
    peak_end = pulse_end[psr_name]
    
    phase = np.linspace(0,360,1024)
    phase_peak = phase[peak_end[0]:peak_end[1]]
    
    fig = plt.figure(figsize=(5*3+2,12))

    grid=plt.GridSpec(4,3,fig,wspace=0,hspace=0)
    ax = []
    for j in range(3):
        axis=plt.subplot(grid[0:1,j]),plt.subplot(grid[1:2,j]),plt.subplot(grid[2:,j]) 
        ax.append(axis)
        
    for index,filename in enumerate(file_list):
        if psr_name == 'J1644-4559':
            peak_end = [750,840]
            phase_peak = phase[peak_end[0]:peak_end[1]]
        channel = filename.split('_')[-1].split('.')[0]
        channel = channel.split('-')[0],channel.split('-')[1].split('M')[0]
        channel = list(channel)
        for c,chn in enumerate(channel):
            chn = int(chn)
            if chn%2==1:
                chn+=1
            channel[c] = chn
        if channel[0] > channel[1]:
            channel.reverse()
        if channel[0] < 1000:
            if psr_name == 'J1644-4559':
                peak_end = [750,840]
                phase_peak = phase[peak_end[0]:peak_end[1]]
            sub = 0
        if channel[0] > 1000 and channel[0] < 2000:
            sub = 1
        if channel[0] > 2500:
            sub = 2
        data = np.load('.\%s_v3\\%s'%(psr_name,filename))
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
        #data[distb] = np.nan
        #data[mannual[psr_name]] = np.nan
        
    
        I_peak = data[:,0,peak_end[0]:peak_end[1]]
        I_prof = np.nanmean(I_peak,axis=0)
        top = max(I_prof)
        I_prof = I_prof/top
        
        Q_peak = np.nanmean(data[:,1,peak_end[0]:peak_end[1]],axis=0)#(w[:,None]*data[:,1,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
        U_peak = np.nanmean(data[:,2,peak_end[0]:peak_end[1]],axis=0)#(w[:,None]*data[:,2,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
        L_prof = np.sqrt(Q_peak**2+U_peak**2)
        
    
        Q_base = np.nanmean(data[:,1,:],axis=0)#(w[:,None]*data[:,1,:]).sum(0)/w.sum(0)
        U_base = np.nanmean(data[:,2,:],axis=0)#(w[:,None]*data[:,2,:]).sum(0)/w.sum(0)
        L = np.sqrt(Q_base**2+U_base**2)
        L_base = ((((L[peak_end[0]-50:peak_end[0]])).sum(-1) + ((L[peak_end[1]:peak_end[1]+50])).sum(-1))/
                        (L[peak_end[0]-50:peak_end[0]].shape[-1]+L[peak_end[1]:peak_end[1]+50].shape[-1]))
        L_prof = L_prof - L_base 
        L_prof = L_prof/top
    
        V_prof = np.nansum(w[:,None]*data[:,3,peak_end[0]:peak_end[1]],axis=0)/w.sum(0)
        V_prof = V_prof/top
    

        u = np.nanmean(data[:,1,:],axis=0)#(w[:,None]*data[:,1,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
        q = np.nanmean(data[:,2,:],axis=0)#(w[:,None]*data[:,2,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)

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

        
        ax[sub][0].errorbar(phase_peak[idx],ang[idx],ang_err[idx],c='k',fmt='.',linewidth=0.2,ms=1.5)
        ax[sub][0].set_ylim(max(min(ang)-5,-180),min(max(ang)+5,180))
        
        ax[sub][0].set_yticks([30,60,90,120,150,180])
        ax[sub][0].locator_params(axis='y',nbins=3)

        
        ax[sub][0].set_xlim(phase_peak[0],phase_peak[-1])
        ax[sub][0].set_title("%d-%d MHz"%(channel[0],channel[1]),fontsize=33)
     
        
        plt.sca(ax[sub][1]) #选择一个子图，进行后续设置
        plt.plot(phase_peak,I_prof,label='I',c='black',linewidth=1)
        plt.plot(phase_peak,L_prof,label='L',c='red',linewidth=1,linestyle='dashed')
        plt.plot(phase_peak,V_prof,label='V',c='blue',linewidth=1,linestyle='dotted')
        plt.locator_params(axis='x',nbins=3)
        plt.xlim(phase_peak[0],phase_peak[-1])
        plt.yticks([0,0.5,1])
        
        if sub==0:
            plt.tick_params(labelleft=True)
           # kwargs = dict(transform=ax[sub][0].transAxes,fontsize=25,ha='right',va='center')
            plt.ylabel('normalised flux',fontsize=25)
        
            
        plt.sca(ax[sub][2])
        norm = mpl.colors.Normalize(vmin=0, vmax=4*top) #如果画出来的图 颜色映射不合适的话，可以手动调整。参考值：vmax取平均脉冲轮廓最大值的4倍
        current_cmap = mpl.cm.binary
        current_cmap.set_bad(color='white')
        I_peak[np.where(w==0)]=np.nan
        plt.imshow(I_peak,norm=norm,cmap=current_cmap,aspect='auto',origin='lower',extent=(phase_peak[0],phase_peak[-1], 0, I_peak.shape[0]))
        plt.tick_params(labelbottom=True)
        plt.locator_params(axis='x',nbins=3)
        #plt.gca().xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
        
        
        
        if sub == 0:
            plt.tick_params(labelleft=True)
            plt.locator_params(axis='y',nbins=3)
            plt.ylabel('pulse number')

        for jj in range(3):
            ax[sub][jj].xaxis.set_minor_locator(AutoMinorLocator(5))
            ax[sub][jj].yaxis.set_minor_locator(AutoMinorLocator(5))
        
    ax[1][-1].set_xlabel('phase/$^\circ$')
    ax[0][0].set_ylabel('PA/$^\circ$')
    if len(file_list) == 2:
        ax[1][1].legend(loc='upper left',bbox_to_anchor=(2.0,1.0),
                      frameon=False,fontsize=25)
    
    else:
        ax[-1][1].legend(loc='upper left',bbox_to_anchor=(1.0,1.0),
                          frameon=False,fontsize=25)
    kwargs = dict(transform=ax[1][0].transAxes,fontsize=35,ha='center',va='bottom')
    ax[1][0].text(x=0.5,y=1.35,s='PSR %s'%psr_name,**kwargs)
    ax[0][0].tick_params(labelleft=True)
    align_yaxis(np.array(ax)[:len(file_list),0])
    align_yaxis(np.array(ax)[:len(file_list),1])
    
    if psr_name == 'J1752-2806':
        for i in range(3):
            for j in range(3):
                ax[i][j].set_xticks([255,265])
    plt.savefig('.\\revised_v3\\profile\\%s_profile.pdf'%(psr_name),bbox_inches='tight',pad_inches = 0)
    plt.show()