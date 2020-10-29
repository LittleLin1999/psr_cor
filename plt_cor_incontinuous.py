# -*- coding: utf-8 -*-
"""
Created on Sun Apr  5 23:05:54 2020

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

plt.rcParams['font.sans-serif']=['Times New Roman'] #设置字体为罗马字体
plt.rcParams['xtick.direction'] = 'in' #坐标轴刻度向内
plt.rcParams['ytick.direction'] = 'in'

dash = (0,(26,14,26,14))
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
       extrema[i, 0] *= 1.1
       extrema[i, 1] *= 1.1
           
       # set axes limits
       [axes[i].set_ylim(*extrema[i]) for i in range(len(extrema))]
       
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
       lowers = extrema[:, 1].min()
       uppers = extrema[:, 0].max()
       
       extrema[:,1] = lowers
       extrema[:,0] = uppers
       # bump by 10% for a margin
       extrema[i, 0] *= 1.1
       extrema[i, 1] *= 1.1
       
       # set axes limits
       [axes[i].set_xlim(*extrema[i]) for i in range(len(extrema))]
       
       
def plot_cor(w,data,data_all,flux,p_end,lag,star_name,freq,bin_num=1024):
    
    phase = np.linspace(0,360,bin_num)
    
    fig=plt.figure(figsize=(9,9))
    grid=plt.GridSpec(18,18,fig,wspace=0.,hspace=0.)
    
    pulse_off11 = data[:,p_end[0][0]-50:p_end[0][0]]
    pulse_off12 = data[:,p_end[0][1]:p_end[0][1]+50]
    pulse_off21 = data[:,p_end[1][0]-50:p_end[1][0]]
    pulse_off22 = data[:,p_end[1][1]:p_end[1][1]+50]
    
    pulse_off = np.concatenate((pulse_off11,pulse_off12,pulse_off21,pulse_off22),
                               axis=1)
    cr_off = cor.weighed_corr_map(w, pulse_off, lag,-10000)
    
    cr_off = np.concatenate([cr_off[:int(cr_off.shape[0]/3),-int(cr_off.shape[1]/3):],
                                    cr_off[-int(cr_off.shape[0]/3):,:int(cr_off.shape[1]/3)]])
    cr_off = np.std(cr_off)
        
    I_prof = data_all[:,0,:].mean(0)#(w[:,None]*data_all[:,0,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
    Q_peak = data_all[:,1,:].mean(0)#(w[:,None]*data[:,1,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
    U_peak = data_all[:,2,:].mean(0)#(w[:,None]*data[:,2,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
    L_prof = np.sqrt(Q_peak**2+U_peak**2)
    L_base = np.nanmean(np.concatenate([L_prof[:p_end[0][0]],L_prof[p_end[1][1]:],
                         L_prof[p_end[0][1]:p_end[1][0]]],axis=0))
    L_prof = L_prof-L_base
    V_prof = (w[:,None]*data_all[:,3,:]).sum(0)/w.sum(0)
    
    top = max(I_prof)
    I_prof = I_prof/top
    L_prof = L_prof/top
    V_prof = V_prof/top
    
    flux_peak = np.concatenate([data[:,p_end[0][0]:p_end[0][1]],data[:,p_end[1][0]:p_end[1][1]]],
                               axis=1)
    
    flux_peak = flux_peak/top
    cr_block = cor.weighed_corr_map(w, flux_peak, lag,abs(cr_off))
    cr = []
    cr.append(cr_block[-(p_end[1][1]-p_end[1][0]):,:(p_end[0][1]-p_end[0][0])])
    cr.append(cr_block[-(p_end[1][1]-p_end[1][0]):,(p_end[0][1]-p_end[0][0]):])
    cr.append(cr_block[:-(p_end[1][1]-p_end[1][0]),:(p_end[0][1]-p_end[0][0])])
    cr.append(cr_block[:-(p_end[1][1]-p_end[1][0]),(p_end[0][1]-p_end[0][0]):])
    #plt.show()
    ax = []
    
    for side,end in enumerate(p_end):

        phase_peak = phase[end[0]:end[1]]
        if side==0:
            axis =  plt.subplot(grid[:6,6:10]),plt.subplot(grid[-4:,:6])
            
        else:
            axis =  plt.subplot(grid[:6,10:]),plt.subplot(grid[6:14,:6])
        ax.append(axis)
        plt.sca(ax[side][0])
        plt.locator_params(axis='y',nbins=4)
        #plt.title(r'%s    %s    $n_{\rm{lag}}$ = %d' % (flux,freq, lag))
    
        plt.plot(phase_peak,I_prof[end[0]:end[1]],label='I',c='black',linewidth=1)
        plt.plot(phase_peak,L_prof[end[0]:end[1]],label='L',c='red',linestyle='dashed',linewidth=1)
        plt.plot(phase_peak,V_prof[end[0]:end[1]],label='V',c='blue',linestyle='dotted',linewidth=1)
        plt.xlim(phase_peak[0],phase_peak[-1])
        plt.tick_params(labelbottom=False)
        plt.xticks(phase[end[0]:end[1]][15:-5])
        plt.locator_params(axis='x',nbins=side+3)
        if side==0:
            plt.ylabel('normalised flux')
            
        
          
        plt.sca(ax[side][1])
        plt.plot(I_prof[end[0]:end[1]],phase_peak,c='black',linewidth=1)
        plt.plot(L_prof[end[0]:end[1]],phase_peak,c='red',linestyle='dashed',linewidth=1)
        plt.plot(V_prof[end[0]:end[1]],phase_peak,c='blue',linestyle='dotted',linewidth=1)
        ax[side][1].invert_xaxis()
        plt.ylim(phase_peak[0],phase_peak[-1])
        plt.locator_params(axis='x',nbins=4)
        plt.yticks(phase[end[0]:end[1]][15:-5])
        plt.locator_params(axis='y',nbins=side+3)
        ax[side][1].yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
        #plt.ylabel('phase/$^\circ$' )
    align_yaxis([ax[0][0],ax[1][0]])
    align_xaxis([ax[0][1],ax[1][1]])
    ax[1][0].legend()
    
    axis_set = {'right':[0,0],'top':[0,1],'left':[1,0],'bottom':[1,1]}
    for key in axis_set:
        val = axis_set[key]
        axis = ax[val[0]][val[1]]
        axis.spines[key].set_linestyle(dash)
        axis.spines[key].set_linewidth(0.5)
        axis.spines[key].set_color('k')
        if val==[1,1]:
            axis.set_xticks([])
        if val==[1,0]:
            axis.set_yticks([])
    d = .01
    kwargs = dict(transform=ax[1][0].transAxes,color='k', clip_on=False)
    #ax[1][0].plot((-d,+d),(-2*d,2*d),**kwargs)
    ax[1][0].plot((-6/8*d,+6/8*d),(1-2*8/6*d,1+2*8/6*d),**kwargs)
    kwargs = dict(transform=ax[1][1].transAxes,color='k', clip_on=False)
    ax[1][1].plot((-d,+d),(-2*d,2*d),**kwargs)
    #ax[1][1].plot((1-d,1+d),(-2*d,+2*d),**kwargs)
    ax_cr = [plt.subplot(grid[6:14,6:10]),plt.subplot(grid[6:14,10:]),
            plt.subplot(grid[-4:,6:10]),plt.subplot(grid[-4:,10:])]
    
    spine_set = [['bottom','right'],['left','bottom'],
                 ['top','right'],['top','left']]
    for b in range(4):
        plt.xticks([])
        plt.yticks([])
        plt.sca(ax_cr[b])
        for key in spine_set[b]:
            ax_cr[b].spines[key].set_linestyle(dash)
            ax_cr[b].spines[key].set_linewidth(0.5)
            ax_cr[b].spines[key].set_color('k')

        phase_peak = phase[p_end[int(b/3)][0]:p_end[int(b/3)][1]]
        
        norm=mpl.colors.Normalize(vmin=-1,vmax=1) #将1设为最红，-1为最蓝，中间的0为白
        cmap = mpl.cm.seismic
        cmap.set_bad(color='#404040',alpha=0.15)
        plt.imshow(cr[b],norm=norm, cmap=cmap, aspect='auto', origin='lower',
                   extent=(phase_peak[0], phase_peak[-1], phase_peak[0], phase_peak[-1]))
    
    for bt in [2,3]:
        ax_cr[bt].set_xticks(phase[p_end[int(bt/3)][0]:p_end[int(bt/3)][1]][15:-5])
        ax_cr[bt].locator_params(axis='x',nbins=bt+1)
        ax_cr[bt].xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))

    kwargs = dict(transform=ax_cr[0].transAxes,color='k', clip_on=False)
    ax_cr[0].plot((-6/4*d,+6/4*d),(-2*d,2*d),**kwargs)
    ax_cr[0].plot((1-6/4*d,1+6/4*d),(1-2*d,1+2*d),**kwargs)
    kwargs = dict(transform=ax_cr[-1].transAxes,color='k', clip_on=False)
    ax_cr[-1].plot((-6/8*d,+6/8*d),(-4*d,4*d),**kwargs)
    #ax_cr[-1].plot((1-d,1+d),(1-2*d,1+2*d),**kwargs)
    add = fig.add_subplot(grid[6:,6:], frameon=False)
    add.set_xlabel('phase/$^\circ$' )
    add.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)

    add = fig.add_subplot(grid[6:,:6], frameon=False)
    add.set_ylabel('phase/$^\circ$',labelpad=12)
    add.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    
    add = fig.add_subplot(grid[:6,6:], frameon=False)
    add.set_title(r'%s    %s    $n_{\rm{lag}}$ = %d' % (flux,freq, lag))
    add.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)

    if flux == r'$|$ V $|$':
        flux = 'abs_V'
    #plt.savefig('%s_%s_lag=%d.pdf'%(flux,freq.replace(' ',''),lag))
    plt.show()
   
    
def plot_all(w,data,peak_end,pulse_end,lags,star_name,freq):
    
    
    for i in lags:
        
        plot_cor(w,data[:,0],data,'I',peak_end,i,star_name,freq,bin_num=1024)
        Q = data[:,1]
        V = data[:,2]
        L = np.sqrt(Q**2+V**2)
        base_L = np.nanmean(np.concatenate([L[:,:peak_end[0][0]],L[:,peak_end[1][1]:],
                                            L[:,peak_end[0][1]:peak_end[1][0]]],axis=-1),axis=-1)
        L = L- base_L[:,None]
        plot_cor(w,L,data,'L',peak_end,i,star_name,freq,bin_num=1024)
       
        plot_cor(w,abs(data[:,-1]),data,r'$|$ V $|$',peak_end,i,star_name,freq,bin_num=1024)
        
    
    

lags = [0,1,2,3,4,5]
for psr_name in folder_name:
    
    file_list = os.listdir(r'.\\'+psr_name+'_v2')
    
    os.chdir('./revised_v2/cor/')
    
    for filename in file_list:
        #filename = 'J1456-6843_704-1727MHz.npy'
        data = np.load('..\..\%s_v2\\%s'%(psr_name,filename))
        
        peak_end = pulse_end[psr_name]
        left=peak_end[0]
        right=peak_end[1]        

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
        
        
        freq = filename.split('.')[0].split('_')[-1]
        freq_chn = []
        freq_chn.append(int(freq.split('-')[0]))
        freq_chn.append(int(freq.split('-')[-1].split('M')[0]))
        for c,channel in enumerate(freq_chn):
            if channel%2 == 1:
                freq_chn[c] += 1
        freq = '%d - %d MHz'%(freq_chn[0],freq_chn[1])
        
        plot_all(w,data,peak_end,pulse_end,lags,psr_name,freq)

    os.chdir('../../')
        
