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

plt.rcParams["xtick.minor.visible"] =  False
plt.rcParams["ytick.minor.visible"] =  False

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
       extrema[i, 0] *= 1.05
       extrema[i, 1] *= 1.05
           
       # set axes limits
       [axes[i].set_ylim(*extrema[i]) for i in range(len(extrema))]
       
def align_xaxis(axes):
       
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
       extrema[i, 0] *= 1.05
       extrema[i, 1] *= 1.05
       
       # set axes limits
       [axes[i].set_xlim(*extrema[i]) for i in range(len(extrema))]
    

def plot_cor(w,data_all,p_end,lags,star_name,freq,bin_num=1024):
    
    phase = np.linspace(0,360,bin_num)
    
    fig=plt.figure(figsize=(14.3,28.4))
    g=plt.GridSpec(108,14,fig,wspace=0.1,hspace=0.4)
    
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
    
    Q = data_all[:,1]
    V = data_all[:,2]
    L = np.sqrt(Q**2+V**2)
    base_L = np.nanmean(np.concatenate([L[:,:p_end[0][0]],L[:,p_end[1][1]:],
                         L[:,p_end[0][1]:p_end[1][0]]],axis=1),axis=1)
    L = L- base_L[:,None]
    L = L/top
    
    I = data_all[:,0]/top    
    abs_V = abs(data_all[:,-1]/top)
    
    columns = {'I':I,'L':L,'V':abs_V}
    ls = {'I':['black','solid'],'L':['red','dashed'],'V':['blue','dotted']}
    profile = {'I':I_prof,'L':L_prof,'V':V_prof}
    
    axis = []
    
    jj = -1
    for key,data in columns.items():
        jj += 1
        
        # profile top
        gg = g[:8,4*jj:(jj+1)*4]
        subgrid = gg.subgridspec(1,3,hspace=0,wspace=0)
   
        ax_IP = plt.subplot(subgrid[:,0])
        plt.sca(ax_IP) 
        plt.locator_params(axis='y',nbins=3)
        plt.locator_params(axis='x',nbins=3)
        plt.plot(phase[p_end[0][0]:p_end[0][1]],profile[key][p_end[0][0]:p_end[0][1]]/max(abs(profile[key])),label=key,c=ls[key][0],
                 linestyle=ls[key][1],linewidth=1)
        plt.xlim(phase[p_end[0][0]],phase[p_end[0][1]])
        plt.yticks([0,0.5,1])
        plt.tick_params(right=False,labelbottom=False)
        #plt.tick_params(which='minor',right=False,labelbottom=False)
        ax_IP.axes.xaxis.set_ticklabels([])
        ax_IP.spines['right'].set_linestyle(dash)
        ax_IP.spines['right'].set_linewidth(0.5)
        ax_IP.spines['right'].set_color('k')
        axis.append(ax_IP)
        
        ax_MP = plt.subplot(subgrid[:,1:])
        plt.sca(ax_MP) 
        plt.locator_params(axis='y',nbins=3)
        plt.locator_params(axis='x',nbins=3)
        plt.plot(phase[p_end[1][0]:p_end[1][1]],profile[key][p_end[1][0]:p_end[1][1]]/max(abs(profile[key])),label=key,c=ls[key][0],
                 linestyle=ls[key][1],linewidth=1)
        plt.xlim(phase[p_end[1][0]],phase[p_end[1][1]])
        plt.tick_params(left=False)
        #plt.tick_params(which='minor',left=False)
        plt.yticks([0,0.5,1])
        ax_MP.axes.xaxis.set_ticklabels([])
        ax_MP.spines['left'].set_linestyle(dash)
        ax_MP.spines['left'].set_linewidth(0.5)
        ax_MP.spines['left'].set_color('k')
        d = 0.03
        kwargs = dict(transform=ax_MP.transAxes,color='k', clip_on=False)
        #ax_MP.plot((-d,+d),(-2*d,2*d),**kwargs)
        #ax_MP.plot((-d,+d),(1-2*d,1+2*d),**kwargs)
        axis.append(ax_MP)
    
        
        align_yaxis(np.array(axis))
        
        if key == 'I':
            ax_IP.tick_params(labelleft=True)
            ax_IP.locator_params(axis='y',nbins=3)
            ax_IP.set_ylabel('normalised flux')
        
        add = fig.add_subplot(g[:8,4*jj:(jj+1)*4], frameon=False)
        title = key
       # if key == 'V':
       #     title = r'$|$ V $|$'   
        add.set_title(title,fontsize=35,ha='center',va='bottom',position=(0.5,1.025))
        add.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
        add.minorticks_off()
        # correlation matrix
        for lag in lags:
            
            gs = g[8+lag*16:24+lag*16,jj*4:(jj+1)*4]
            
            tl = False
            
            if lag==lags[-1]:
                tl = True
        
            cor_block(gs,data,p_end,lag,tl=tl)
    
    align_yaxis(np.array(axis))    
    
    for lag in lags:
        add = fig.add_subplot(g[8+lag*16:24+lag*16,0:4], frameon=False)
        add.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False,
                        labelleft=True)
        add.minorticks_off()
        add.set_ylabel(r'$\mathrm{\mathit{n}_{\rm{lag}}}$=%d'%lag, fontsize=35,)
        add.yaxis.set_label_coords(-0.05, 0.5)
            
    for j in np.arange(6):
        
        sd_grid = g[8+j*16:8+(j+1)*16,-2:]
        side = sd_grid.subgridspec(3,1,hspace=0,wspace=0)
        
        axis = []
        
        
        ax_MP = plt.subplot(side[:-1,:])
        plt.sca(ax_MP) 
        plt.locator_params(axis='y',nbins=2)
        plt.locator_params(axis='x',nbins=3)
        for key in profile.keys():
            plt.plot(profile[key][p_end[1][0]:p_end[1][1]],
                     phase[p_end[1][0]:p_end[1][1]],label=key,c=ls[key][0],
                     linestyle=ls[key][1],linewidth=1)
            
        plt.ylim(phase[p_end[1][0]],phase[p_end[1][1]])
        plt.tick_params(bottom=False,labelright=True)
        #plt.tick_params(axis='x', which='minor', bottom=False)
        
        plt.yticks([225,250])
        ax_MP.axes.xaxis.set_ticklabels([])

        ax_MP.spines['bottom'].set_linewidth(0.5)
        ax_MP.spines['bottom'].set_color('k')
        ax_MP.spines['bottom'].set_linestyle(dash)
        
        
        d = 0.03
        

                
        ax_IP = plt.subplot(side[-1,:])
        plt.sca(ax_IP) 
        plt.locator_params(axis='y',nbins=1)
        plt.locator_params(axis='x',nbins=3)
        for key in profile.keys():
            plt.plot(profile[key][p_end[0][0]:p_end[0][1]],
                     phase[p_end[0][0]:p_end[0][1]],label=key,c=ls[key][0],
                     linestyle=ls[key][1],linewidth=1)
            
        plt.ylim(phase[p_end[0][0]],phase[p_end[0][1]])
        plt.xticks([0,0.5,1])
        plt.yticks([70])
        plt.tick_params(top=False,labelright=True)
        plt.tick_params(axis='x', which='minor', top=False)
        ax_IP.axes.xaxis.set_ticklabels([])
        
        ax_IP.spines['top'].set_linestyle(dash)
        ax_IP.spines['top'].set_linewidth(0.5)
        ax_IP.spines['top'].set_color('k')

        kwargs = dict(transform=ax_IP.transAxes,color='k', clip_on=False)
        #ax_IP.plot((-d,+d),(1-2*d,1+2*d),**kwargs)
        #ax_IP.plot((1-d,1+d),(1-2*d,1+2*d),**kwargs)
        
        axis.append(ax_MP)
        axis.append(ax_IP)
        
        align_xaxis(np.array(axis))
        
        if j==0:
            ax_MP.legend(loc='lower left', bbox_to_anchor=(0.05,1.02),borderaxespad = 0.0,ncol=1,frameon=False)
        
    frame = fig.add_subplot(g[:-4,:], frameon=False)
    frame.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    frame.minorticks_off()
    frame.xaxis.set_label_position("bottom")
    frame.set_xlabel('phase/$^\circ$',fontsize=35)
    frame.yaxis.set_label_position("right")
    frame.set_ylabel('phase/$^\circ$',fontsize=35)
    frame.yaxis.set_label_coords(1.08, 0.5)
    frame.xaxis.set_label_coords(0.5,-0.025,)
    
    ca = fig.add_subplot(g[-1,:-2], frameon=False)
    pos = ca.get_position()
    pos2 = [pos.x0, pos.y0-0.015, pos.width, pos.height]
    ca.set_position(pos2)
    plt.sca(ca)
    norm=mpl.colors.Normalize(vmin=-1,vmax=1) #将1设为最红，-1为最蓝，中间的0为白
    cmap = mpl.cm.seismic
    plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
             cax=ca, orientation='horizontal')
    
    al = fig.add_subplot(111,frameon=False)
    al.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    al.minorticks_off()
    al.set_title('PSR %s %s'%(star_name,freq),fontsize=35,
                 ha='center',va='bottom',position=(0.5,1.035))
    plt.savefig('%s_%s.pdf'%(star_name,freq.replace(' ','')),
                             bbox_inches='tight',pad_inches = 0)
    plt.show()

def cor_block(gs,data,p_end,lag,tl=False):
    
    ssgrid = gs.subgridspec(3,3,hspace=0,wspace=0)

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
        

    
    flux_peak = np.concatenate([data[:,p_end[0][0]:p_end[0][1]],data[:,p_end[1][0]:p_end[1][1]]],
                               axis=1)
    
    flux_peak = flux_peak
    cr_block = cor.weighed_corr_map(w, flux_peak, lag,abs(cr_off))
    cr = []
    cr.append(cr_block[-(p_end[1][1]-p_end[1][0]):,:(p_end[0][1]-p_end[0][0])])
    cr.append(cr_block[-(p_end[1][1]-p_end[1][0]):,(p_end[0][1]-p_end[0][0]):])
    cr.append(cr_block[:-(p_end[1][1]-p_end[1][0]),:(p_end[0][1]-p_end[0][0])])
    cr.append(cr_block[:-(p_end[1][1]-p_end[1][0]),(p_end[0][1]-p_end[0][0]):])
    

    ax_cr = [plt.subplot(ssgrid[:2,0]),plt.subplot(ssgrid[:2,1:]),
            plt.subplot(ssgrid[-1,0]),plt.subplot(ssgrid[-1,1:])]
    
    spine_set = [['bottom','right'],['left','bottom'],
                 ['top','right'],['top','left']]
    xt = {0:[70],1:[225,250],2:[70],3:[225,250]}
    yt = {0:[225,250],1:[225,250],2:[70],3:[70]}
    
    phase = np.linspace(0,360,1024)
    phase_IP = phase[p_end[0][0]:p_end[0][1]]
    phase_MP = phase[p_end[1][0]:p_end[1][1]]
    phase_peak = {0:(phase_IP[0],phase_IP[-1],phase_MP[0],phase_MP[-1]),
                  1:(phase_MP[0],phase_MP[-1],phase_MP[0],phase_MP[-1]),
                  2:(phase_IP[0],phase_IP[-1],phase_IP[0],phase_IP[-1]),
                  3:(phase_MP[0],phase_MP[-1],phase_IP[0],phase_IP[-1])}
    
    ax_cr[0].tick_params(bottom=False,right=False)
    #ax_cr[0].tick_params(which='minor', bottom=False,right=False)
    ax_cr[1].tick_params(bottom=False,left=False)
    #ax_cr[1].tick_params(which='minor', bottom=False,left=False)
    ax_cr[2].tick_params(top=False,right=False)
    #ax_cr[1].tick_params(which='minor', top=False,right=False)
    ax_cr[3].tick_params(top=False,left=False)
    #ax_cr[3].tick_params(which='minor', top=False,left=False)
    
    
    for b in range(4):
        plt.sca(ax_cr[b])
        for key in spine_set[b]:
            ax_cr[b].spines[key].set_linestyle(dash)
            ax_cr[b].spines[key].set_linewidth(0.5)
            ax_cr[b].spines[key].set_color('k')
            
            plt.xticks(xt[b])
            plt.yticks(yt[b])
            
    
        
        norm=mpl.colors.Normalize(vmin=-1,vmax=1) #将1设为最红，-1为最蓝，中间的0为白
        cmap = mpl.cm.seismic
        cmap.set_bad(color='#404040',alpha=0.15)
        plt.imshow(cr[b],norm=norm, cmap=cmap, aspect='auto', origin='lower',
                   extent=phase_peak[b])
    

    
    if tl == True:
        
        for bt in [2,3]:
            
            ax_cr[bt].tick_params(labelbottom=True)
            ax_cr[bt].locator_params(axis='x',nbins=bt-1)
            
        
    d = 0.03
    kwargs = dict(transform=ax_cr[0].transAxes,color='k', clip_on=False)
    #ax_cr[0].plot((-3/2*d,+3/2*d),(-d,d),**kwargs)
    #ax_cr[0].plot((1-3/2*d,1+3/2*d),(1-d,1+d),**kwargs)
    kwargs = dict(transform=ax_cr[-1].transAxes,color='k', clip_on=False)
    #ax_cr[-1].plot((-3/4*d,+3/4*d),(-2*d,2*d),**kwargs)
    #ax_cr[-1].plot((1-3/4*d,1+3/4*d),(1-2*d,1+2*d),**kwargs)

    
    
def plot_all(w,data,peak_end,pulse_end,lags,star_name,freq):
    
    
    plot_cor(w,data,peak_end,lags,star_name,freq,bin_num=1024)
       
        
    
    

lags = [0,1,2,3,4,5]
for psr_name in folder_name:
    
    file_list = os.listdir(r'.\\'+psr_name+'_v3')
    
    os.chdir('./revised_v3/cor/')
    
    file_list = file_list[::-1]
    for filename in file_list:
        #filename = 'J1456-6843_704-1727MHz.npy'
        data = np.load('..\..\%s_v3\\%s'%(psr_name,filename))
        
        peak_end = pulse_end[psr_name]
        left=peak_end[0]
        right=peak_end[1]        

        pulse_off = np.concatenate([data[:,0,:peak_end[0][0]],
                            data[:,0,peak_end[0][1]:peak_end[1][0]],
                            data[:,0,peak_end[1][1]:]],
                            axis=-1)
        base_data = np.nanmean(pulse_off,axis=-1)[:,None]
        data = data - base_data[:,:,None]
       
        pulse_off_V = np.concatenate([data[:,0,peak_end[0][0]-100:peak_end[0][0]],
                        data[:,0,peak_end[0][1]:peak_end[1][0]],
                        data[:,0,peak_end[1][1]:peak_end[1][1]+100]],
                        axis=-1)
        base_data_V = np.nanmean(pulse_off_V,axis=-1)
        data[:,-1] = data[:,-1] - base_data_V[:,None]
        
        w = np.ones((data.shape[0],))
        mean = data[:,0,:].mean(-1)
        distb = np.where(mean<=0)
        w[distb] = 0
        w[mannual[psr_name]] = 0
        data[distb] = np.nan
        data[mannual[psr_name]] = np.nan
        
        
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
        
        plot_all(w,data,peak_end,pulse_end,lags,psr_name,freq)

    os.chdir('../../')
        
