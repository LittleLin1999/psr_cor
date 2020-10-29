# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 16:26:55 2020

@author: Littl
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np




def weighed_corr(w1,w2,l1,l2,lag):
    
    n = len(l1)
    
 
    s1 = l1[lag:]
    s2 = l2[:n-lag]
    
    ww1 = w1[lag:]
    ww2 = w2[:n-lag]

    w = ww1*ww2
    
    y1 = w*s1
    y2 = w*s2
    
    avg1 = np.sum(y1)/np.sum(w)
    avg2 = np.sum(y2)/np.sum(w)

    cov = np.sum((y1-avg1)*(y2-avg2))
    
    var1 = np.sum((y1-avg1)**2)
    var2 = np.sum((y2-avg2)**2)
    
    
    coeff = cov/np.sqrt(var1*var2)

    return coeff

    
def weighed_corr_map(w,l,lags,sigma):

    ll = l.T
    cr = np.zeros((len(ll),len(ll)))
    #delta = np.zeros((len(ll),len(ll)))
    for j in range(len(ll)):
        for k in range(len(ll)):
            coe = weighed_corr(w,w,ll[j], ll[k], lags)
            if abs(coe) >= 3*sigma:
                cr[j][k] = coe
            else:
                cr[j][k] = np.nan
            #delta[j][k]=delta_c
    return cr#,delta

def plot_cor(w,data,flux,peak_end,lag,star_name,freq,bin_num=1024):
    
    phase = np.linspace(0,360,bin_num)
    
    flux_peak = data[:,0,peak_end[0]:peak_end[1]]
    
    I_prof = data[:,0,peak_end[0]:peak_end[1]].mean(0)
    top = max(I_prof)
    I_prof = I_prof/top
    flux_peak = flux_peak/top
    
    phase_peak = phase[peak_end[0]:peak_end[1]]
    
    pulse_off1 = data[:,0,peak_end[0]-50:peak_end[0]]
    pulse_off2 = data[:,0,peak_end[1]:peak_end[1]+50]
    
    pulse_off = np.concatenate((pulse_off1,pulse_off2),axis=1)
    
    
    cr_off = weighed_corr_map(w, pulse_off, lag,-10000)
    
    
    cr_off = np.concatenate([cr_off[:int(cr_off.shape[0]/3),-int(cr_off.shape[1]/3):],
                                    cr_off[-int(cr_off.shape[0]/3):,:int(cr_off.shape[1]/3)]])
 
    cr_off = np.std(cr_off)
    
    
    cr = weighed_corr_map(w, flux_peak, lag,-1000)
    
 
    
    plt.rcParams['figure.figsize']=(9,9)
    
    
    
    fig=plt.figure()
    
    grid=plt.GridSpec(18,18,fig,wspace=0,hspace=0)
    
    ax1 =  plt.subplot(grid[:6,6:])
    plt.sca(ax1)
    plt.locator_params(axis='y',nbins=4)
    plt.title(r'%s    %s    $n_{\rm{lag}}$ = %d' % (star_name,freq, lag))
    
   
    
    plt.plot(phase_peak,I_prof,label='I',c='black',linewidth=1)
    plt.xlim(phase_peak[0],phase_peak[-1])
    plt.tick_params(labelbottom=False)
    plt.legend()
    plt.ylabel('normalised flux')
    
    ax2 =  plt.subplot(grid[6:,:6])
    plt.sca(ax2)
    plt.plot(I_prof,phase_peak,c='black',linewidth=1)
   
    ax2.invert_xaxis()
    plt.ylim(phase_peak[0],phase_peak[-1])
    plt.locator_params(axis='x',nbins=4)
    plt.ylabel('phase/$^\circ$' )
    
    ax3 = plt.subplot(grid[6:,6:])
    
    plt.sca(ax3)
    
    norm=mpl.colors.Normalize(vmin=-1,vmax=1) #将1设为最红，-1为最蓝，中间的0为白
    cmap = mpl.cm.seismic
    cmap.set_bad(color='#404040',alpha=0.15)
    corr_map = plt.imshow(cr,norm=norm, cmap=cmap, aspect='auto', origin='lower',
           extent=(phase_peak[0], phase_peak[-1], phase_peak[0], phase_peak[-1]))
    
    
    cbaxes = ax3.add_axes([0, 1, 0.03, 0.8])
    plt.colorbar(corr_map,orientation="horizontal",cax=cbaxes)
    plt.xlabel('phase/$^\circ$' )
    plt.yticks([])
    

    #fig = plt.figure(figsize=(18,0.5))
    #grid = plt.GridSpec(1,18,fig,wspace=0,hspace=0)
    #ax = plt.subplot(grid[:,:])
    #plt.colorbar(corr_map,cax=ax,orientation='horizontal')
    #plt.savefig('colorbar%s.pdf'%(flux), bbox_inches = 'tight')
    return abs(cr_off)

peak_end = [1400,2000]
data = np.load('J1644-4559_ntsc.npy')

base_data = ((((data[:,:,:peak_end[0]])).sum(-1) + ((data[:,:,peak_end[1]:])).sum(-1))/
                (data[:,:,:peak_end[0]].shape[-1]+data[:,:,peak_end[1]:].shape[-1]))
data = data - base_data[:,:,None]
   
w = np.ones((data.shape[0],))
mean = data[:,0,:].mean(-1)
distb = np.where(mean<=0)
w[distb] = 0

prof = data.mean(0)[0]
base  = (data[:,:,:peak_end[0]].mean(-1))[:,0]
w[np.where(base>(np.max(prof)/4))]=0

for lag in [0]:
    plot_cor(w,data,'I',peak_end,lag,'J1644-4559','1400 MHz',bin_num=4096)