# -*- coding: utf-8 -*-
"""
Created on Sat Oct 10 14:55:26 2020

@author: Littl
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

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

plt.rcParams['xtick.labelbottom'] = True
plt.rcParams['xtick.labeltop'] = True
plt.rcParams['ytick.labelleft'] = True
plt.rcParams['ytick.labelright'] = True

plt.rcParams["xtick.minor.visible"] =  True
plt.rcParams["ytick.minor.visible"] =  True

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






def cor_matrix(w,flux,peak_end,lags,bin_num=1024):
    

    
    flux_peak = flux[:,peak_end[0]:peak_end[1]]
    
    prof = flux[:,peak_end[0]:peak_end[1]].mean(0)
    top = max(prof)
    prof = prof/top
    flux_peak = flux_peak/top
    
    
    pulse_off1 = flux[:,peak_end[0]-50:peak_end[0]]
    pulse_off2 = flux[:,peak_end[1]:peak_end[1]+50]
    
    pulse_off = np.concatenate((pulse_off1,pulse_off2),axis=1)
    
    
    cr_set = []
    for lag in lags:
        cr_off = weighed_corr_map(w, pulse_off, lag,-10000)
    
    
        cr_off = np.concatenate([cr_off[:int(cr_off.shape[0]/3),-int(cr_off.shape[1]/3):],
                                    cr_off[-int(cr_off.shape[0]/3):,:int(cr_off.shape[1]/3)]])
 
        cr_off = np.std(cr_off)
    
    
        cr = weighed_corr_map(w, flux_peak, lag, abs(cr_off))
    
        cr_set.append(cr)
        
    return cr_set

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
#%% ntsc
lags = [0,1,2,3,4,5]
    
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
    
cor_ntsc = cor_matrix(w,data[:,0,:],peak_end,lags,4096)


I_ntsc = w[:,None]*data[:,0,peak_end[0]:peak_end[1]]
I_ntsc = I_ntsc.mean(0)
I_ntsc = I_ntsc/max(I_ntsc)

#%% km

lags = [0,1,2,3,4,5]
    
peak_end = [150,250]
data = np.load('J1644-4559_km.npy')

base_data = ((((data[:,:,:peak_end[0]])).sum(-1) + ((data[:,:,peak_end[1]:])).sum(-1))/
                (data[:,:,:peak_end[0]].shape[-1]+data[:,:,peak_end[1]:].shape[-1]))
data = data - base_data[:,:,None]
   
w = np.ones((data.shape[0],))
mean = np.nanmean(data[:,0,:],axis=-1)
distb = np.where(mean<=0)
w[distb] = 0

prof = np.nanmean(data,axis=0)[0]
base  = np.nanmean(data[:,:,:peak_end[0]],axis=-1)[:,0]
w[np.where(base>(np.max(prof)/4))]=0
    
cor_km = cor_matrix(w,data[:,0,:],peak_end,lags,4096)


I_km = w[:,None]*data[:,0,peak_end[0]:peak_end[1]]
I_km = np.nanmean(I_km,axis=0)
I_km = I_km/max(I_km)

#%% parkes1
peak_end = [750,840]
data = np.load('J1644-4559_1399-1799MHz.npy')
base_data = ((((data[:,:,:peak_end[0]])).sum(-1) + ((data[:,:,peak_end[1]:])).sum(-1))/
                (data[:,:,:peak_end[0]].shape[-1]+data[:,:,peak_end[1]:].shape[-1]))
data = data - base_data[:,:,None]
   
w = np.ones((data.shape[0],))
mean = np.nanmean(data[:,0,:],axis=-1)
distb = np.where(mean<=0)
w[distb] = 0

prof = np.nanmean(data,axis=0)[0]
base  = np.nanmean(data[:,:,:peak_end[0]],axis=-1)[:,0]
w[np.where(base>(np.max(prof)/4))]=0
    
cor_parkes1 = cor_matrix(w,data[:,0,:],peak_end,lags,1024)


I_p1 = w[:,None]*data[:,0,peak_end[0]:peak_end[1]]
I_p1 = np.nanmean(I_p1,axis=0)
I_p1 = I_p1/max(I_p1)

#%% parkes2
peak_end = [750,840]
data = np.load('J1644-4559_2172-2302MHz.npy')
base_data = ((((data[:,:,:peak_end[0]])).sum(-1) + ((data[:,:,peak_end[1]:])).sum(-1))/
                (data[:,:,:peak_end[0]].shape[-1]+data[:,:,peak_end[1]:].shape[-1]))
data = data - base_data[:,:,None]
   
w = np.ones((data.shape[0],))
mean = np.nanmean(data[:,0,:],axis=-1)
distb = np.where(mean<=0)
w[distb] = 0

prof = np.nanmean(data,axis=0)[0]
base  = np.nanmean(data[:,:,:peak_end[0]],axis=-1)[:,0]
w[np.where(base>(np.max(prof)/4))]=0
    
cor_parkes2 = cor_matrix(w,data[:,0,:],peak_end,lags,1024)


I_p2 = w[:,None]*data[:,0,peak_end[0]:peak_end[1]]
I_p2 = np.nanmean(I_p2,axis=0)
I_p2 = I_p2/max(I_p2)


#%% plot

I_prof = [I_ntsc, I_p1, I_km, I_p2]


plt.rcParams['figure.figsize']=(16.3,28.4)
fig=plt.figure()
grid=plt.GridSpec(108,4,fig,wspace=0.02,hspace=0.2)

phase_1024 = np.linspace(0,360,1024)
phase_4096 = np.linspace(0,360,4096)
phase_peak = [  phase_4096[1400:2000],
                phase_1024[750:840],
                phase_1024[150:250],
                phase_1024[750:840]  ]

title=['1400-1800 MHz',  # ntsc
       '1400-1800 MHz',
       '2302-2173 MHz',  # yn      
       '2302-2173 MHz']
name=['Haoping','Parkes','Yunnan','Parkes']

axis = []
for lag in np.arange(6):
 
    block = cor_ntsc[lag],cor_parkes1[lag] ,cor_km[lag], cor_parkes2[lag]
    
    for jj in range(len(block)):
        
        ax1 = plt.subplot(grid[:8,jj:(jj+1)])
        plt.sca(ax1)
        plt.locator_params(axis='y',nbins=3)
        plt.locator_params(axis='x',nbins=3)
        plt.plot(phase_peak[jj],I_prof[jj],'k',linewidth=1)
        plt.xlim(phase_peak[jj][0],phase_peak[jj][-1])
        plt.title("%s \n %s"%(name[jj],title[jj]),fontsize=20)
        ax1.axes.xaxis.set_ticklabels([])
    
        axis.append(ax1)
    
        #plt.title(title,fontsize=35,ha='center',va='bottom',position=(0.5,1.025))
        if jj == 0:
            ax1.yaxis.set_label_position("left")
            ax1.set_ylabel('normalised flux')
            ax1.tick_params(labelright=False)
            
        else:
            ax1.axes.yaxis.set_ticklabels([])
        
        ax = plt.subplot(grid[8+lag*16:24+lag*16,jj:(jj+1)])
        plt.sca(ax)
    
        norm=mpl.colors.Normalize(vmin=-1,vmax=1) #将1设为最红，-1为最蓝，中间的0为白
        cmap = mpl.cm.seismic
        cmap.set_bad(color='#404040',alpha=0.15)
        corr_map = plt.imshow(block[jj],norm=norm, cmap=cmap, aspect='equal', origin='lower',
                              extent=(phase_peak[jj][0], phase_peak[jj][-1], phase_peak[jj][0], phase_peak[jj][-1]))

        
        ax.tick_params(labeltop=False)
        ax.tick_params(labelleft=False)
        
        plt.locator_params(axis='x',nbins=3)
        plt.locator_params(axis='y',nbins=3)
        if jj == 0:
            plt.ylabel(r'$\mathrm{\mathit{n}_{\rm{lag}}}$=%d'%lag,fontsize=35)
        
        if jj != (len(block)-1):
            ax.axes.yaxis.set_ticklabels([])
           
        
        if lag != lags[-1]:
            ax.axes.xaxis.set_ticklabels([])
        

add = fig.add_subplot(grid[:-4,:], frameon=False)
add.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
add.minorticks_off()
add.xaxis.set_label_position("bottom")
add.set_xlabel('phase/$^\circ$',fontsize=35)
add.yaxis.set_label_position("right")
add.set_ylabel('phase/$^\circ$',fontsize=35)
add.yaxis.set_label_coords(1.08, 0.5)
            
axis = np.array(axis)
align_yaxis(axis)

ca = fig.add_subplot(grid[-1,:], frameon=False)
pos = ca.get_position()
pos2 = [pos.x0, pos.y0-0.015, pos.width, pos.height]
ca.set_position(pos2)
plt.sca(ca)
plt.colorbar(corr_map,cax=ca,orientation='horizontal')

al = fig.add_subplot(111,frameon=False)
al.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
al.minorticks_off()
#al.set_title('PSR %s %s'%(star_name,freq),fontsize=35,
#             ha='center',va='bottom',position=(0.5,1.035))

plt.savefig('cor_km_ntsc.pdf',
            bbox_inches='tight',pad_inches = 0)

plt.show()
  

    
