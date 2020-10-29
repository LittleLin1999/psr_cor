# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 15:34:48 2020

@author: Littl
"""

import astropy.io.fits as fits
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

plt.rcParams['font.sans-serif']=['Times New Roman'] 
plt.rcParams['xtick.direction'] = 'in' 
plt.rcParams['ytick.direction'] = 'in'
#plt.rcParams['savefig.dpi'] = 500 

pulse_end = { 'J0738-4042':[590,800],'J0837-4135':[655,730], 'J0942-5552':[330,480],
             'J1327-6222':[65,130],'J1644-4559':[750,840],'J1401-6357':[80,135],
                                  'J1752-2806':[710,765],
            'J1901-0906':[710,880],'J1456-6843':[520,720],
            'J1651-4246':[375,750],'J1825-0935':[170,750]
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
plt.rcParams['font.size']=20

plt.rcParams['xtick.bottom'] = True
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.left'] = True
plt.rcParams['ytick.right'] = True

plt.rcParams['xtick.labelbottom'] = True
plt.rcParams['xtick.labeltop'] = False
plt.rcParams['ytick.labelleft'] = True
plt.rcParams['ytick.labelright'] = False

plt.rcParams["xtick.minor.visible"] =  True
plt.rcParams["ytick.minor.visible"] =  True



import os
file_list = os.listdir(r'./calibPrimit/semester1')
g = len(file_list)//3
for k in range(g):
    file_group = file_list[k*3:(k+1)*3]
    file_group[-2], file_group[-1] = file_group[-1], file_group[-2]
    print(file_group)
    freq = []
    intensity = []
    
    for file in file_group:
        pulsar_name = file.split('.')[0].split('_')[0]
        hdul = fits.open(r'./calibPrimit/semester1/%s'%file)
        
        dataval_B=hdul['SUBINT'].data['DATA'][:,1,:,:]#读取BB偏振所有子积分、所有通道、所有bin
        dataval_A=hdul['SUBINT'].data['DATA'][:,0,:,:]#读取AA偏振所有子积分、所有通道、所有bin
        dataval_CR=hdul['SUBINT'].data['DATA'][:,2,:,:]#读取CR偏振所有子积分、所有通道、所有bin
        dataval_CI=hdul['SUBINT'].data['DATA'][:,3,:,:]#读取CI偏振所有子积分、所有通道、所有bin
        
        #读取AA偏振的scale 和 offset
        n_chn=dataval_A.shape[1]
        scl_A=hdul['SUBINT'].data['DAT_SCL'][:,0:n_chn]
        offs_A=hdul['SUBINT'].data['DAT_OFFS'][:,0:n_chn]
        #读取BB偏振的scale和offset
        scl_B=hdul['SUBINT'].data['DAT_SCL'][:,n_chn:2*n_chn]
        offs_B=hdul['SUBINT'].data['DAT_OFFS'][:,n_chn:2*n_chn]
        #读取CR偏振的scale和offset
        scl_CR=hdul['SUBINT'].data['DAT_SCL'][:,2*n_chn:3*n_chn]
        offs_CR=hdul['SUBINT'].data['DAT_OFFS'][:,2*n_chn:3*n_chn]
        #读取CI偏振的scale和offset
        scl_CI=hdul['SUBINT'].data['DAT_SCL'][:,3*n_chn:4*n_chn]
        offs_CI=hdul['SUBINT'].data['DAT_OFFS'][:,3*n_chn:4*n_chn]
        
        wts=hdul['SUBINT'].data['DAT_WTS'][:]
        
        #outval格式（subs，channels，bins）
        outval_A=(dataval_A*scl_A[:,:,None]+offs_A[:,:,None])*wts[:,:,None]
        outval_B=(dataval_B*scl_B[:,:,None]+offs_B[:,:,None])*wts[:,:,None]
        outval_CR=(dataval_CR*scl_CR[:,:,None]+offs_CR[:,:,None])*wts[:,:,None]
        outval_CI=(dataval_CI*scl_CI[:,:,None]+offs_CI[:,:,None])*wts[:,:,None]
        
        peak_end = pulse_end[pulsar_name]
        I = outval_A+outval_B
        pulse_off = np.concatenate([I[0,:,:peak_end[0]],I[0,:,peak_end[1]:]],axis=1)
        rms = np.std(pulse_off,axis=1)
        I = I/rms[:,None]
        i = I[0,:,peak_end[0]:peak_end[1]]
        base_i = ((I[0,:,peak_end[0]-70:peak_end[0]].sum(-1)+
                  I[0,:,peak_end[1]:peak_end[1]+70].sum(-1))/
                  (I[0,:,peak_end[0]-70:peak_end[0]].shape[-1]+
                   I[0,:,peak_end[1]:peak_end[1]+70].shape[-1]))[:,None]
        i = i - base_i
        #i[np.where(wts[0])==0] = np.nan
        f = hdul['SUBINT'].data['DAT_FREQ'][0]
        if f[0] < f[1]:
            f = f[::-1]
            i = i[::-1]
        freq.append(f)
        intensity.append(i)
    
    freq_lm = np.concatenate(freq[1:],axis=0)
    intensity_lm = np.concatenate(intensity[1:],axis=0)
    
    freq_h = freq[0]
    intensity_h = intensity[0]
    
    phase = np.linspace(0,360,1024)
    phase = phase[peak_end[0]:peak_end[1]]
    
    
    current_cmap = mpl.cm.binary
    current_cmap.set_bad(color='white')
    
    
    flm_b = [min([freq_lm[0],freq_lm[-1]]),max([freq_lm[0],freq_lm[-1]])]
    fh_b = [min([freq_h[0],freq_h[-1]]),max([freq_h[0],freq_h[-1]])]
   
    fig = plt.figure(figsize=(12,5))
    grid=plt.GridSpec(3,1,fig,wspace=0,hspace=0.2)
    ax=plt.subplot(grid[0:1,0])
    ax2=plt.subplot(grid[1:,0])
    #fig,(ax,ax2) = plt.subplots(2, 1, sharex=True)
    
    ax.imshow(intensity_h[:,:],cmap=current_cmap,aspect='auto',origin='upper',extent=(phase[0],phase[-1], fh_b[0], fh_b[-1]))
    ax2.imshow(intensity_lm[:,:],cmap=current_cmap,aspect='auto',origin='upper',extent=(phase[0],phase[-1], flm_b[0], flm_b[-1]))

   
    ax.spines['bottom'].set_visible(False)
    ax2.spines['top'].set_visible(False)
    ax.set_xticks([])
    #ax.tick_params(labeltop='off') # don't put tick labels at the top
    ax2.xaxis.tick_bottom()
    
    d = .01
    kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    ax.plot((1-d,1+d),(-2*d,+2*d), **kwargs) # top-left diagonal
    ax.plot((-d,d),(-2*d,+2*d), **kwargs) # top-right diagonal
    kwargs.update(transform=ax2.transAxes) # switch to the bottom axes
    ax2.plot((-d,d),(1-d,1+d), **kwargs) # bottom-right diagonal
    ax2.plot((1-d,1+d),(1-d,1+d), **kwargs) # bottom-left diagonal

    ax2.axhline(1344,c='k',linestyle='dashed')
    ax2.axhline(3008,c='k',linestyle='dashed')
    ax.axhline(3264,c='k',linestyle='dashed')

    al = fig.add_subplot(111,frameon=False)
    al.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    al.minorticks_off()
    al.set_title('PSR %s'%(pulsar_name),fontsize=25,
                 ha='center',va='bottom',x=0.5,y=1.035)
    
    al.set_xlabel('phase/$^\circ$')
    al.set_ylabel('frequency/MHz')
    al.yaxis.set_label_coords(-0.12,0.5)
    plt.savefig('%s.pdf'%pulsar_name,
                bbox_inches='tight',pad_inches = 0)
    plt.show()
    