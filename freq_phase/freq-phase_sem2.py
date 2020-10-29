# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 15:34:48 2020

@author: Littl
"""

import astropy.io.fits as fits
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

plt.rcParams['figure.figsize']=(12,5)

pulse_end = { 'J0738-4042':[420,625],'J0837-4135':[655,730], 'J0942-5552':[330,480],
             'J1327-6222':[65,130],'J1644-4559':[750,840],'J1401-6357':[90,135],
                                  'J1752-2806':[710,765],
            'J1901-0906':[710,880],'J1456-6843':[520,720],
            'J1651-4246':[400,750]
        }

band = [ '704-1728 MHz', '1728-3008 MHz' '3007-4032 MHz']
band = [1728,3008]
import os
file_list = os.listdir(r'./calibPrimit/semester2')

#%%
for file in file_list:
    
    pulsar_name = file.split('.')[0].split('_')[0]
    print(pulsar_name)
    hdul = fits.open(r'./calibPrimit/semester2/%s'%file)
    
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
    base_i = ((I[0,:,:peak_end[0]].sum(-1)+
              I[0,:,peak_end[1]:].sum(-1))/
              (I[0,:,:peak_end[0]].shape[-1]+
               I[0,:,peak_end[1]:].shape[-1]))[:,None]
    i = i - base_i
    #i[np.where(i.mean(-1)==0)] = np.nan
    freq = hdul['SUBINT'].data['DAT_FREQ'][0]
    phase = np.linspace(0,360,1024)
    phase = phase[peak_end[0]:peak_end[1]]
    
    if freq[0] > freq[1]:
        og = 'upper'
    else:
        og = 'lower'
    
    current_cmap = mpl.cm.binary
    current_cmap.set_bad(color='white')
    
    f_b = [min([freq[0],freq[-1]]),max([freq[0],freq[-1]])]
    plt.imshow(i,cmap=current_cmap,aspect='auto',origin=og,extent=(phase[0],phase[-1], f_b[0], f_b[-1]))
    for bound in band:
        plt.axhline(bound,c='k',linestyle='dashed')
    plt.xlabel('phase/$^\circ$')
    plt.ylabel('frequency/MHz')
    plt.title('PSR %s'%file.split('.')[0],y=1.035,fontsize=25)
    plt.savefig('%s.pdf'%file.split('.')[0],
                bbox_inches='tight',pad_inches = 0)
    plt.show()