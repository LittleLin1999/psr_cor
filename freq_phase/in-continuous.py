# -*- coding: utf-8 -*-
"""
Created on Fri Apr  3 21:45:27 2020

@author: Littl
"""

#%%

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


band = [ '704-1728 MHz', '1728-3008 MHz' '3007-4032 MHz']

pulse_end = [[170,240],[610,750]]
band = [1728,3008]

#%%

file = 'J1825-0935.calibPrmit'

pulsar_name = file.split('.')[0].split('_')[0]
hdul = fits.open(r'./calibPrimit/%s'%file)

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


I = outval_A+outval_B

I[0][np.where(I[0].mean(-1)==0)] = np.nan

pulse_off = np.concatenate([I[0,:,:170],
                            I[0,:,pulse_end[0][1]:pulse_end[1][0]],
                            I[0,:,750:]],
                            axis=-1)
rms = np.std(pulse_off,axis=1)
I = I/rms[:,None]
i = I[0,:,:]
#np.nanmean(pulse_off,axis=-1)[:,None]
pulse_off = np.concatenate([I[0,:,:170],
                            I[0,:,pulse_end[0][1]:pulse_end[1][0]],
                            I[0,:,750:]],
                            axis=-1)
base_i = np.nanmean(pulse_off,axis=-1)[:,None]
#((I[0,:,peak_end[0]-70:peak_end[0]].sum(-1)+
#          I[0,:,peak_end[1]:peak_end[1]+70].sum(-1))/
#          (I[0,:,peak_end[0]-70:peak_end[0]].shape[-1]+
#           I[0,:,peak_end[1]:peak_end[1]+70].shape[-1]))[:,None]
i = i - base_i

freq = hdul['SUBINT'].data['DAT_FREQ'][0]


phase = np.linspace(0,360,1024)
phase = [phase[pulse_end[0][0]:pulse_end[0][1]],phase[pulse_end[1][0]:pulse_end[1][1]]]

if freq[0] > freq[1]:
    og = 'upper'
else:
    og = 'lower'

current_cmap = mpl.cm.binary
current_cmap.set_bad(color='white')

fig = plt.figure(figsize=(12,5))
grid=plt.GridSpec(1,3,fig,wspace=0.05,hspace=0.)
ax=plt.subplot(grid[:,0])
ax2=plt.subplot(grid[:,1:])


ax.imshow(i[:,pulse_end[0][0]:pulse_end[0][1]],cmap=current_cmap,aspect='auto',origin=og,extent=(phase[0][0],phase[0][-1], freq[0], freq[-1]))
ax2.imshow(i[:,pulse_end[1][0]:pulse_end[1][1]],cmap=current_cmap,aspect='auto',origin=og,extent=(phase[1][0],phase[1][-1], freq[0], freq[-1]))

plt.subplots_adjust(wspace=0.05)
ax.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax.yaxis.tick_left()
#ax.tick_params(labeltop='off') # don't put tick labels at the top
ax2.yaxis.tick_right()

d = .01
kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
ax.plot((1-2*d,1+2*d),(-d,+d), **kwargs) # top-left diagonal
ax.plot((1-2*d,1+2*d),(1-d,1+d), **kwargs) # bottom-left diagonal
kwargs.update(transform=ax2.transAxes) # switch to the bottom axes
ax2.plot((-d,d),(-d,+d), **kwargs) # top-right diagonal
ax2.plot((-d,d),(1-d,1+d), **kwargs) # bottom-right diagonal

for bound in band:
    ax.axhline(bound,c='k',linestyle='dashed')
    ax2.axhline(bound,c='k',linestyle='dashed')
    
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