# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 15:47:15 2020

@author: Littl
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.rcParams['font.sans-serif']=['Times New Roman'] #设置字体为罗马字体
plt.rcParams['xtick.direction'] = 'in' #坐标轴刻度向内
plt.rcParams['ytick.direction'] = 'in'
#plt.rcParams['savefig.dpi'] = 500 #保存图片的分辨率
plt.rcParams['font.size']=20


def plot_profile(data,peak_end,w,star_name,freq,bin_num = 1024):
    
    phase = np.linspace(0,360,bin_num)
    phase_peak = phase[peak_end[0]:peak_end[1]]
    
    
    I_peak = w[:,None]*data[:,0,peak_end[0]:peak_end[1]]
    I_prof = I_peak.mean(0)
    
    I_peak = I_peak/max(I_prof)
    I_prof = I_prof/max(I_prof)
    #Q_peak = data[:,1,peak_end[0]:peak_end[1]].mean(0)#(w[:,None]*data[:,1,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
    #U_peak = data[:,2,peak_end[0]:peak_end[1]].mean(0)#(w[:,None]*data[:,2,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
    #L_prof = np.sqrt(Q_peak**2+U_peak**2)
    
    #Q_base = data[:,1,:].mean(0)#(w[:,None]*data[:,1,:]).sum(0)/w.sum(0)
    #U_base = data[:,2,:].mean(0)#(w[:,None]*data[:,2,:]).sum(0)/w.sum(0)
    #L = np.sqrt(Q_base**2+U_base**2)
    #L_base = ((((L[peak_end[0]-50:peak_end[0]])).sum(-1) + ((L[peak_end[1]:peak_end[1]+50])).sum(-1))/
    #                    (L[peak_end[0]-50:peak_end[0]].shape[-1]+L[peak_end[1]:peak_end[1]+50].shape[-1]))
    #L_prof = L_prof - L_base 
    
    #V_prof = (w[:,None]*data[:,3,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
   
    
    
    plt.rcParams['figure.figsize']=(5,12) #设置图像的 宽长比
    plt.figure()
    
    grid=plt.GridSpec(4,1,wspace=0,hspace=0)
    ax1=plt.subplot(grid[0:1,0])
    ax2=plt.subplot(grid[1:4,0]) #上方第一个子图画平均脉冲，下方第二个子图画单个脉冲，两个图所占的高度比为1:3
   
    plt.sca(ax1) #选择第一个子图，进行后续设置
    plt.plot(phase_peak,I_prof,label='I',c='black',linewidth=1)
    #plt.plot(phase_peak,L_prof,label='L',c='orange',linewidth=1)
    #plt.plot(phase_peak,V_prof,label='V',c='blue',linewidth=1)
    plt.legend(frameon=False)
    plt.ylabel('normalised flux')
    plt.xlim(phase_peak[0],phase_peak[-1])
    plt.tick_params(labelbottom=False)
    plt.title('%s %s'%(star_name,freq))
    
    plt.sca(ax2)
    norm = mpl.colors.Normalize(vmin=0, vmax=2) #如果画出来的图 颜色映射不合适的话，可以手动调整。参考值：vmax取平均脉冲轮廓最大值的4倍
    plt.imshow(w[:,None]*I_peak,norm=norm,cmap='binary',aspect='auto',origin='lower',extent=(phase_peak[0],phase_peak[-1], 0, I_peak.shape[0]))
    plt.xlabel('phase/$^\circ$')
    plt.ylabel('pulse number')
    
    
    plt.subplots_adjust(wspace=0,hspace=0) #使得两个子图共享横坐标
    
    plt.savefig('%s_%s_profile.pdf'%(star_name,freq),bbox_inches = 'tight')
    plt.show()

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

plot_profile(data,peak_end,w,star_name='J1644-4559',freq='1400 MHz',bin_num = 4096)

