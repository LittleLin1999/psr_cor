# -*- coding: utf-8 -*-
"""
Created on Mon Mar  9 15:43:53 2020

@author: Littl
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os

plt.rcParams['font.sans-serif']=['Times New Roman'] #设置字体为罗马字体
plt.rcParams['xtick.direction'] = 'in' #坐标轴刻度向内
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['savefig.dpi'] = 500 #保存图片的分辨率

starname = 'J1401-6357'
filelist = os.listdir('.\\%s'%starname)

band1 = filelist[:3]
band1.append(filelist[-1]) 
band2 = filelist[3:8]
band3 = filelist[8:-1]
#%%
def plot_profile(data,peak_end,w,star_name,freq,bin_num = 1024):
    
    phase = np.linspace(0,360,bin_num)
    phase_peak = phase[peak_end[0]:peak_end[1]]
    
    
    I_peak = data[:,0,peak_end[0]:peak_end[1]]
    I_prof = I_peak.sum(0)/w.sum(0)
    
    Q_peak = data[:,1,peak_end[0]:peak_end[1]].sum(0)/w.sum(0)
    U_peak = data[:,2,peak_end[0]:peak_end[1]].sum(0)/w.sum(0)
    L_prof = np.sqrt(Q_peak**2+U_peak**2)
    
    V_prof = data[:,3,peak_end[0]:peak_end[1]].sum(0)/w.sum(0)
   
    
    
    plt.rcParams['figure.figsize']=(5,12) #设置图像的 宽长比
    plt.figure()
    
    grid=plt.GridSpec(4,1,wspace=0,hspace=0)
    ax1=plt.subplot(grid[0:1,0])
    ax2=plt.subplot(grid[1:4,0]) #上方第一个子图画平均脉冲，下方第二个子图画单个脉冲，两个图所占的高度比为1:3
   
    plt.sca(ax1) #选择第一个子图，进行后续设置
    plt.plot(phase_peak,I_prof,label='I',c='black',linewidth=1)
    plt.plot(phase_peak,L_prof,label='L',c='orange',linewidth=1)
    plt.plot(phase_peak,V_prof,label='V',c='blue',linewidth=1)
    plt.legend()
    plt.ylabel('normalised flux')
    plt.xlim(phase_peak[0],phase_peak[-1])
    plt.tick_params(labelbottom=False)
    plt.title('%s %s'%(star_name,freq))
    
    plt.sca(ax2)
    norm = mpl.colors.Normalize(vmin=0, vmax=4*max(I_prof)) #如果画出来的图 颜色映射不合适的话，可以手动调整。参考值：vmax取平均脉冲轮廓最大值的4倍
    plt.imshow(w[:,None]*I_peak,norm=norm,cmap='jet',aspect='auto',origin='lower',extent=(phase_peak[0],phase_peak[-1], 0, I_peak.shape[0]))
    plt.xlabel('phase(°)')
    plt.ylabel('pulse number')
    
    
    plt.subplots_adjust(wspace=0,hspace=0) #使得两个子图共享横坐标

    plt.figure(1)

#%%

merge = []
for i in band1:
    
    a = np.load('.\\%s\\%s'%(starname,i),allow_pickle=True)
    
    merge.append(a)
    
merge = np.mean(merge,axis=0)

np.save('%s_704-1727MHz'%(starname),merge)

plot_profile(merge,[300,800],np.ones((merge.shape[0],)),starname,'n',bin_num = 1024)

merge = []
for i in band2:
    
    a = np.load('.\\%s\\%s'%(starname,i))
    merge.append(a)
    
merge = np.mean(merge,axis=0)
np.save('%s_1727-3007MHz'%(starname),merge)



merge = []
for i in band3:
    
    a = np.load('.\\%s\\%s'%(starname,i))
    merge.append(a)
    
merge = np.mean(merge,axis=0)
np.save('%s_3007-4031MHz'%(starname),merge)


#plot_profile(merge,[750,840],np.ones((merge.shape[0],)),starname,'n',bin_num = 1024)