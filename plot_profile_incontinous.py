# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 09:02:36 2020

@author: Littl
"""
import os

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

folder_name = ['J1825-0935']
pulse_end = { 'J1825-0935':[[165,235],[610,750]]}

mannual =  { 'J1825-0935':[2330,2331,2333,2335,2447,2448,3284]}

abc = ['a','b','c']
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
plt.rcParams["ytick.minor.visible"] =  True
def align_yaxis(axes):
      
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
    
    
    file_list = os.listdir(r'.\\'+psr_name+'_v3')
    
    #file_list[0], file_list[1], file_list[-1] = file_list[-1], file_list[0], file_list[1]
    file_list = file_list[::-1]
    
    
    peak_end = pulse_end[psr_name]
    
    phase = np.linspace(0,360,1024)
    
    
    fig = plt.figure(figsize=(4*3+2,12))

    grid=plt.GridSpec(4,9,fig,wspace=0.1,hspace=0)
    ax = []
    for j in range(len(file_list)):
        axis=[[plt.subplot(grid[0,j*3]),plt.subplot(grid[1,j*3]),plt.subplot(grid[2:4,j*3])], #上方第一个子图画平均脉冲，下方第二个子图画单个脉冲，两个图所占的高度比为1:3
             [plt.subplot(grid[0,j*3+1:(j+1)*3]),plt.subplot(grid[1,j*3+1:(j+1)*3]),plt.subplot(grid[2:4,j*3+1:(j+1)*3])]] 
        ax.append(axis)
    ax = np.array(ax)
    for index,filename in enumerate(file_list):
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
            sub = 0
        if channel[0] > 1000 and channel[0] < 2000:
            sub = 1
        if channel[0] > 2500:
            sub = 2
        data = np.load('.\%s_v3\\%s'%(psr_name,filename))
        
        
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
    
        data[np.where(w==0)] = np.nan
            
        I_prof = np.nanmean(data[:,0],axis=0)
        top = max(I_prof)
        I_prof = I_prof/top
        
        Q_prof = np.nanmean(data[:,1,:],axis=0)#(w[:,None]*data[:,1,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
        U_prof = np.nanmean(data[:,2,:],axis=0)#(w[:,None]*data[:,2,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
        L_prof = np.sqrt(Q_prof**2+U_prof**2)
        
    
        Q_base = np.nanmean(data[:,1,:],axis=0)#(w[:,None]*data[:,1,:]).sum(0)/w.sum(0)
        U_base = np.nanmean(data[:,2,:],axis=0)#(w[:,None]*data[:,2,:]).sum(0)/w.sum(0)
        L = np.sqrt(Q_base**2+U_base**2)
        L_base = np.nanmean(np.concatenate([L[:peak_end[0][0]],
                            L[peak_end[0][1]:peak_end[1][0]],
                            L[peak_end[1][1]:]],
                            axis=-1),axis=0)
        
        L_prof = L_prof - L_base 
        L_prof = L_prof/top
    
        V_prof = np.nanmean(data[:,3,:],axis=0)
        V_prof = V_prof/top
    

        u = np.nanmean(data[:,1,:],axis=0)#(w[:,None]*data[:,1,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
        q = np.nanmean(data[:,2,:],axis=0)#(w[:,None]*data[:,2,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
        
        ang=np.arctan2(u, q)/np.pi*180./2.
        ang[ang<0] += 180
        
        u_err = np.nanstd(np.concatenate([u[:peak_end[0][0]],
                            u[peak_end[0][1]:peak_end[1][0]],
                            u[peak_end[1][1]:]],
                            axis=-1))
        q_err = np.nanstd(np.concatenate([q[:peak_end[0][0]],
                            q[peak_end[0][1]:peak_end[1][0]],
                            q[peak_end[1][1]:]],
                            axis=-1))
        
        ang_err2 = ((u*q_err)**2 + (q*u_err)**2)/(q**2 + u**2)**2
        ang_err = np.sqrt(ang_err2)/np.pi*180/2
     
        stdL = np.nanstd(np.concatenate([L[:peak_end[0][0]],
                            L[peak_end[0][1]:peak_end[1][0]],
                            L[peak_end[1][1]:]],
                            axis=-1))
        snrL = np.abs(L/stdL)

        
        for side,end in enumerate(peak_end[::-1]):
            
            side = abs(side-1)
            
            
            phase_peak = phase[end[0]:end[1]]


            I_sub = I_prof[end[0]:end[1]]
            L_sub = L_prof[end[0]:end[1]]
            V_sub = V_prof[end[0]:end[1]]
            
            
            ang_sub = ang[end[0]:end[1]]       
            ang_err_sub = ang_err[end[0]:end[1]]
            snrL_sub = snrL[end[0]:end[1]]
            idx = np.where(snrL_sub>3)
            
        
            
            ax[sub][side][0].errorbar(phase_peak[idx],ang_sub[idx],ang_err_sub[idx],c='k',fmt='.',linewidth=0.2,ms=1.5)
            ax[sub][side][0].set_ylim(max(min(ang)-5,-180),min(max(ang)+5,180))
        
            ax[sub][side][0].set_yticks([30,60,90,120,150,180])
            ax[sub][side][0].locator_params(axis='y',nbins=3)

        
            ax[sub][side][0].set_xlim(phase_peak[0],phase_peak[-1])
            #ax[sub][side][0].set_title("%d-%d MHz"%(channel[0],channel[1]),fontsize=33)
    

            plt.sca(ax[sub][side][1]) #选择一个子图，进行后续设置
            plt.plot(phase_peak,I_sub,label='I',c='black',linewidth=1)
            plt.plot(phase_peak,L_sub,label='L',c='red',linewidth=1,linestyle='dashed')
            plt.plot(phase_peak,V_sub,label='V',c='blue',linewidth=1,linestyle='dotted')
            plt.yticks([0,0.5,1])
            #plt.locator_params(axis='x',nbins=3)
            plt.xlim(phase_peak[0],phase_peak[-1])
            plt.tick_params(labelbottom=False)
            
            if sub==0 and side==0:
                plt.ylabel('normalised flux',fontsize=25)
            
            else:
                
                plt.yticks([])
               # if side == 1:
               #     plt.title('%d - %d MHz'%(channel[0],channel[1]))

            plt.subplots_adjust(wspace=0.1)
            
            plt.sca(ax[sub][side][2])
            norm = mpl.colors.Normalize(vmin=0, vmax=4*top) #如果画出来的图 颜色映射不合适的话，可以手动调整。参考值：vmax取平均脉冲轮廓最大值的4倍
            current_cmap = mpl.cm.binary
            current_cmap.set_bad(color='white')
            I_peak = data[:,0,end[0]:end[1]]
            plt.imshow(I_peak,norm=norm,cmap=current_cmap,aspect='auto',origin='lower',extent=(phase_peak[0],phase_peak[-1], 0, I_peak.shape[0]))
            
            #plt.xticks(phase_peak[begin:back])
            
            #plt.locator_params(axis='x',nbins=Nbin)
            #plt.gca().xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
            if sub == 0 and side==0:
                #plt.xlabel('phase(°)')
                plt.ylabel('pulse number')
            else:
                plt.yticks([])
    
            #plt.subplots_adjust(wspace=0,hspace=0) #使得两个子图共享横坐标
        
        dash = (0,(26,14,26,14))
        d = .01
        for a,axis in enumerate(ax[sub][0]):
            axis.spines['right'].set_visible(False)

            #axis.spines['right'].set_linewidth(0.5)
            #axis.spines['right'].set_color('green')
            axis.yaxis.tick_left()
            kwargs = dict(transform=axis.transAxes, color='k', clip_on=False)
            if a == 2:
                axis.plot((1-2*d,1+2*d),(-d,+d), **kwargs) # top-left diagonal
                axis.plot((1-2*d,1+2*d),(1-d,1+d), **kwargs) # bottom-left diagonal
            else:
                axis.plot((1-2*d,1+2*d),(-2*d,+2*d), **kwargs) # top-left diagonal
                axis.plot((1-2*d,1+2*d),(1-2*d,1+2*d), **kwargs) # bottom-left diagonal
        for a,axis in enumerate(ax[sub][1]):
            axis.spines['left'].set_visible(False)
            
            #axis.spines['left'].set_linewidth(0.5)
            #axis.spines['left'].set_color('green')
            axis.yaxis.tick_right()
            kwargs.update(transform=axis.transAxes) # switch to the bottom axes
            if a==2:
                axis.plot((-d,d),(-d,+d), **kwargs) # top-right diagonal
                axis.plot((-d,d),(1-d,1+d), **kwargs) # bottom-right diagonal
            else:
                axis.plot((-d,+d),(-2*d,+2*d), **kwargs) # top-left diagonal
                axis.plot((-d,+d),(1-2*d,1+2*d), **kwargs) # bottom-left diagonal
            
        ax[sub][1][0].set_title("%d-%d MHz"%(channel[0],channel[1]),
                                                    x=-0.3,y=1.05,
                                                    horizontalalignment='left',
                                                    fontsize=30) 
    plt.sca(ax[-1][-1][1])
    plt.legend(loc='upper left', bbox_to_anchor=(1,1.05),borderaxespad = 0.02,ncol=1,frameon=False,
               fontsize=25)
       
    ax[0][0][0].set_ylabel(r'PA/$^\circ$')
    
    
    for i in range(3):
        for jj in range(3):
                ax[i][0][jj].set_xticks([70])
                ax[i][1][jj].set_xticks([225,250])
        ax[0][0][i].tick_params(labelleft=True,labelright=False)
        ax[i][0][-1].tick_params(labelbottom=True)
        ax[i][1][-1].tick_params(labelbottom=True)
        

    
    
    align_yaxis(np.array(ax)[:,:,1].flatten())
    align_yaxis(np.array(ax)[:,:,0].flatten())
    
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.minorticks_off()
    plt.title('PSR %s'%psr_name,y=1.1,fontsize=35)
    
    ax[1][-1][-1].set_xlabel('phase/$^\circ$',ha='center',va='top',x=0.25,y=-0.5)
    plt.savefig('.\\revised_v3\\profile_v3\\%s_profile.pdf'%(psr_name),bbox_inches='tight',pad_inches = 0)
    plt.show()