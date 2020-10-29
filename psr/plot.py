# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 00:40:19 2020

@author: Littl
"""

import matplotlib.pyplot as plt
import matplotlib as mpl



import numpy as np

import sys
sys.path.append('../')

import psr.cor as cor
import psr.drifting as drifting

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

flux = ['I','L','V']

def plot_profile(data,peak_end,star_name,freq,bin_num = 1024):
    
    phase = np.linspace(0,360,bin_num)
    phase_peak = phase[peak_end[0]:peak_end[1]]
    
    
    I_peak = data[:,0,peak_end[0]:peak_end[1]]
    I_prof = I_peak.mean(0)
    
    Q_peak = np.nanmean(data[:,1,peak_end[0]:peak_end[1]],axis=0)#(w[:,None]*data[:,1,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
    U_peak = np.nanmean(data[:,2,peak_end[0]:peak_end[1]],axis=0)#(w[:,None]*data[:,2,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
    L_prof = np.sqrt(Q_peak**2+U_peak**2)
    
    Q_base = np.nanmean(data[:,1,:],axis=0)#(w[:,None]*data[:,1,:]).sum(0)/w.sum(0)
    U_base = np.nanmean(data[:,2,:],axis=0)#(w[:,None]*data[:,2,:]).sum(0)/w.sum(0)
    L = np.sqrt(Q_base**2+U_base**2)
    L_base = np.concatenate([L[peak_end[0]-50:peak_end[0]] , L[peak_end[1]:peak_end[1]+50]])
    L_base = np.nanmean(L_base)
                   
    L_prof = L_prof - L_base 
    
    V_prof = np.nanmean(data[:,3,peak_end[0]:peak_end[1]],axis=0)
   
    
    
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
    plt.imshow(I_peak,norm=norm,cmap='jet',aspect='auto',origin='lower',extent=(phase_peak[0],phase_peak[-1], 0, I_peak.shape[0]))
    plt.xlabel('phase(°)')
    plt.ylabel('pulse number')
    
    
    plt.subplots_adjust(wspace=0,hspace=0) #使得两个子图共享横坐标
    
    plt.savefig('%s_%s_profile.pdf'%(star_name,freq))
    plt.show()
    
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
       
def plot_cor(w,data_all,peak_end,lags,star_name,freq,bin_num=1024):
    
    phase = np.linspace(0,360,bin_num)
    
    
    
    I_prof = np.nanmean(data_all[:,0,peak_end[0]:peak_end[1]],axis=0)#(w[:,None]*data_all[:,0,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
    Q_peak = np.nanmean(data_all[:,1,peak_end[0]:peak_end[1]],axis=0)#(w[:,None]*data[:,1,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
    U_peak = np.nanmean(data_all[:,2,peak_end[0]:peak_end[1]],axis=0)#(w[:,None]*data[:,2,peak_end[0]:peak_end[1]]).sum(0)/w.sum(0)
    L_prof = np.sqrt(Q_peak**2+U_peak**2)
    V_prof = np.nanmean(data_all[:,3,peak_end[0]:peak_end[1]],axis=0)
    
    top = max(I_prof)
    I_prof = I_prof/top
    L_prof = L_prof/top
    V_prof = V_prof/top
    
    
    Q = data_all[:,1]
    U = data_all[:,2]
    L = np.sqrt(Q**2+U**2)
    base_L = np.concatenate((L[:,:peak_end[0]],L[:,peak_end[1]:]),axis=1)
   # print(base_L.shape,base_L)
    base_L = np.nanmean(base_L,axis=1)
   # print(base_L.shape)
    L = L- base_L[:,None]
    L = L/top
    
    I = data_all[:,0]/top
    abs_V = abs(data_all[:,-1]/top)
    
    
    phase_peak = phase[peak_end[0]:peak_end[1]]
    
    columns = {'I':I,'L':L,'V':abs_V}
    ls = {'I':['black','solid'],'L':['red','dashed'],'V':['blue','dotted']}
    profile = {'I':I_prof,'L':L_prof,'V':V_prof}
    
    plt.rcParams['figure.figsize']=(14.3,28.4)
    fig=plt.figure()
    grid=plt.GridSpec(108,14,fig,wspace=0.1,hspace=0.4)
    
    for j in np.arange(6):
        ax =  plt.subplot(grid[8+j*16:8+(j+1)*16,-2:])
        plt.sca(ax)
        plt.plot(I_prof,phase_peak,c=ls['I'][0],linewidth=1,label='I')
        plt.plot(L_prof,phase_peak,c=ls['L'][0],linestyle=ls['L'][1],linewidth=1,label='L')
        plt.plot(V_prof,phase_peak,c=ls['V'][0],linestyle=ls['V'][1],linewidth=1,label='V')
        
        plt.ylim(phase_peak[0],phase_peak[-1])
        plt.locator_params(axis='x',nbins=3)
        plt.locator_params(axis='y',nbins=3)
        
        
        ax.axes.xaxis.set_ticklabels([])
        ax.tick_params(labelleft=False,labelright=True)
        
        if j==0:
            plt.legend(loc='lower left', bbox_to_anchor=(0.05,1.02),borderaxespad = 0.0,ncol=1,frameon=False)
        
    jj = -1
    axis = []
    for key,data in columns.items():
        
        jj += 1
        
        pulse_off1 = data[:,peak_end[0]-50:peak_end[0]]
        pulse_off2 = data[:,peak_end[1]:peak_end[1]+50]
    
        pulse_off = np.concatenate((pulse_off1,pulse_off2),axis=1)
        
        ax1 = plt.subplot(grid[:8,4*jj:(jj+1)*4])
        plt.sca(ax1)
        plt.locator_params(axis='y',nbins=3)
        plt.locator_params(axis='x',nbins=3)
        plt.plot(phase_peak,profile[key]/max(abs(profile[key])),label=key,c=ls[key][0],
                 linestyle=ls[key][1],linewidth=1)
        plt.xlim(phase_peak[0],phase_peak[-1])
        ax1.axes.xaxis.set_ticklabels([])
        
        axis.append(ax1)
        
        title = key
       # if key == "V":
       #     title = r'$|$ V $|$'
        plt.title(title,fontsize=35,ha='center',va='bottom',position=(0.5,1.025))
        if jj == 0:
            ax1.yaxis.set_label_position("left")
            ax1.set_ylabel('normalised flux')
            ax1.tick_params(labelright=False)
        else:
            ax1.axes.yaxis.set_ticklabels([])
        for index,lag in enumerate(lags):
            cr_off = cor.weighed_corr_map(w, pulse_off, lag,-10000)
    
            cr_off = np.concatenate([cr_off[:int(cr_off.shape[0]/3),-int(cr_off.shape[1]/3):],
                                    cr_off[-int(cr_off.shape[0]/3):,:int(cr_off.shape[1]/3)]])
  
            cr_off = np.std(cr_off)
    
            data_peak = data[:,peak_end[0]:peak_end[1]]
            cr = cor.weighed_corr_map(w, data_peak, lag,abs(cr_off))
            
            ax = plt.subplot(grid[8+index*16:24+index*16,jj*4:(jj+1)*4])
    
            plt.sca(ax)
        
            norm=mpl.colors.Normalize(vmin=-1,vmax=1) #将1设为最红，-1为最蓝，中间的0为白
            cmap = mpl.cm.seismic
            cmap.set_bad(color='#404040',alpha=0.15)
            corr_map = plt.imshow(cr,norm=norm, cmap=cmap, aspect='equal', origin='lower',
                                  extent=(phase_peak[0], phase_peak[-1], phase_peak[0], phase_peak[-1]))

            ax.axes.yaxis.set_ticklabels([])
            ax.tick_params(labeltop=False)
            plt.locator_params(axis='x',nbins=3)
            plt.locator_params(axis='y',nbins=3)
            if jj == 0:
                plt.ylabel(r'$\mathrm{\mathit{n}_{\rm{lag}}}$=%d'%lag,fontsize=35)
            
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
    
    ca = fig.add_subplot(grid[-1,:-2], frameon=False)
    pos = ca.get_position()
    pos2 = [pos.x0, pos.y0-0.015, pos.width, pos.height]
    ca.set_position(pos2)
    plt.sca(ca)
    plt.colorbar(corr_map,cax=ca,orientation='horizontal')
    
    al = fig.add_subplot(111,frameon=False)
    al.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    al.minorticks_off()
    al.set_title('PSR %s %s'%(star_name,freq),fontsize=35,
                 ha='center',va='bottom',position=(0.5,1.035))
    
    plt.savefig('%s_%s.pdf'%(star_name,freq.replace(' ','')),
                bbox_inches='tight',pad_inches = 0)
    
    plt.show()
  
    
    
def plot_cor_freq(w1,w2,data1,data2,peak_end,star_name,freq1,freq2,lag,bin_num=1024):
    
    phase = np.linspace(0,360,bin_num)
    phase_peak = phase[peak_end[0]:peak_end[1]]
    
    I1_peak = data1[:,peak_end[0]:peak_end[1]]
    I2_peak = data2[:,peak_end[0]:peak_end[1]]
    
    plt.rcParams['figure.figsize']=(6,5)
    plt.figure()
    
    cmap='seismic' #设置颜色映射为 红白蓝
    norm=mpl.colors.Normalize(vmin=-1,vmax=1) #将1设为最红，-1为最蓝，中间的0为白
    cr = cor.weighed_corr_map_freq(w1,w2,I1_peak,I2_peak,lag)
    
    plt.imshow(cr, norm=norm, cmap=cmap, aspect='auto', origin='lower',
           extent=(phase_peak[0], phase_peak[-1], phase_peak[0], phase_peak[-1]))
    plt.colorbar()
    plt.title('%s %s vs.%s cross-correlation : lag = %d' % (star_name,freq1,freq2,lag))
    plt.xlabel('%s phase  (°)'%freq1 )
    plt.ylabel('%s phase  (°)'%freq2 )
    plt.savefig('%s %s vs %s_corr_lag=%d.png'%(star_name,freq1,freq2,lag))
    plt.show()        
    
def plot_cor_noise(auto,w1,w2,data1,data2,flux,peak_end,pulse_end,star_name,freq,lag,bin_num=1024):
    
    left = peak_end[0]
    right = peak_end[1]
    
    noil = pulse_end[0] - left #!!!设置扣噪声的区域
    noir = pulse_end[1] - left
    
    peak1=data1[:,left:right]
    peak2=data2[:,left:right]
    
    phase=np.linspace(0,360,1024)
    phase_peak=phase[left:right]

    #用于生成去噪声区的分界线
    corll = [(noil + left) * 360 / 1024, (noil + left) * 360 / 1024] #!!!
    corlr = [(noil + left) * 360 / 1024, (noir + left) * 360 / 1024] #!!!
    corrr = [(noir + left) * 360 / 1024, (noir + left) * 360 / 1024] #!!!
    
    #计算每个周期噪声的方差
    noise1 = np.zeros(data1.shape[0]) #!!!
    noise2 = np.zeros(data2.shape[0]) #!!!
    for i in range(data1.shape[0]): #!!!
        noise1[i] = (np.sum((data1[i, left-20:left])**2) + np.sum((data1[i, right:right+20])**2) )/ (60) #!!!
    for i in range(data2.shape[0]): #!!!
        noise2[i] = (np.sum((data2[i, left-20:left])**2) + np.sum((data2[i, right:right+20])**2)) / (60) #!!!
       # noise2[i] = (np.sum((data2[i, 0:left]) * (data2[i, 0:left])) + np.sum((data2[i, right:]) * (data2[i, right:]))) / (
       #         left + 1024 - right) #!!!
    
    if auto == True and lag == 0: #!!!
        flagnoi = 0 #!!!
    else: #!!!
        flagnoi = 1 #!!!
    
    cr, delta=cor.weighed_corr_map_noi(w1, w2, peak1, peak2, lag, noise1, noise2, noil, noir, flagnoi) #!!!
    norm=mpl.colors.Normalize(vmin=1,vmax=-1)
    current_cmap = mpl.cm.seismic
    current_cmap.set_bad(color='yellow',alpha=0.3)
    plt.imshow(cr,norm=norm,cmap=current_cmap, aspect='auto', origin='lower',
               extent=(phase_peak[0], phase_peak[-1], phase_peak[0], phase_peak[-1]))
    plt.plot(corll, corlr, c='g', linestyle='--', linewidth=0.5 ) #!!!扣噪声区的界限
    plt.plot(corlr, corrr, c='g', linestyle='--', linewidth=0.5 ) #!!!
    plt.plot(corrr, corlr, c='g', linestyle='--', linewidth=0.5 ) #!!!
    plt.plot(corlr, corll, c='g', linestyle='--', linewidth=0.5 ) #!!!
    plt.colorbar()
    
    plt.title('%s-%s %s cross-correlation (reduce noise) : lag = %d'%(star_name,freq,flux,lag))
    plt.xlabel('phase (°)')
    plt.ylabel('phase (°)')
    plt.savefig('%s_%s_corr_noise_lag=%d.png'%(star_name,freq,lag))
    plt.show()

    plt.imshow(delta,cmap='jet', aspect='auto', origin='lower',
               extent=(phase_peak[0], phase_peak[-1], phase_peak[0], phase_peak[-1]))
    plt.plot(corll, corlr, c='g', linestyle='--', linewidth=0.5 ) #!!!扣噪声区的界限
    plt.plot(corlr, corrr, c='g', linestyle='--', linewidth=0.5 ) #!!!
    plt.plot(corrr, corlr, c='g', linestyle='--', linewidth=0.5 ) #!!!
    plt.plot(corlr, corll, c='g', linestyle='--', linewidth=0.5 ) #!!!
    plt.colorbar()
    
    #plt.title('%s-%s %s cross-correlation (reduce noise) error : lag = %d'%(star_name,freq,flux,lag))
    plt.xlabel('phase (°)')
    plt.ylabel('phase (°)')
    

    plt.savefig('%s_%s_corr_noise_lag=%d.png'%(star_name,freq,lag))
    plt.show()
        
def plot_LRFS(w,data,peak_end,star_name,freq,bin_num=1024):
    
    
    
    phase = np.linspace(0,360,bin_num)
    phase_peak = phase[peak_end[0]:peak_end[1]]
    
    LRFS,phase_angle,f = drifting.LRFS(w[:,None]*data[:,peak_end[0]:peak_end[1]])
    
    
    
    mmp = 5
    for i in range(mmp):
        spec = LRFS.mean(1)[300:]
        max_i = np.argmax(spec) + 300
        fil = np.concatenate([spec[max_i-10:max_i-5],spec[max_i+5:max_i+10]])
        LRFS[max_i] = np.full(LRFS[max_i].shape,np.mean(fil))
    
    
    
    plt.rcParams['figure.figsize']=(10,12)
    
    
    gs = plt.GridSpec(5,3,wspace=0,hspace=0)
    ax1 = plt.subplot(gs[1:-1,1:])
    ax4 = plt.subplot(gs[0,1:])
    ax3 = plt.subplot(gs[1:-1,0])
    ax2 = plt.subplot(gs[-1,1:])
    
   
    
    spec = LRFS[:503].mean(1)
    ll = []
    fs = []
    for s in np.arange(5,len(spec)-5):
        sp = np.mean(LRFS[:503][s-3:s+3],axis=0)
        fq = np.mean(f[s-3:s+3])
        ll.append(sp)
        fs.append(fq)
    ll = np.array(ll)
    fs = np.array(fs)
    plt.sca(ax1)
    
    plt.imshow(ll[:],interpolation='sinc',cmap='hot',aspect='auto',origin='lower',extent=[phase_peak[0],phase_peak[-1],fs[0],fs[-1]])
    #plt.xlabel('longitude/$^\circ$')
    #plt.yticks([])
    plt.tick_params(labeltop=False,labelbottom=False,labelleft=False,
                    labelright=False)

    plt.sca(ax2)    
    
    plt.plot(phase_peak,phase_angle,'.',color='black',markersize=3)
    #plt.yticks([120,0,-120])
    plt.ylabel('FFT phase/$^\circ$')
    ax2.yaxis.set_label_coords(1.3,0.5)
    plt.xlabel('phase/$^\circ$')
    plt.tick_params(labeltop=False,labelbottom=True,labelleft=False,
                    labelright=True)
    
   
    plt.sca(ax3)    
    #plt.xticks([])
    plt.ylim(fs[0],fs[-1])
    plt.yticks([x*0.025 for x in range(int(0.125//0.025)+2)])
    plt.ylabel('frequency $P_1/P_3$')
  
    ss = ll[:,:].mean(-1)
    plt.plot(ss[:]/max(ss[:]),fs[:],color='red',linewidth=0.7,zorder=0)
    plt.tick_params(labelright=False,labelbottom=False)
    plt.gca().invert_xaxis()
    
    
    print(fs[np.argmax(ss)],1/fs[np.argmax(ss)])
    
    d_ss = ss/max(ss)
    #f_base = np.mean(np.concatenate([d_ss[np.argmax(d_ss)-30:np.argmax(d_ss)-10],d_ss[np.argmax(d_ss)+10:np.argmax(d_ss)+30]]))
    f_std = np.std(np.concatenate([d_ss[np.argmax(d_ss)-30:np.argmax(d_ss)-10],d_ss[np.argmax(d_ss)+10:np.argmax(d_ss)+30]]))    
    bulk = np.concatenate([d_ss[np.argmax(d_ss)-30:np.argmax(d_ss)-1],d_ss[np.argmax(d_ss)+1:np.argmax(d_ss)+30]])
    f_bulk = np.concatenate([fs[np.argmax(d_ss)-30:np.argmax(d_ss)-1],fs[np.argmax(d_ss)+1:np.argmax(d_ss)+30]])
    err_pos = np.where(abs(1-bulk-f_std)<0.1)
    
    err  =max(f_bulk[err_pos[0]]-fs[np.argmax(ss)])
    print(err)
   
    print(1/(fs[np.argmax(ss)]+err)-1/fs[np.argmax(ss)],1/(fs[np.argmax(ss)]-err)-1/fs[np.argmax(ss)])
    
    
    
    plt.sca(ax4)
    I_prof = data[:,peak_end[0]:peak_end[1]].mean(0)
    plt.plot(phase_peak,I_prof/max(I_prof),color='k')
    plt.xlim(phase_peak[0],phase_peak[-1])
    plt.tick_params(labelbottom=False,labelleft=False,labelright=True,labeltop=False)
    plt.yticks([0,0.5,1])
    plt.ylabel('normalised flux')
    ax4.yaxis.set_label_coords(1.2,0.5)
    #plt.ylim(min(I_prof),max(I_prof))
    
    plt.title('%s %s'%(star_name,freq),fontsize=35,y=1.1)
    
    plt.figure(1)
    plt.savefig('%s_%s_LRFS.pdf'%(star_name,freq.replace(' ','')),
                bbox_inches='tight',pad_inches = 0)
    plt.show()
    
def plot_sliding_LRFS(w,data,peak_end,block,star_name,freq,bin_num=1024):
    
    S,f = drifting.sliding_LRFS(data[:,peak_end[0]:peak_end[1]],block)
    
        
    plt.rcParams['figure.figsize']=(12,8)
    plt.figure()
    
    gs = plt.GridSpec(4,5)
    ax1 = plt.subplot(gs[:,1:])
    ax2 = plt.subplot(gs[:,0])
    
    plt.sca(ax1)
    plt.imshow(S[:],aspect='auto',cmap='hot',origin='lower',extent=[0,S.shape[1],0,f[-1]])
    plt.xlabel('block')
    #plt.ylabel('Frequency $P_1/P_3$')
    plt.yticks([])
    plt.title('Time varying LRFS of %s %s'%(star_name,freq))
    
    plt.sca(ax2)    
    plt.xticks([])
    plt.ylim(0,f[-1])
    plt.ylabel('frequency $P_1/P_3$')
    plt.plot(-S.mean(1),f,color='black',linewidth=0.2,zorder=0)

    plt.subplots_adjust(wspace =0, hspace =0)
    
    plt.figure(1)
    plt.savefig('%s_%s_sliding LRFS.png'%(star_name,freq))
    plt.show()
    
    
def plot_all(w,data,peak_end,pulse_end,lags,star_name,freq):
    
    #plot_profile(data,peak_end,star_name,freq)

    #plot_LRFS(w,data[:,0],peak_end,star_name,freq)
    
    #plot_sliding_LRFS(w,data[:,0],peak_end,512,star_name,freq,bin_num=1024)
    
   
    plot_cor(w,data,peak_end,lags,star_name,freq,bin_num=1024)

    

        