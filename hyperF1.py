# -*- coding: utf-8 -*-
"""
Created on Tue Mar 10 10:04:08 2020

@author: Littl
"""
plt.rcParams['font.sans-serif']=['Times New Roman'] #设置字体为罗马字体
plt.rcParams['xtick.direction'] = 'in' #坐标轴刻度向内
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['savefig.dpi'] = 500 #保存图片的分辨率

import scipy.special as sc
import numpy as np
import mpmath as mm

def distribution(n,r,rho):
    
    par00 = np.log(n-2)+(n/2-2)*np.log(1-r**2)+(n/2-1/2)*np.log(1-rho**2)
    par01 = 1/2*np.log(2*np.pi)-(n-3/2)*np.log(1-rho*r)
    par0 = (par00-par01)
    
    par10 = sc.gammaln(n-1)
    par11 = sc.gammaln(n-1/2)
    par1 = (par10-par11)
    
    par2 = sc.hyp2f1(1/2,1/2,n-1/2,1/2+1/2*rho*r)
    
    print(par0,par1,par2)
    return np.e**((par0+par1)+np.log(par2))

def f(n,r):
    par0 = n*np.log(2)+(n/2-1.5)*np.log(1-r**2)+2*sc.gammaln(n/2+1)
    par1 = np.log(np.pi)+2*np.log(n)+sc.gammaln(n-1)
    return np.e**(par0-par1)

def sigma(n,rho):
    
    inner0 = n/2*mm.log(1-rho**2)
    
    inner1 = mm.log(mm.hyp3f2(3/2,n/2,n/2,1/2,n/2+1,rho**2))
  
    inner2 = -mm.log(n)
    inner = mm.mp.e**(inner0+inner1+inner2)
    
    outer0 = n/2*np.log(1-rho**2)+np.log(abs(rho))
    outer1 = 2*sc.gammaln(n/2+1/2)
    outer2 = mm.log(mm.hyp2f1(n/2+1/2,n/2+1/2,n/2+1,rho**2))
    outer3 = -sc.gammaln(n/2)
    outer4 = -sc.gammaln(n/2+1)
    outer = outer0+outer1+outer2+outer3+outer4
    outer = 2*outer
    
    return float(mm.sqrt(inner - mm.mp.e**(outer)))


import os

folder_name = ['J0738-4042','J0837-4135', 'J0942-5552',
               'J1327-6222','J1401-6357','J1644-4559',
                           'J1752-2806','J1456-6843',
               'J1901-0906',
               'J1651-4246']
pulse_end = { 'J0738-4042':[590,800],'J0837-4135':[655,730], 'J0942-5552':[330,480],
             'J1327-6222':[65,130],'J1644-4559':[750,840],'J1401-6357':[80,135],
                                  'J1752-2806':[710,765],
            'J1901-0906':[710,880],'J1456-6843':[520,780],
            'J1651-4246':[400,750]
        }
mannual =  { 
        'J0738-4042':[],'J0837-4135':[2339],
        'J0942-5552':[204,206,800,801,802,803,1201,1550,1551,1628,1629,1634,1635,1663]+list(range(5100,5200)),
        'J1327-6222':[1866,2327,2328],'J1401-6357':[],'J1644-4559':[],
        'J1651-4246':[],'J1752-2806':[],'J1825-0935':[2330,2331,2333,2335,2447,2448,3284],
        'J1901-0906':[],'J1456-6843':[]
        
        
        }

import matplotlib.pyplot as plt
for psr_name in folder_name:
    
    file_list = os.listdir('./'+psr_name+'_v2')
    
    #os.chdir('plot_error')
    

    data = np.load('.\%s_v2\%s'%(psr_name,file_list[0]))
    peak_end = pulse_end[psr_name]
    base_data = ((data[:,:,:peak_end[0]].sum(-1) + data[:,:,peak_end[1]:].sum(-1))/
                    (data[:,:,:peak_end[0]].shape[-1]+data[:,:,peak_end[1]:].shape[-1]))
    data = data - base_data[:,:,None]
    left=peak_end[0]
    right=peak_end[1]

    
    w = np.ones((data.shape[0],))
    mean = data[:,0].mean(1)
    distb = np.where(mean<=0)
    w[distb] = 0
    w[mannual[psr_name]] = 0

    N = np.sum(w)
    
    s = np.linspace(0,0.6,100)
    if psr_name == 'J1456-6843':
        s = np.linspace(0,0.5,100)
  
    sig = []
    for i in s:
        sig.append(sigma(int(N),i))
    
    plt.figure()
    plt.plot(s,s)
    plt.fill_between(s,s+sig,s-sig,alpha=0.3)
    plt.title('The distribution of the cross-correlation coefficient for %s'%psr_name)
    plt.xlabel('Cross-correlation coefficient')
    plt.ylabel('Cross-correlation coefficient')
    plt.savefig('./plot_error/error_%s'%psr_name)
    plt.show()
    
    plt.figure()
    plt.fill_between(s,100*np.array(sig/s),-100*np.array(sig/s),alpha=0.3)
    plt.title('The distribution of the cross-correlation coefficient for %s'%psr_name)
    plt.xlabel('Cross-correlation coefficient')
    plt.ylabel('$(\Delta{r}/r)/\%$')
    plt.savefig('./plot_error/frac_err_%s'%psr_name)
    plt.show()