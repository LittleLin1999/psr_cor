# -*- coding: utf-8 -*-
"""
Created on Sun Oct  4 10:38:15 2020

@author: Littl
"""

import astropy.io.fits as fits
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pandas as pd

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

plt.rcParams['figure.figsize']=(12,12)

P = { 'J0738-4042':0.3749,'J0837-4135':0.7516, 'J0942-5552':0.6643,
             'J1327-6222':0.5299,'J1644-4559':0.4551,'J1401-6357':0.8428,
              'J1825-0935':0.7690,  'J1752-2806':0.5626,
            'J1456-6843':0.2634,
            'J1651-4246':0.8441
        }

Pdot = { 'J0738-4042':1.38,'J0837-4135':3.54, 'J0942-5552':22.8,
             'J1327-6222':18.8,'J1644-4559':20.1,'J1401-6357':16.9,
              'J1825-0935':52.4,  'J1752-2806':8.13,
            'J1456-6843':0.0990,
            'J1651-4246':4.76
        }

cor_length = {'J0738-4042':5,'J0837-4135':1,'J1327-6222':12,
              'J1401-6357':4,'J1752-2806':3,'J1644-4559':14}

plt.figure()
for i,l in cor_length.items():
    plt.plot(P[i],l,marker='*',ms=20)
plt.xlabel(r'$P/s$')
plt.ylabel(r'$n_{\rm lag}$',fontsize=30)


plt.figure()
for i,l in cor_length.items():
    plt.plot(Pdot[i],l,marker='*',ms=20)
plt.xlabel(r'$\dot{P}/s$')
plt.ylabel(r'$n_{\rm lag}$',fontsize=30)