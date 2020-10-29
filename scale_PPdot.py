# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 16:41:05 2020

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

ng = pd.read_csv('./revised_v2/negative_points.csv')

con = np.where(ng['VMIN']<0)
ng = ng.iloc[con]
    
marker = {'I':'o','L':"^",'V':"s"}
color = {0:'k',1:'red',2:'blue'}  #related to freq
ffreq = {0:'low freq',1:'middle freq',2:'high freq'}  #related to freq
#%%
fig1 = plt.figure(figsize=(12,24))
ax = fig1.add_subplot(211, projection='3d')
bx = fig1.add_subplot(212, projection='3d')

X = np.arange(-0.7, 0, 0.05)
Y = np.arange(0, 2, 0.05)
X, Y = np.meshgrid(X, Y)
Z = np.zeros(X.shape)
ax.plot_surface(X, Y, Z, alpha=0.1)
bx.plot_surface(X, Y, Z, alpha=0.1)


for psr_name, p in P.items():
    con = np.where(ng['PSR']==psr_name) 
    data = ng.iloc[con]
    pdot = Pdot[psr_name]
    
    for index, row in data.iterrows():
        cc = (int(row['FREQ'].split('-')[0]))//1300
        if cc==1 and row['INTENSITY']=='I':
            print(p,pdot,abs(row['DIST']))
            ax.bar3d(np.log10(p),np.log10(pdot),0,dx=0.01,dy=0.05,dz=np.log10(abs(row['VMIN'])))
            bx.bar3d(x=np.log10(p),y=np.log10(pdot),z=0,dx=0.01,dy=0.05,dz=np.log10(abs(row['DIST'])))


ax.set_xlabel('$log_{10}(P/s)$',labelpad=20)
ax.set_ylabel('$log_{10}(\dot{P} / \\times 10^{-15})$',labelpad =20)
ax.set_zlabel('$log_{10}(Abs ~of ~the ~most~ negative ~c_{ij}/^\circ $)',labelpad =20)
ax.set_title('values of the most negative $c_{ij}$ for I at the middle frequency band')
ax.set_xlim(X.min(),X.max())
ax.set_ylim(Y.min(),Y.max())


bx.set_xlabel('$log_{10}(P/s)$',labelpad=20)
bx.set_ylabel('$log_{10}(\dot{P} / \\times 10^{-15})$',labelpad =20)
bx.set_zlabel('$log_{10}(seperation~ of~ the~ most~ negative~ c_{ij}/^\circ) $',labelpad =20)
bx.set_title('seperation of the most negative $c_{ij}$ for I at the middle frequency band')
bx.set_xlim(X.min(),X.max())
bx.set_ylim(Y.min(),Y.max())

#plt.savefig('neg_I.png',dpi=200)
plt.show()

#%%
fig1 = plt.figure(figsize=(12,24))
ax = fig1.add_subplot(211, projection='3d')
bx = fig1.add_subplot(212, projection='3d')

X = np.arange(-0.6, 0, 0.05)
Y = np.arange(0, 2, 0.05)
X, Y = np.meshgrid(X, Y)
Z = np.zeros(X.shape)
ax.plot_surface(X, Y, Z, alpha=0.1)
bx.plot_surface(X, Y, Z, alpha=0.1)


for psr_name, p in P.items():
    con = np.where(ng['PSR']==psr_name) 
    data = ng.iloc[con]
    pdot = Pdot[psr_name]
    
    for index, row in data.iterrows():
        cc = (int(row['FREQ'].split('-')[0]))//1300
        if cc==1 and row['INTENSITY']=='L':
            print(p,pdot,abs(row['DIST']))
            ax.bar3d(np.log10(p),np.log10(pdot),0,dx=0.01,dy=0.05,dz=np.log10(abs(row['VMIN'])))
            bx.bar3d(x=np.log10(p),y=np.log10(pdot),z=0,dx=0.01,dy=0.05,dz=np.log10(abs(row['DIST'])))


ax.set_xlabel('$log_{10}(P/s)$',labelpad=20)
ax.set_ylabel('$log_{10}(\dot{P} / \\times 10^{-15})$',labelpad =20)
ax.set_zlabel('$log_{10}(Abs ~of ~the ~most~ negative ~c_{ij}/^\circ $)',labelpad =20)
ax.set_title('values of the most negative $c_{ij}$ for L at the middle frequency band')
ax.set_xlim(X.min(),X.max())
ax.set_ylim(Y.min(),Y.max())


bx.set_xlabel('$log_{10}(P/s)$',labelpad=20)
bx.set_ylabel('$log_{10}(\dot{P} / \\times 10^{-15})$',labelpad =20)
bx.set_zlabel('$log_{10}(seperation~ of~ the~ most~ negative~ c_{ij}/^\circ) $',labelpad =20)
bx.set_title('seperation of the most negative $c_{ij}$ for L at the middle frequency band')
bx.set_xlim(X.min(),X.max())
bx.set_ylim(Y.min(),Y.max())

#plt.savefig('neg_L.png',dpi=200)
plt.show()

#%%
fig1 = plt.figure(figsize=(12,24))
ax = fig1.add_subplot(211, projection='3d')
bx = fig1.add_subplot(212, projection='3d')

X = np.arange(-0.7, 0, 0.05)
Y = np.arange(-1.5, 2, 0.05)
X, Y = np.meshgrid(X, Y)
Z = np.zeros(X.shape)
ax.plot_surface(X, Y, Z, alpha=0.1)
bx.plot_surface(X, Y, Z, alpha=0.1)


for psr_name, p in P.items():
    con = np.where(ng['PSR']==psr_name) 
    data = ng.iloc[con]
    pdot = Pdot[psr_name]
    
    for index, row in data.iterrows():
        cc = (int(row['FREQ'].split('-')[0]))//1300
        if cc==1 and row['INTENSITY']=='V':
            print(p,pdot,abs(row['DIST']))
            ax.bar3d(np.log10(p),np.log10(pdot),0,dx=0.01,dy=0.05,dz=np.log10(abs(row['VMIN'])))
            bx.bar3d(x=np.log10(p),y=np.log10(pdot),z=0,dx=0.01,dy=0.05,dz=np.log10(abs(row['DIST'])))


ax.set_xlabel('$log_{10}(P/s)$',labelpad=20)
ax.set_ylabel('$log_{10}(\dot{P} / \\times 10^{-15})$',labelpad =20)
ax.set_zlabel('$log_{10}(Abs ~of ~the ~most~ negative ~c_{ij}/^\circ $)',labelpad =20)
ax.set_title('values of the most negative $c_{ij}$ for |V| at the middle frequency band')
ax.set_xlim(X.min(),X.max())
ax.set_ylim(Y.min(),Y.max())


bx.set_xlabel('$log_{10}(P/s)$',labelpad=20)
bx.set_ylabel('$log_{10}(\dot{P} / \\times 10^{-15})$',labelpad =20)
bx.set_zlabel('$log_{10}(seperation~ of~ the~ most~ negative~ c_{ij}/^\circ) $',labelpad =20)
bx.set_title('seperation of the most negative $c_{ij}$ for |V| at the middle frequency band')
bx.set_xlim(X.min(),X.max())
bx.set_ylim(Y.min(),Y.max())

#plt.savefig('neg_V.png',dpi=200)
plt.show()

#%%
dc = pd.read_csv('./revised_v2/decoherence.csv')
X = np.arange(-0.7, 0, 0.05)
Y = np.arange(-1.5, 2, 0.05)
X, Y = np.meshgrid(X, Y)
Z = np.zeros(X.shape)

fig1 = plt.figure(figsize=(12,36))
ax = fig1.add_subplot(311)#, projection='3d')
bx = fig1.add_subplot(312)#, projection='3d')
cx = fig1.add_subplot(313)#, projection='3d')
'''
ax.plot_surface(X, Y, Z, alpha=0.1,shade=True)
bx.plot_surface(X, Y, Z, alpha=0.1,shade=True)
cx.plot_surface(X, Y, Z, alpha=0.1,shade=True)
    '''
for psr_name, p in P.items():
    con = np.where(dc['PSR']==psr_name)   
    data = dc.iloc[con]
    pdot = Pdot[psr_name]

    
    for index, row in data.iterrows():
        cc = (int(row['FREQ'].split('-')[0]))//1300
        if cc==1:
            
            if row['INTENSITY']=='I':
                ax.scatter(p/2/pdot,row['LENGTH'])
                #ax.bar3d(x=np.log10(p),y=np.log10(pdot),z=0,dx=0.01,dy=0.05,dz=np.log10(row['LENGTH']),zorder=2)
                
            if row['INTENSITY']=='L':
                bx.plot(p/2/pdot,row['LENGTH'])
                #bx.bar3d(x=np.log10(p),y=np.log10(pdot),z=0,dx=0.01,dy=0.05,dz=np.log10(row['LENGTH']))
                
            if row['INTENSITY']=='V':
                cx.plot(p/2/pdot,row['LENGTH'])
                #cx.bar3d(x=np.log10(p),y=np.log10(pdot),z=0,dx=0.01,dy=0.05,dz=np.log10(row['LENGTH']))
                
'''
ax.set_xlabel('$log_{10}(P/s)$',labelpad=20)
ax.set_ylabel('$log_{10}(\dot{P} / \\times 10^{-15})$',labelpad =20)
ax.set_zlabel('$log_{10}(Average~ decoherence~ length/^\circ)$ ',labelpad =20)
ax.set_title('average decoherence length for I at the middle frequency band')
ax.set_xlim(X.min(),X.max())
ax.set_ylim(Y.min(),Y.max())

bx.set_xlabel('$log_{10}(P/s)$',labelpad=20)
bx.set_ylabel('$log_{10}(\dot{P} / \\times 10^{-15})$',labelpad =20)
bx.set_zlabel('$log_{10}(Average~decoherence~length/^\circ$) ',labelpad =20)
bx.set_title('average decoherence length for L at the middle frequency band')
bx.set_xlim(X.min(),X.max())
bx.set_ylim(Y.min(),Y.max())

cx.set_xlabel('$log_{10}(P/s)$',labelpad=20)
cx.set_ylabel('$log_{10}(\dot{P} / \\times 10^{-15})$',labelpad =20)
cx.set_zlabel('$log_{10}(Average ~decoherence~ length/^\circ$) ',labelpad =20)
cx.set_title('average decoherence length for |V| at the middle frequency band')
cx.set_xlim(X.min(),X.max())
cx.set_ylim(Y.min(),Y.max())
'''

#plt.savefig('decoherence.png',dpi=200)
plt.show()

