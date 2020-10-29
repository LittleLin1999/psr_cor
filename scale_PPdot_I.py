# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 16:48:45 2020

@author: Littl
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 16:41:05 2020

@author: Littl
"""

import astropy.io.fits as fits
import matplotlib.pyplot as plt
import matplotlib as mpl
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

plt.rcParams['figure.figsize']=(12,5)

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


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(211)
labels = []
for psr_name,p in Pdot.items():

    con = np.where(ng['PSR']==psr_name)   
    data = ng.iloc[con]
    
    
    for index, row in data.iterrows():
        cc = (int(row['FREQ'].split('-')[0]))//1300
        if cc == 0:
            line = plt.plot(p,row['VMIN'],marker[row['INTENSITY']],color=color[cc],
                 markersize=10,markerfacecolor='none',label=row['INTENSITY']+' '+ffreq[cc])
        labels.append(line[0].get_label())
labels = set(labels)
plt.xlabel('$\dot{P}$')
plt.ylabel('Minimum $c_{ij}$')
plt.legend(loc='upper left',labels=labels,bbox_to_anchor=(1, 0.5))


ax = fig.add_subplot(212)
labels = [] 
for psr_name,p in P.items():
    
    con = np.where(ng['PSR']==psr_name)   
    data = ng.iloc[con]

    for index, row in data.iterrows():
        cc = (int(row['FREQ'].split('-')[0]))//1300
        if cc == 0:
            line = plt.plot(p,row['VMIN'],marker[row['INTENSITY']],color=color[cc],
                 markersize=10,markerfacecolor='none',label=row['INTENSITY']+' '+ffreq[cc])
        labels.append(line[0].get_label())
labels = set(labels)

plt.xlabel('$P$')
plt.ylabel('Minimum $c_{ij}$')

plt.subplots_adjust(hspace=0.5)
#plt.savefig('./revised_v2/PPdot_relation/Min.pdf',bbox_inches='tight',pad_inches = 0)
plt.show()

#%%
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(211)
for psr_name,p in Pdot.items():
    

    con = np.where(ng['PSR']==psr_name)   
    data = ng.iloc[con]

    for index, row in data.iterrows():
        
        cc = (int(row['FREQ'].split('-')[0]))//1300
        if cc == 0:
            plt.plot(p,row['DIST'],marker[row['INTENSITY']],color=color[cc],
                 markersize=10,markerfacecolor='none',label=row['INTENSITY'])
plt.xlabel('$\dot{P}$')
plt.ylabel('distance/$^\circ$')
plt.legend(loc='upper left',labels=labels,bbox_to_anchor=(1, 0.5))
    


ax = fig.add_subplot(212)  
for psr_name,p in P.items():
    
    con = np.where(ng['PSR']==psr_name)   
    data = ng.iloc[con]

    for index, row in data.iterrows():
        cc = (int(row['FREQ'].split('-')[0]))//1300
        if cc == 0:
            plt.plot(p,row['DIST'],marker[row['INTENSITY']],color=color[cc],
                 markersize=10,markerfacecolor='none',label=row['INTENSITY'])
plt.xlabel('$P$')
plt.ylabel('distance/$^\circ$')

plt.subplots_adjust(hspace=0.5)
#plt.savefig('./revised_v2/PPdot_relation/distance.pdf',bbox_inches='tight',pad_inches = 0)
plt.show()

#%%

dc = pd.read_csv('./revised_v2/decoherence.csv')

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(211)
for psr_name,p in Pdot.items():
    
    con = np.where(dc['PSR']==psr_name)   
    data = dc.iloc[con]

    for index, row in data.iterrows():
        
        cc = (int(row['FREQ'].split('-')[0]))//1300
        if cc == 0:
            plt.plot(p,row['LENGTH'],marker[row['INTENSITY']],color=color[cc],
                     markersize=10,markerfacecolor='none',label=row['INTENSITY'])
plt.xlabel('$\dot{P}$')
plt.ylabel('mean decoherence length/$^\circ$')
plt.legend(loc='upper left',labels=labels,bbox_to_anchor=(1, 0.5))

    



ax = fig.add_subplot(212)  
for psr_name,p in P.items():
    
    con = np.where(dc['PSR']==psr_name)   
    data = dc.iloc[con]

    for index, row in data.iterrows():
        
        cc = (int(row['FREQ'].split('-')[0]))//1300
        if cc == 0:
            plt.plot(p,row['LENGTH'],marker[row['INTENSITY']],color=color[cc],
                 markersize=10,markerfacecolor='none',label=row['INTENSITY'])
plt.xlabel('$P$')
plt.ylabel('mean decoherence length/$^\circ$')

plt.subplots_adjust(hspace=0.5)
#plt.savefig('./revised_v2/PPdot_relation/decoherence.pdf',bbox_inches='tight',pad_inches = 0)
plt.show()
plt.show()