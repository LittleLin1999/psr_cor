# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 14:30:51 2020

@author: Littl
"""
import os
import numpy as np

dirct = '/home/jjc/parkes/'

obs1 = ['J0837-4135','J1327-6222','1752-2808','J1651-4246']
freq_bin = ['low','middle','high']

obs2 = ['J0942-5552','J1644-4559','J0738-4042','J1825-0935',
        'J1901-0906','J1401-6357']

def get_file_list1(psr_name,freq_bin_index,num):

    file_list = []
    path = '%s%s'%(dirct,psr_name)
    

    file_path = '%s/%s_calibrated'%(path,freq_bin[freq_bin_index])

        
    for a in np.arange(num):
        b="%04d"%a
    
        file = "%s/%s_%s_%s.calibPrm"%(file_path,psr_name,freq_bin[freq_bin_index], b)
        file_list.append(file)
            
    return file_list            

def get_file_list2(psr_name,num):
    
    file_list = []
    path = '%s%s'%(dirct,psr_name)
    

        
    for a in np.arange(num):
        b="%04d"%a
    
        file = "%s/%s_%s.calibPrm"%(path,psr_name, b)
        file_list.append(file)
            
    return file_list  