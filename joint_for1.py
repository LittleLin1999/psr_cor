# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 14:03:11 2020

@author: Littl
"""

import numpy as np

from astropy.io import fits

freq_band = ['low','middle','high']
file_num = {


'J0837-4135':32
}




dirct = '/home/jjc/parkes/'
def get_file_list1(psr_name,freq_file,num):

    file_list = []
    path = '%s%s'%(dirct,psr_name)
    

    file_path = '%s/%s_calibrated'%(path,freq_file)

    begin = 1    
    for a in np.arange(begin,num+1):
        b="%04d"%a
    
        file = "%s/%s_%s_%s.calibPrm"%(file_path,psr_name,freq_file, b)
        file_list.append(file)
    return file_list            

for pname in file_num:
    filenum = file_num[pname]
    for ff in freq_band:
        filename_list = get_file_list1(pname,ff,filenum)
       
        hdul = fits.open(filename_list[0])
        nbin = hdul['SUBINT'].header['NBIN']
        freq = hdul['SUBINT'].data['DAT_FREQ']
   
        freq_block = 256
        freq_num = freq.shape[1]//freq_block -1
        freq_info = []
        freq_bin = []
        for i in np.arange(0,freq_num):
            freq_info.append([freq[0,freq_block*i],freq[0,freq_block*(i+1)]])
            freq_bin.append([freq_block*(i),freq_block*(i+1)]) # from low to high
            freq_info.append([freq[0,freq_block*(freq_num)],freq[0,-1]])
            freq_bin.append([freq_block*(freq_num),-1])
   
        # data = [[]]*(freq_num+1) 不能用*4：*是指针引用
        data = []
        for filename in filename_list :
            subfile = np.load('%s.npy'%(filename.split('/')[-1].split('.')[0])) # (nchn,subint, pol, bin)
            data.append(subfile)
    
        data = np.concatenate(data, axis=1) # (chn,subint,pol,bin)    
        print(data.shape)
        for j in range(data.shape[0]):
            np.save('%s_profile_%d-%dMHz.npy'% (pname,freq_info[j][0],freq_info[j][1]), data[j])    