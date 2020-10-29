# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import astropy.io.fits as fits

name='J0837-4135_high_0012.calibPrm' #文件名
hdul=fits.open(name)
#hdul['SUBINT'].data['DATA'].shape:(100,4,3328,1024)
dataval_B=hdul['SUBINT'].data['DATA'][:,1,:,:]#读取BB偏振所有子积分、所有通道、所有bin
dataval_A=hdul['SUBINT'].data['DATA'][:,0,:,:]#读取AA偏振所有子积分、所有通道、所有bin
dataval_CR=hdul['SUBINT'].data['DATA'][:,2,:,:]#读取CR偏振所有子积分、所有通道、所有bin
dataval_CI=hdul['SUBINT'].data['DATA'][:,3,:,:]#读取CI偏振所有子积分、所有通道、所有bin
#读取AA偏振的scale 和 offset
n_chn=dataval_A.shape[1]
scl_A=hdul['SUBINT'].data['DAT_SCL'][:,0:n_chn]
offs_A=hdul['SUBINT'].data['DAT_OFFS'][:,0:n_chn]
#读取BB偏振的scale和offset
scl_B=hdul['SUBINT'].data['DAT_SCL'][:,n_chn:2*n_chn]
offs_B=hdul['SUBINT'].data['DAT_OFFS'][:,n_chn:2*n_chn]
#读取CR偏振的scale和offset
scl_CR=hdul['SUBINT'].data['DAT_SCL'][:,2*n_chn:3*n_chn]
offs_CR=hdul['SUBINT'].data['DAT_OFFS'][:,2*n_chn:3*n_chn]
#读取CI偏振的scale和offset
scl_CI=hdul['SUBINT'].data['DAT_SCL'][:,3*n_chn:4*n_chn]
offs_CI=hdul['SUBINT'].data['DAT_OFFS'][:,3*n_chn:4*n_chn]

wts=hdul['SUBINT'].data['DAT_WTS'][:]

#outval格式（subs，channels，bins）
outval_A=(dataval_A*scl_A[:,:,None]+offs_A[:,:,None])*wts[:,:,None]
outval_B=(dataval_B*scl_B[:,:,None]+offs_B[:,:,None])*wts[:,:,None]
outval_CR=(dataval_CR*scl_CR[:,:,None]+offs_CR[:,:,None])*wts[:,:,None]
outval_CI=(dataval_CI*scl_CI[:,:,None]+offs_CI[:,:,None])*wts[:,:,None]

#stokes par : bin freq
i=((outval_A+outval_B))[:,:,:].mean(1)