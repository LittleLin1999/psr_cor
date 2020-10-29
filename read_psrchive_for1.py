import numpy as np
import psrchive
from astropy.io import fits
#import matplotlib.pyplot as plt
import gc
freq_bin = ['low','middle','high']
file_num = {
'J0837-4135':32,
'J1327-6222':74,
'J1752-2806':70,
'J1651-4246':40
}


pulseon_dict = {

'J0837-4135':[[600,800]],
'J1327-6222':[[30,200]],
'1752-2806':[[685,790]],
'J1651-4246':[[300,800]]
}



dirct = '/home/jjc/parkes/'
def get_file_list1(psr_name,freq_file,num):

    file_list = []
    path = '%s%s'%(dirct,psr_name)
    

    file_path = '%s/%s_calibrated'%(path,freq_file)

        
    for a in np.arange(num+1):
        b="%04d"%a
    
        file = "%s/%s_%s_%s.calibPrm"%(file_path,psr_name,freq_file, b)
        file_list.append(file)
            
    return file_list            


def read_single_file(filename, slice_bl,freq_info,freq_bin):
    print('open file %s'%filename)
    arch = psrchive.Archive_load(filename)
    arch.remove_baseline() # remove baseline
    arch.convert_state('Stokes')
    data = arch.get_data() # (subint, pol, chan, bin)
    #nchan = data.shape[2]
    #print(data.shape)
    #freq = arch.get_frequencies()
    wts = arch.get_weights()
    zap_channel = np.all(wts  == 0, axis=0) # channels zapped in the ar file
    zap_subint = np.all(wts == 0, axis=1) # sub-integrations zapped in the ar file
    data[:,:,zap_channel] = np.nan # set the zapped channels to NaN
   # print(data.shape)
    noise = data[:, :, :, slice_bl] # get the pulse off data
    var_file = noise[~zap_subint].var(axis=(0,3)) # RMS of the Stokes vector (pol, chan)
    rms_file = np.sqrt(var_file)
    inorm = rms_file[0]
    data /= inorm[None, None, :, None]
    #profile_file = np.nanmean(data, axis=2)
    profile_file = []
    for k in range(len(freq_bin)):
        kk = freq_bin[k]
        #print('%d-%dMHz'%(freq_info[k][0],freq_info[k][1]))
        profile_file.append(np.nanmean(data[:, :,kk[0]:kk[1], :], axis=2))
    np.save('%s.npy'%(filename.split('/')[-1].split('.')[0]), profile_file)
       
    del profile_file
    gc.collect()
        
   # return(var_file, np.count_nonzero(~zap_subint)) # profiel_file.shape=(nchn,subint, pol, bin), var_file.shape=(pol, chan)


def read_files(filename_list,pulsarname,freq_info,freq_bin):
    


    # get the pulse off region
    slice_bl = np.repeat(True, nbin)
    
    
    for r in pulseon_dict[pulsarname]:
        slice_bl[r[0]-30:r[1]+50] = False
        
        '''
        for i in range(r[0] - 100, r[1] + 100):
            if i <= (nbin - 1):
                j = i
            else:
                j = i - nbin
            slice_bl[j] = False
        '''    
 #   print(slice_bl.shape)

   # profile = []
   # var = []
   # nsub_valid = []
    
    for filename in filename_list:
        
        read_single_file(filename, slice_bl,freq_info,freq_bin)
       # var_file, nsub_valid_file = read_single_file(filename, slice_bl,freq_info,freq_bin)
       # profile.append(profile_file)
        
      # var.append(var_file)
      # nsub_valid.append(nsub_valid_file)
        
       # del profile_file
      #  del var_file
      #  del nsub_valid_file
       # gc.collect()
        
   # var = np.array(var) # (file, pol, chan)
   # profile = np.concatenate(profile, axis=1) # (chn, subint, pol, bin)
    #nsub_valid = np.array(nsub_valid)
    #var_mean = ((var.T*nsub_valid).sum(2)/nsub_valid.sum()).T # (pol, chan)
    #rms_mean = np.sqrt(var_mean)
    '''
    for j in range(profile.shape[0]):
        np.save('%s_profile_%d-%dMHz.npy'%(pulsarname,freq_info[j][0],freq_info[j][1]), profile[j])
     '''   
    #np.save('%s_rms.npy'%pulsarname, rms_mean)
    '''
    plt.imshow(profile[:,0], aspect='auto')
    plt.show()
    plt.close()
    '''

if __name__ == '__main__':
   #filename_list = sys.argv[1:]
   for pname in file_num:
       print(pname)
       
       filenum = file_num[pname]
       for ff in freq_bin:
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
           print(freq_info)
          # print(freq_bin)
           hdul.close()
           read_files(filename_list,pname,freq_info,freq_bin)
