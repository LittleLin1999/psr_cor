#画出相关系数图，保存相关系数的数据
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import xlwt
plt.rcParams['font.sans-serif']=['Times New Roman']
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['savefig.dpi'] = 500

#要设置：脉冲存在的bins以及频段（高中低全），星的名称，lag的集合，还有相关图的峰对应横纵坐标大致所在的bin，相关图的颜色映射
left=600
right=800 #J0738-4042
#left=600
#right=750 #J1825-0935
chn_1='all'
chn_2='all'
star='J0738-4042'
lags=[0,1,2,3,9,200]

#互相关
def weighed_corr(w1, w2, l1, l2, lag, noise1, noise2, flagnoi): #!!!两个noise对应两个相位处的噪声序列,flag决定要不要扣噪声
    n = len(l1)

    s1 = l1[lag:]
    s2 = l2[:n - lag]

    ww1 = w1[lag:]
    ww2 = w2[:n - lag]
    w = ww1 * ww2

    y1 = w * s1
    y2 = w * s2

    avg1 = np.sum(y1) / np.sum(w)
    avg2 = np.sum(y2) / np.sum(w)

    cov = np.sum((y1 - avg1) * (y2 - avg2))

    var1 = np.sum((y1 - avg1) ** 2)
    var2 = np.sum((y2 - avg2) ** 2)

    if flagnoi == 1:
        noise1 = noise1[lag:] * w
        noise2 = noise2[:n - lag] * w
        var01 = var1 - np.sum(noise1) * (1 - 1 / np.sum(w))
        var02 = var2 - np.sum(noise2) * (1 - 1 / np.sum(w))
        if var01 > 0 and var02 > 0: #!!! 万一出现根号里面小于0,就强行不扣噪声
            var1 = var01
            var2 = var02

    coeff = cov / np.sqrt(var1 * var2)

    return coeff

def weighed_corr_map(w1, w2, l1, l2, lag, noise1, noise2, noil, noir, flagnoi): #!!! noil和noir是扣噪声区域的左右边界,单位是bin,并且以脉冲起始位置left作为0
    ll1 = l1.T
    ll2 = l2.T
    #cr,n分别记录全部相关图/负相关部分的相关图
    cr = np.zeros((len(ll1), len(ll2)))
    n = np.zeros((len(ll1), len(ll2)))
    for j in range(len(ll1)):
        for k in range(len(ll2)):
            if flagnoi == 0 and j == k:
                flagnoi = 0
            elif j < noil or j > noir or k < noil or k > noir:
                flagnoi = 0
            else:
                flagnoi = 1
            coe = weighed_corr(w1, w2, ll1[j], ll2[k], lag, noise1, noise2, flagnoi)
            cr[j][k] = coe
            if coe < 0:
                n[j][k] = coe
    return cr, n

I1=np.load('E:/pulsar/correlation/%s/npy/%s_I.npy'%(star,chn_1))
I2=np.load('E:/pulsar/correlation/%s/npy/%s_I.npy'%(star,chn_2))
#计算每个周期噪声的方差
noise1 = np.zeros(I1.shape[0])
noise2 = np.zeros(I2.shape[0])
for i in range(I1.shape[0]):
    noise1[i] = (np.sum((I1[i, 0:left]) * (I1[i, 0:left])) + np.sum((I1[i, right:]) * (I1[i, right:]))) / (left + 1024 - right)
for i in range(I2.shape[0]):
    noise2[i] = (np.sum((I2[i, 0:left]) * (I2[i, 0:left])) + np.sum((I2[i, right:]) * (I2[i, right:]))) / (
                left + 1024 - right)

peak1=I1[:,left:right]
peak2=I2[:,left:right]
#计算每个周期内，脉冲部分的平均值
peak_ave1=peak1.mean(1)
peak_ave2=peak2.mean(1)
#设置初始权重为1，再将有问题的周期权重设为0
w1=np.ones(I1.shape[0])
w2=np.ones(I2.shape[0])
thre=0.01
for j in range(I1.shape[0]):
    if peak_ave1[j]<thre:
        w1[j]=0
for j in range(I2.shape[0]):
    if peak_ave2[j]<thre:
        w2[j]=0
phase=np.linspace(0,360,1024)
phase_peak=phase[left:right]

#相关图的峰对应横纵坐标大致所在的bin（减去脉冲起始的值），色图使用的色标，相关性正负(p正n负)，去噪声区
pha_peak1=683-left
pha_peak2=725-left
cmap='seismic'
pn='n'
noil = 648 - left #!!!设置扣噪声的区域
noir = 759 - left
#J0738-4042
'''pha_peak1=654-left
pha_peak2=697-left
cmap='seismic'
pn='p'
''' #J1825-0935
#检索相关图上峰的bin的范围
pha_range=15
#用于生成去噪声区的分界线
corll = [(noil + left) * 360 / 1024, (noil + left) * 360 / 1024]
corlr = [(noil + left) * 360 / 1024, (noir + left) * 360 / 1024]
corrr = [(noir + left) * 360 / 1024, (noir + left) * 360 / 1024]
#计算相关图
workbook=xlwt.Workbook()
norm=mpl.colors.Normalize(vmin=1,vmax=-1)
for lag in lags:
    print(lag)
    if chn_1 == chn_2 and lag == 0:
        flagnoi = 0
    else:
        flagnoi = 1
    cr, n=weighed_corr_map(w1, w2, peak1, peak2, lag, noise1, noise2, noil, noir, flagnoi)
    plt.imshow(cr, norm=norm, cmap=cmap, aspect='auto', origin='lower',
               extent=(phase_peak[0], phase_peak[-1], phase_peak[0], phase_peak[-1]))
    plt.plot(corll, corlr, c='g', linestyle='--', linewidth=0.5 ) #!!!扣噪声区的界限
    plt.plot(corlr, corrr, c='g', linestyle='--', linewidth=0.5 )
    plt.plot(corrr, corlr, c='g', linestyle='--', linewidth=0.5 )
    plt.plot(corlr, corll, c='g', linestyle='--', linewidth=0.5 )
    plt.colorbar()
    plt.title('%s-%s cross-cor : lag = %d'%(chn_1,chn_2,lag))
    plt.xlabel('phase %s (°)'%chn_2)
    plt.ylabel('phase %s (°)'%chn_1)
    plt.savefig('E:/pulsar/correlation/%s/correlation_map/%s-%s/cor_all_lag=%d'%(star,chn_1,chn_2,lag))
    plt.clf()

    '''plt.imshow(n, cmap='hot', aspect='auto', origin='lower',
               extent=(phase_peak[0], phase_peak[-1], phase_peak[0], phase_peak[-1]))
    plt.colorbar()
    plt.title('%s-%s cross-cor_negative : lag = %d' % (chn_1, chn_2, lag))
    plt.xlabel('phase 1 (°)')
    plt.ylabel('phase 2 (°)')
    plt.savefig('E:/pulsar/correlation/%s/correlation_map/%s-%s/cor_n_lag=%d' % (star, chn_1, chn_2, lag))
    plt.clf()'''

    #将相关数据存入表格中
    sheet_cr=workbook.add_sheet('all_lag=%d'%lag)
    sheet_n=workbook.add_sheet('n_lag=%d'%lag)
    sheet_cr.write(0,1,'phase_%s(bins)'%chn_2)
    sheet_cr.write(1,0,'phase_%s(bins)'%chn_1)
    sheet_n.write(0,1,'phase_%s(bins)'%chn_2)
    sheet_n.write(1,0,'phase_%s(bins)'%chn_1)
    for j in range(right-left):
        sheet_cr.write(0,j+2,left+j)
        sheet_cr.write(j+2,0,left+j)
        sheet_n.write(0, j + 2, left + j)
        sheet_n.write(j + 2, 0, left + j)
    for j in range(right-left):
        for k in range(right-left):
            sheet_cr.write(j+2,k+2,cr[j,k])
            sheet_n.write(j+2,k+2,n[j,k])

    #检索相关图上两个峰的位置
    if pn=='n':
        cr=-n
    jx0=pha_peak1
    jy0=pha_peak2
    cor_max=cr[jx0,jy0]
    for jx in range(pha_peak1-pha_range,pha_peak1+pha_range):
        for jy in range(pha_peak2-pha_range,pha_peak2+pha_range):
            if jx>jy:
                continue
            if cr[jx,jy]>cr[jx0,jy0]:
                jx0=jx
                jy0=jy
                cor_max=cr[jx,jy]
    print('lag=',lag,' first phase/bin:',round((jx0+left)*360/1024,1),'/',jx0+left)
    print('second phase/bin:',round((jy0+left)*360/1024,1),'/',jy0+left)
    print('correlation=',cor_max)

    jx0=pha_peak2
    jy0=pha_peak1
    cor_max=cr[jx0,jy0]
    for jx in range(pha_peak2-pha_range,pha_peak2+pha_range):
        for jy in range(pha_peak1-pha_range,pha_peak1+pha_range):
            if jx<jy:
                continue
            if cr[jx,jy]>cr[jx0,jy0]:
                jx0=jx
                jy0=jy
                cor_max=cr[jx,jy]
    print('lag=',-lag,' first phase/bin:',round((jy0+left)*360/1024,1),'/',jy0+left)
    print('second phase/bin:',round((jx0+left)*360/1024,1),'/',jx0+left)
    print('correlation=',cor_max)

workbook.save('E:/pulsar/correlation/%s/correlation_map/%s-%s/cor.xls'%(star,chn_1,chn_2))