import sys
import numpy as np
import numpy
sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages')
sys.path.append('/home/aklantbhowmick/anaconda3/lib/python3.7/site-packages/scalpy/')
sys.path.append('/home/aklantbhowmick/anaconda3/envs/nbodykit-env/lib/python3.6/site-packages/')
import matplotlib as mpl

import arepo_package
import scipy.interpolate
from matplotlib import pyplot as plt



def savitzky_golay(y, window_size, order, deriv=0, rate=1):

    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def mean_plot(x,y,xscl,yscl,nbins):
    #nbins = 5
    if(yscl==True):
        y=np.log10(y)
    if(xscl==True):
        x=np.log10(x)
    n, _ = np.histogram(x, bins=nbins)
    sy, _ = np.histogram(x, bins=nbins, weights=y)
    sy2, _ = np.histogram(x, bins=nbins, weights=y*y)
    mean = sy / n
    #std = np.sqrt(sy2/n - mean*mean)
    std=1/numpy.sqrt(n)
    #plt.plot(x, y, 'bo')
    #plt.errorbar((_[1:] + _[:-1])/2, mean,std, color='blue', label = 'z = 8')
    #mean= savitzky_golay(mean, 11, 3)
    #print (_[1:] + _[:-1])/2
    #print mean
    mask=(std/mean)*100<100.
    
    x=((_[1:] + _[:-1])/2)[mask]
    y=mean[mask]
    yul=y+std[mask]
    yll=y-std[mask]
    return x,y#,plt.errorbar(x,y,color=colour,linewidth=3,label=labl),plt.fill_between(x,yll, yul,color=colour,alpha=0.5)
    
    #legend(frameon = False, loc = 'upper right', ncol = 1)
   
    #legend(frameon = False, loc = 'upper right', ncol = 2)
    #plt.show()
    #yscale('log')
    #xlabel(xlabl,fontsize=20)
    #ylabel(ylabl,fontsize=20)
    
def mean_plot_t(x,y,colour,labl,xlabl,ylabl,xscl,yscl,nbins):
    #nbins = 5
    if(yscl==True):
        y=log10(y)
    if(xscl==True):
        x=log10(x)
    n, _ = np.histogram(x, bins=nbins)
    sy, _ = np.histogram(x, bins=nbins, weights=y)
    sy2, _ = np.histogram(x, bins=nbins, weights=y*y)
    mean = sy / n
    #std = np.sqrt(sy2/n - mean*mean)
    std=1/numpy.sqrt(n)
    #plt.plot(x, y, 'bo')
    #plt.errorbar((_[1:] + _[:-1])/2, mean,std, color='blue', label = 'z = 8')
    #mean= savitzky_golay(mean, 11, 3)
    #print (_[1:] + _[:-1])/2
    #print mean
    mask=(std/mean)*100<2
    
    x=((_[1:] + _[:-1])/2)[mask]
    y=mean[mask]
    yul=y+std[mask]
    yll=y-std[mask]
    return y,x#plt.errorbar(y,x,color=colour,linewidth=3,label=labl)#,plt.fill_between(x,yll, yul,color=colour,alpha=0.5)
    
    #legend(frameon = False, loc = 'upper right', ncol = 1)
   
    #legend(frameon = False, loc = 'upper right', ncol = 2)
    #plt.show()
    #yscale('log')
    #xlabel(xlabl,fontsize=20)
    #ylabel(ylabl,fontsize=20)


basePath='/ufrc/lblecha/aklantbhowmick/arepo_runs_aklant/L25_n256/output/'

basePath='/n/ghernquist/Illustris/Runs/Illustris-1/'

run='L205n2500TNG'

basePath='/n/hernquistfs3/IllustrisTNG/Runs/'+run+'/output/'

subhalo_property='SubhaloBHMass'
desired_redshift=0.
SubhaloBHMass,output_redshift=arepo_package.get_subhalo_property(basePath, subhalo_property, desired_redshift, list_all=True)
subhalo_property='SubhaloMass'
SubhaloMass,output_redshift=arepo_package.get_subhalo_property(basePath, subhalo_property, desired_redshift, list_all=True)


subhalo_property='SubhaloMassType'
SubhaloMassType,output_redshift=arepo_package.get_subhalo_property(basePath, subhalo_property, desired_redshift, list_all=True)
SubhaloStellarMass=SubhaloMassType[:,4]

Subhalo5Mass=SubhaloMassType[:,5]


f,ax=plt.subplots(figsize=(12,10))

opacity=1
NBINS=20

mask=(SubhaloMass>1e8/1e10) & (SubhaloBHMass>0)
colormap='Reds_r'
opacity=1
panel2=ax.hist2d(numpy.log10(SubhaloMass[mask]*1e10),numpy.log10(SubhaloBHMass[mask]*1e10), bins=(NBINS,NBINS), norm=mpl.colors.LogNorm(),cmap=colormap,alpha=opacity)

mean_SubhaloMass,mean_SubhaloBHMass=mean_plot(SubhaloMass[mask]*1e10,SubhaloBHMass[mask]*1e10,True,True,13)
ax.plot(mean_SubhaloMass,mean_SubhaloBHMass,linewidth=2)

ax.tick_params(labelsize=30)

ax.set_xlabel('$\log_{10}M_h[M_{\odot}/h]$',fontsize=30)
ax.set_ylabel('$\log_{10}M_{bh}[M_{\odot}/h]$',fontsize=30)


numpy.save(run+'mean_BHM_HM.npy',[mean_SubhaloMass,mean_SubhaloBHMass])


plt.savefig(run+'BHM_HM.pdf',bbox_inches='tight')
f,ax=plt.subplots(figsize=(12,10))
colormap='Blues_r'
opacity=1
NBINS=20

mask=(SubhaloStellarMass>1e8/1e10) & (SubhaloBHMass>0)
panel2=ax.hist2d(numpy.log10(SubhaloMass[mask]*1e10),numpy.log10(SubhaloStellarMass[mask]*1e10), bins=(NBINS,NBINS), norm=mpl.colors.LogNorm(),cmap=colormap,alpha=opacity)

mean_SubhaloMass,mean_SubhaloBHMass=mean_plot(SubhaloMass[mask]*1e10,SubhaloStellarMass[mask]*1e10,True,True,13)
ax.plot(mean_SubhaloMass,mean_SubhaloBHMass,linewidth=2)

numpy.save(run+'mean_SM_HM.npy',[mean_SubhaloMass,mean_SubhaloBHMass])

ax.tick_params(labelsize=30)
ax.set_xlabel('$\log_{10}M_h[M_{\odot}/h]$',fontsize=30)
ax.set_ylabel('$\log_{10}M_{*}[M_{\odot}/h]$',fontsize=30)

plt.savefig(run+'SM_HM.pdf',bbox_inches='tight')
