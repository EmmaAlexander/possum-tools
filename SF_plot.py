import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits
import math as mth
import warnings
warnings.filterwarnings('ignore')
from scipy import spatial
from astropy.cosmology import Planck18
from astropy import units as u


srcs=['1705-2445','1709-2226','1712-2435','1713-2502']
zs=[0.036292,0.030955,0.06063,0.0286]
marker=['o','o','o','o']

cmap=plt.cm.twilight

rootdirectory='/Volumes/TARDIS/Work/askap/'
fig_directory='/Users/emma/OneDrive/PhD/thesis/Figures/'


plt.rcParams.update({'lines.linewidth':0.3})
plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'cmu serif'

SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 12
plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize


plt.figure(dpi= 400)
for i in range(0,len(srcs)):
	src=srcs[i]
	z=zs[i]

	SF=np.loadtxt(rootdirectory+src+'/'+src+'SF_witherrs_new.txt')
	pix=SF[0,:]
	arcsec_per_kpc=Planck18.arcsec_per_kpc_comoving(z) # scale in arcsec per kpc
	pix_size=4*u.arcsec #arcsec /pix
	scale=pix_size/arcsec_per_kpc
	kpc=scale*pix/u.kpc
	markersizes=50*SF[3,:]/np.nanmax(SF[3,:])
	color=cmap((0.5+i)/(len(srcs)+0.5))

	#rrors=np.sqrt(np.square(SF[2,:])+np.square(41.6*SF[4,:]))

	#plt.errorbar(kpc,SF[1,:],SF[2,:],color=color,linestyle='none')
	plt.scatter(kpc,SF[1,:],marker='s',s=12,color=color,label=src)
	plt.scatter(kpc,SF[4,:],marker='x',s=12,color='k',label=src)
	#plt.errorbar(kpc,SF[1,:],errors,linestyle='none',color=color)
	#plt.plot(kpc,np.log(SF[3,:]),color=color)
	#medianline=np.nanmedian(SF[3,:])
	#medianstd=np.nanstd(SF[3,:])/len(SF[3,:])

	#lower=np.log10(medianline-medianstd)
	#upper=np.log10(medianline+medianstd)
	#plt.fill_between([0,1500],lower,upper,alpha=0.2,color=color)
	#plt.axhline(y=np.log10(medianline),color=color)
plt.xlabel('Scale [kpc]')
plt.xlim(0,200)
plt.ylim(10,200000)
plt.yscale('log')
#plt.axhline(y=10.8,linestyle='dotted',color='b')
plt.ylabel(r'SF [rad$^{2}$ m$^{-4}$]')
plt.legend(loc='lower right')
plt.savefig(fig_directory+'Oph_SF_newerrs.png',dpi=300,transparent='True')

plt.figure(dpi= 400)
for i in range(0,len(srcs)):
	src=srcs[i]
	z=zs[i]

	SF=np.loadtxt(rootdirectory+src+'/'+src+'SF_witherrs.txt')
	pix=SF[0,:]
	pix_size=4*u.arcsec #arcsec /pix
	scale=pix_size/arcsec_per_kpc
	kpc=scale*pix/u.kpc
	arcmins=4*SF[0,:]/20
	markersizes=50*SF[3,:]/np.nanmax(SF[3,:])
	color=cmap((1+i)/(len(srcs)+1))

	plt.scatter(arcmins,SF[1,:],marker=marker[i],s=markersizes,color=color,label=src)

plt.xlabel('Scale [arcminutes]')
#plt.xlim(0,500)
plt.xlim(0,50)
plt.ylim(20,200000)
plt.yscale('log')
#plt.axhline(y=10.8,linestyle='dotted',color='b')
plt.ylabel(r'SF [rad$^{2}$ m$^{-4}$]')
plt.legend(loc='lower right')
plt.savefig(fig_directory+'Oph_SF_arcmin.png',dpi=300,transparent='True')

