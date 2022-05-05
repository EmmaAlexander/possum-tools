#from functions import fit, peakfit, quadfit, image_rms, cube_rms,fitsopen, finitelims,vmax,get_rm_scale,burnkfit,weighted_sample_avg_std
from functions import *
import argparse
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs import utils
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import Ellipse
import warnings
warnings.filterwarnings('ignore')
import os
from os import path
from matplotlib import colors
import cmasher as cmr

rootdirectory='/Volumes/TARDIS/Work/askap/leakage_corrected/'
defaultfile='/Users/emma/GitHub/possum-tools/DataProcess/pilot_sources.txt'
fig_directory='/Users/emma/OneDrive/postdoc/pilotI/'

sig_crop_pix=20
dpi=300

possum_freq=1.367490741
emu_freq=0.94349

#set parameters for plots
plasmacmap=plt.cm.plasma
plasmacmap.set_bad('grey',0.5)
wisecmap=plt.cm.Greys
wisecmap.set_bad('white',1)
rmcmap=plt.cm.seismic
rmcmap.set_bad('gray',0.5)
magmacmap=plt.cm.magma
magmacmap.set_bad('gray',0.5)
inferno=plt.cm.inferno
inferno.set_bad('gray',0.5)
viridis=plt.cm.viridis
viridis.set_bad('gray',0.5)
twilight=plt.cm.twilight
twilight.set_bad('gray',0.5)
cividis=plt.cm.cividis
cividis.set_bad('white',0.5)
cmap_redblue=cmr.get_sub_cmap('twilight_shifted', 0.1, 0.9)
cmap_redblue.set_bad('gray',0.5)
bluecmap=plt.cm.Blues
bluecmap.set_bad('white',0)
redcmap=plt.cm.Reds
redcmap.set_bad('white',0)
yellowcmap=plt.cm.YlOrBr_r
yellowcmap.set_bad('white',0)
purplecmap=plt.cm.Purples
purplecmap.set_bad('white',0)
greyscale=plt.cm.gray
greyscale.set_bad('black',1)
cmap1=cmr.get_sub_cmap('plasma', 0, 0.333)
cmap1.set_bad('white',0)
cmap2=cmr.get_sub_cmap('plasma', 0.333, 0.666)
cmap2.set_bad('white',0)
cmap3=cmr.get_sub_cmap('plasma', 0.666, 1)
cmap3.set_bad('white',0)

plt.rcParams.update({'lines.linewidth':0.3})

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'cmu serif'

SMALL_SIZE = 9
MEDIUM_SIZE = 10
BIGGER_SIZE = 11

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

def main(args):
	###########################################################################
	#read in sources
	print('Reading in source list')
	sources=np.loadtxt(args.filename,dtype='str')
	nsrc=sources.shape[0]
	###########################################################################
	for i in range(1,nsrc):
		try:
			src=sources[i,0]
			print('Working on source {} ({}/{})'.format(src,i+1,nsrc))
			directory=rootdirectory+src+'/'
			#Read in values form source list
			POSSUMSB=sources[i,3]
			EMUSB=sources[i,4]
			LAS=float(sources[i,5]) #arcmin
			FOV=float(sources[i,7]) #degrees
			z=float(sources[i,8]) 
			EMUrms=float(sources[i,9]) 
			possum_rms=float(sources[i,10]) 
			rm_scale=float(sources[i,11]) 

			xmin=float(sources[i,12]) 
			ymin=float(sources[i,13]) 
			xmax=float(sources[i,14]) 
			ymax=float(sources[i,15]) 

			if np.isfinite(rm_scale)==False:
				rm_scale_calc=True
			else:
				rm_scale_calc=False
			
			#Read in data
			print("Reading in the fits files")

			#print(directory+src+'_POSSUM.fits')
			#continuum image 

			possum_im,possum_header=fitsopen(directory+src+'_POSSUM_smooth13.fits')
			ithresh=3*possum_rms
			wcs=WCS(possum_header)

			emu_im,emu_head=fitsopen(directory+src+'_EMU_smooth13.fits')

			#other plot params
			#contours
			contourmults = [1,2,4,8,16,32,64,128,256]
			#contourmults = [1,np.sqrt(2),2,2*np.sqrt(2),4,4*np.sqrt(2),8,8*np.sqrt(2),16,16*np.sqrt(2),32,32*np.sqrt(2),64,64*np.sqrt(2),128,128*np.sqrt(2),256]
			i_contours = [ithresh * i for i in contourmults]
			neg_contours=[-1*ithresh,np.sqrt(2)*ithresh,2*ithresh]

			specind=np.log10(np.divide(possum_im,emu_im))/np.log10(np.divide(possum_freq,emu_freq))
			specind=np.ma.masked_where(possum_im<ithresh,specind)



			plt.figure(dpi= dpi)
			ax=plt.subplot(projection=wcs)
			c=ax.imshow(specind, origin='lower', cmap=viridis,vmin=-1.5,vmax=0)#norm=colors.LogNorm(vmin=0.001,vmax=0.01))#,vmin=0,vmax=vmax(average_p_map_masked,99.99))
			#ax.contour(cont_im,levels=i_contours,colors='black')
			ax.contour(possum_im,levels=i_contours,colors='black')
			plt.xlabel('Right Ascension [J2000]')
			plt.ylabel('Declination [J2000]')
			plt.gca().set_aspect("equal")
			#smoothed beam
			#ellipse= Ellipse(xy=(posx,posy),width=b_smooth,height=b_smooth,angle=0,edgecolor='black',fc='none',lw=0.5)
			#ax.add_patch(ellipse)
			cbar=plt.colorbar(c,fraction=0.046, pad=0.04)
			cbar.ax.set_ylabel('Spectral index')
			plt.xlim(xmin,xmax)
			plt.ylim(ymin,ymax)
			#plt.show()
			plt.savefig(fig_directory+src+'_specind.png',dpi=dpi,bbox_inches='tight',transparent=True)

		except:
			print("uh oh")
		
###############################################################################
if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("-f","--filename",help='Name of file with source list',
		default=defaultfile)
	args = ap.parse_args()
	main(args)