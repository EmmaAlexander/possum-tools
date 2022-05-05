#------------------------imports------------------------
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle
import warnings
warnings.filterwarnings('ignore')
from astropy import units as u
from astropy.cosmology import Planck18
import math
import img_scale
from functions import *
import argparse
from astropy.wcs import utils
import os
from os import path
import cmasher as cmr
import glob
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel
from astropy.coordinates import SkyCoord



defaultfile='/Users/emma/GitHub/possum-tools/DataProcess/pilot_sources.txt'
xraysrcs=['2129-5053','2144-5637','2052-5246']
WISEfolder='/Volumes/NARNIA/pilot_cutouts/WISEfiles/'
dpi=300

#set parameters for plots
# get colour maps and set bad values to be transparent
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
greyscale_r=plt.cm.gray_r
greyscale_r.set_bad('white',1)
cmap1=cmr.get_sub_cmap('plasma', 0, 0.333)
cmap1.set_bad('white',0)
cmap2=cmr.get_sub_cmap('plasma', 0.333, 0.666)
cmap2.set_bad('white',0)
cmap3=cmr.get_sub_cmap('plasma', 0.666, 1)
cmap3.set_bad('white',0)

#perceptively uniform cmap that will be used for contours
twilightcmap = plt.cm.get_cmap('twilight')

# contour levels
#the step parameter is the factor of 2^step each contour goes up by
# so use step=1 for contours which double each time
contourexps=np.arange(start=0,stop=32,step=0.5)
contourmults=np.power(2,contourexps)

contourexps_1=np.arange(start=0,stop=32,step=1)
contourmults_1=np.power(2,contourexps_1)

contourexps_2=np.arange(start=0,stop=32,step=2)
contourmults_2=np.power(2,contourexps_2)

plt.rcParams.update({'lines.linewidth':0.3})

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'cmu serif'

SMALL_SIZE = 10
MEDIUM_SIZE = 12
BIGGER_SIZE = 14

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
beampad=5

def main(args):
	###########################################################################
	#read in sources
	print('Reading in source list')
	sources=np.loadtxt(args.filename,dtype='str')
	nsrc=sources.shape[0]
	###########################################################################
	sourcecount=0
	for i in range(nsrc-2,nsrc-1):
		if 1+1==2:
			src=sources[i,0]
			print(src)
			ra=sources[i,1]
			dec=sources[i,2]
			coords=SkyCoord(ra,dec)
			POSSUMSB=sources[i,3]
			EMUSB=sources[i,4]
			LAS=float(sources[i,5]) #arcmin
			FOV=float(sources[i,7]) #degrees
			z=float(sources[i,8]) 
			EMUrms=float(sources[i,9]) 
			POSSUMrms=float(sources[i,10]) 
			#rmscale=float(sources[i,11]) 
			xmin=float(sources[i,12])
			ymin=float(sources[i,13])
			xmax=float(sources[i,14])
			ymax=float(sources[i,15])

			ratio=(ymax-ymin)/(xmax-xmin)

			scalebarkpc=float(sources[i,16])

			if z!=None:
				dist = Planck18.comoving_distance(z) #distance in MPc
				scale=Planck18.arcsec_per_kpc_comoving(z) # scale in arcsec per kpc
				sb_length_arcsec=scale*scalebarkpc*u.kpc


			imagefolder='/Volumes/NARNIA/pilot_cutouts/'+src+'/'

			if src in xraysrcs:
				purple_im_folder=imagefolder+'/xmm/pps/'
				xrayims=glob.glob(purple_im_folder+'*M1S*IMAGE_8000.FTZ')
				purple_im=xrayims[0]
			else:
				purple_im=None

			#blue
			try:
				ra_str=str(coords.ra.deg)[0:5]
				dec_str=str(coords.dec.deg)[0:5]
				#WISEim=glob.glob(WISEfolder+'*-w1-*ra{}*dec{}*'.format(ra_str,dec_str))
				#blue_im=WISEim[0]
				blue_im=imagefolder+'2200-5611_wise.fits'
				blue_thresh=5 #threshold for lowest contour (typically 3x rms noise level)
				blue_label=r'WISE 3.4 $\mu$m'
			except:
				try:
					WISEim=glob.glob(imagefolder+'*WISE_3.4.fits'.format(ra_str,dec_str))
					blue_im=WISEim[0]
					blue_thresh=5
					blue_label=r'WISE 3.4 $\mu$m'
				except:
					print("uh oh")


			greyscale_im=imagefolder+src+'_DSS.fits' #fits data file to load and display in greyscale
			greyscale_label='DSS'

			# get DES data 
			
			R,Rhead=fitsopen(glob.glob(imagefolder+'DES*_i*.fits')[0])
			G,Ghead=fitsopen(glob.glob(imagefolder+'DES*_r*.fits')[0])
			B,Bhead=fitsopen(glob.glob(imagefolder+'DES*_g*.fits')[0])
			DES_wcs=WCS(Rhead)

			#bins=np.arange(0,200,1)
			#plt.hist(np.ndarray.flatten(R),alpha=0.3,color='r',bins=bins)
			#plt.hist(np.ndarray.flatten(G),alpha=0.3,color='g',bins=bins)
			#plt.hist(np.ndarray.flatten(B),alpha=0.3,color='b',bins=bins)
			#plt.show()
			#break

			#R=R-np.nanmean(R)
			#G=G-np.nanmean(G)
			#B=B-np.nanmean(B)
			#R=R/2.

			img = np.zeros((R.shape[0], R.shape[1], 3), dtype=float)
			img[:,:,0] = img_scale.linear(R, scale_min=0, scale_max=1000)
			img[:,:,1] = img_scale.linear(G, scale_min=0, scale_max=1000)
			img[:,:,2] = img_scale.linear(B, scale_min=0, scale_max=1000)

			#R2,Rhead2=fitsopen(glob.glob(imagefolder+'DES*_i*.fits')[1])
			#G2,Ghead2=fitsopen(glob.glob(imagefolder+'DES*_r*.fits')[1])
			#B2,Bhead2=fitsopen(glob.glob(imagefolder+'DES*_g*.fits')[1])
			#DES_wcs2=WCS(Rhead2)

			#img2 = np.zeros((R2.shape[0], R2.shape[1], 3), dtype=float)
			#img2[:,:,0] = img_scale.linear(R2, scale_min=0.1, scale_max=500)
			#img2[:,:,1] = img_scale.linear(G2, scale_min=0.1, scale_max=500)
			#img2[:,:,2] = img_scale.linear(B2, scale_min=0.1, scale_max=500)




			plt.figure(dpi= dpi,figsize=(8.25,8.25*ratio))

			# load the data for the background greyscale 
			greyscale_data,greyscale_header=fitsopen(greyscale_im)
			wcs=WCS(greyscale_header)
			pix_scale=greyscale_header['CDELT2']*3600*u.arcsec
			# plot the greyscale
			ax=plt.subplot(projection=wcs)
			ax.imshow(greyscale_data,origin='lower',cmap=greyscale,vmin=vmax(greyscale_data,5),vmax=vmax(greyscale_data,99.9))
			ax.imshow(img,transform=ax.get_transform(DES_wcs),origin='lower')
			#ax.imshow(img2,transform=ax.get_transform(DES_wcs2),origin='lower')

			#purple data
			if purple_im!=None:
				purple_data,purple_header=fitsopen(purple_im)
				purple_wcs=WCS(purple_header)
				kernel = Gaussian2DKernel(x_stddev=1)
				smoothed=convolve(purple_data,kernel)
				purple_contours = [1,2,4,8,16,32,64]
				purple_thresh=1
				ax.contour(smoothed, transform=ax.get_transform(purple_wcs), colors=[twilightcmap(0.4)],levels=purple_contours,linewidths=1)

				# calculate transparency array for map (alphas) by normalising the data to be between 0 and 1
				# any data below the threshold value will have alpha=0
				# if you want to saturate the scale at a value lower than the maximum, change np.nanmax(purple_data) in the line below
				#purple_alphas=np.divide(purple_data-purple_thresh,np.nanmax(purple_data)-purple_thresh)
				#purple_alphas=np.where(purple_alphas<0,0,purple_alphas)
				#ax2=plt.imshow(purple_data,origin='lower',transform=ax.get_transform(purple_wcs),cmap=purplecmap,alpha=purple_alphas)    
				purple_label='XMM-Newton'
				if purple_label!=None:
					ax.plot(-100,-100,'-',c=twilightcmap(0.4),label=purple_label,linewidth=3)

			#blue data
			if blue_im!=None:
				blue_data,blue_header=fitsopen(blue_im)
				blue_wcs=WCS(blue_header)
				blue_contours = [blue_thresh * i for i in contourmults_1]
				ax.contour(blue_data, transform=ax.get_transform(blue_wcs), colors=[twilightcmap(0.2)],levels=blue_contours[:1],linewidths=0.1)

				# calculate transparency array for map (alphas) by normalising the data to be between 0 and 1
				blue_alphas=np.divide(blue_data-blue_thresh,np.nanmax(blue_data)-blue_thresh)
				blue_alphas=np.where(blue_alphas<0,0,blue_alphas)
				blue_alphas=np.where(blue_alphas>1,1,blue_alphas)
				ax3=plt.imshow(blue_data,origin='lower',transform=ax.get_transform(blue_wcs),cmap=bluecmap,alpha=blue_alphas,vmin=-10*np.nanmax(blue_data),vmax=10*np.nanmax(blue_data)) 

				if blue_label!=None:
					ax.plot(-100,-100,'-',c=twilightcmap(0.2),label=blue_label,linewidth=3)


			#add EMU data
			if EMUSB!='NaN':
				red_data,red_header=fitsopen(imagefolder+src+'_EMU.fits')
				red_label='ASKAP 944 MHz'
				red_wcs=WCS(red_header)
				red_thresh=3.*EMUrms
				red_contours = [red_thresh * i for i in contourmults_1]
				ax.contour(red_data, transform=ax.get_transform(red_wcs), colors=[twilightcmap(0.7)],levels=red_contours[:1],linewidths=0.1)
				# calculate transparency array for map (alphas) by normalising the data to be between 0 and 1
				red_alphas=np.divide(red_data-red_thresh,vmax(red_data,99.9)-red_thresh)
				red_alphas=np.where(red_alphas<0,0,red_alphas)
				red_alphas=np.where(red_alphas>1,1,red_alphas)
				red_alphas=0.333*np.sqrt(red_alphas)
				ax5=plt.imshow(red_data,origin='lower',transform=ax.get_transform(red_wcs),cmap=redcmap,alpha=red_alphas,vmin=-10*np.nanmax(red_data),vmax=10*np.nanmax(red_data))

				if red_label!=None:
					ax.plot(-100,-100,'-',c=twilightcmap(0.75),label=red_label,linewidth=3)

				r_bmaj=red_header['BMAJ']*u.degree
				r_bmin=red_header['BMIN']*u.degree
				r_bpa=red_header['BPA']
				r_bmaj_pix=r_bmaj/pix_scale
				r_bmin_pix=r_bmin/pix_scale

				lowerleft=utils.skycoord_to_pixel(utils.pixel_to_skycoord(xmin,ymin,red_wcs),wcs)
				upperright=utils.skycoord_to_pixel(utils.pixel_to_skycoord(xmax,ymax,red_wcs),wcs)

				greyxmin=lowerleft[0]
				greyymin=lowerleft[1]
				greyxmax=upperright[0]
				greyymax=upperright[1]	

			#add possum data
			if POSSUMSB!='NaN' and POSSUMSB!='10035':
				yellow_data,yellow_header=fitsopen(imagefolder+src+'_POSSUM.fits')
				yellow_label='ASKAP 1368 MHz'
				yellow_thresh=3.*POSSUMrms
				yellow_wcs=WCS(yellow_header)
				yellow_contours = [yellow_thresh * i for i in contourmults_1]
				ax.contour(yellow_data, transform=ax.get_transform(yellow_wcs), colors=[twilightcmap(0.85)],levels=yellow_contours[:1],linewidths=0.1)

				# calculate transparency array for map (alphas) by normalising the data to be between 0 and 1
				yellow_alphas=np.divide(yellow_data-yellow_thresh,vmax(yellow_data,99.9)-yellow_thresh)
				yellow_alphas=np.where(yellow_alphas<0,0,yellow_alphas)
				yellow_alphas=np.where(yellow_alphas>1,1,yellow_alphas)
				yellow_alphas=0.333*np.sqrt(yellow_alphas)
				ax5=plt.imshow(yellow_data,origin='lower',transform=ax.get_transform(yellow_wcs),cmap=yellowcmap,alpha=yellow_alphas,vmin=-10*np.nanmax(yellow_data),vmax=10*np.nanmax(yellow_data))

				if yellow_label!=None:
					ax.plot(-100,-100,'-',c='gold',label=yellow_label,linewidth=3)

				y_bmaj=yellow_header['BMAJ']*u.degree
				y_bmin=yellow_header['BMIN']*u.degree
				y_bpa=yellow_header['BPA']
				y_bmaj_pix=y_bmaj/pix_scale
				y_bmin_pix=y_bmin/pix_scale	

				lowerleft=utils.skycoord_to_pixel(utils.pixel_to_skycoord(xmin,ymin,yellow_wcs),wcs)
				upperright=utils.skycoord_to_pixel(utils.pixel_to_skycoord(xmax,ymax,yellow_wcs),wcs)

				greyxmin=lowerleft[0]
				greyymin=lowerleft[1]
				greyxmax=upperright[0]
				greyymax=upperright[1]


			if z!=None: 
				sb_length_pix= sb_length_arcsec/pix_scale
				scalebar=Rectangle(xy=(greyxmax-sb_length_pix-10,greyymin+5),width=sb_length_pix,height=sb_length_pix/40.,edgecolor='none',fc='white',alpha=1)
				ax.add_patch(scalebar)
				try:
					scaletext='{} kpc'.format(int(scalebarkpc))
					plt.annotate(xy=(greyxmax-(sb_length_pix/2.)-10,greyymin+5+(sb_length_pix/12.)),text=scaletext,c='white',ha="center")  
				except:
					print("shrug")
			plt.gca().set_aspect("equal")
			plt.xlabel('Right Ascension [J2000]')
			plt.ylabel('Declination [J2000]')
			ax.set_xlim(greyxmin,greyxmax)
			ax.set_ylim(greyymin,greyymax)
			#plt.legend(loc=3)
			plt.savefig('/Volumes/NARNIA/pilot_cutouts/'+src+'_DES.png',dpi=dpi,transparent=True,bbox_inches='tight')
			plt.savefig('/Users/emma/OneDrive/PhD/thesis/Figures/new/'+src+'_mwl_RBG.png',dpi=dpi,transparent=True,bbox_inches='tight')

			##############################################################################################
			
			plt.figure(dpi= dpi,figsize=(8.25,8.25*ratio))
			# plot the greyscale
			ax=plt.subplot(projection=wcs)
			ax.imshow(greyscale_data,origin='lower',cmap=greyscale_r,vmin=vmax(greyscale_data,1),vmax=vmax(greyscale_data,99.9))
			ax.imshow(R,transform=ax.get_transform(DES_wcs),origin='lower',vmin=vmax(R,1),vmax=vmax(R,99.5),cmap=greyscale_r)

			#purple data
			if purple_im!=None:

				#ax.contour(smoothed, transform=ax.get_transform(purple_wcs), colors=[twilightcmap(0.4)],levels=purple_contours,linewidths=0.2)
				kernel = Gaussian2DKernel(x_stddev=2)
				smoothed=convolve(purple_data,kernel)
				# calculate transparency array for map (alphas) by normalising the data to be between 0 and 1
				# any data below the threshold value will have alpha=0
				# if you want to saturate the scale at a value lower than the maximum, change np.nanmax(purple_data) in the line below
				purple_alphas=np.divide(smoothed-purple_thresh,vmax(smoothed,99.9)-purple_thresh)
				purple_alphas=np.where(purple_alphas<0,0,purple_alphas)
				purple_alphas=np.where(purple_alphas>1,1,purple_alphas)
				purple_alphas=0.5*np.sqrt(purple_alphas)
				ax2=plt.imshow(purple_data,origin='lower',transform=ax.get_transform(purple_wcs),cmap=cmap1,alpha=purple_alphas,vmin=-10*np.nanmax(purple_data),vmax=10*np.nanmax(purple_data))    
				purple_label='XMM-Newton'
				ax.contour(smoothed, transform=ax.get_transform(purple_wcs), colors=[twilightcmap(0.4)],levels=purple_contours,linewidths=0.8)

				if purple_label!=None:
					ax.plot(-100,-100,'-',c=twilightcmap(0.4),label=purple_label,linewidth=3)


			#blue data
			if blue_im!=None:
				blue_data,blue_header=fitsopen(blue_im)
				blue_wcs=WCS(blue_header)
				blue_contours = [blue_thresh * i for i in contourmults]
				ax.contour(blue_data, transform=ax.get_transform(blue_wcs), colors=[twilightcmap(0.2)],levels=blue_contours,linewidths=0.8)

				# calculate transparency array for map (alphas) by normalising the data to be between 0 and 1
				blue_alphas=np.divide(blue_data-blue_thresh,np.nanmax(blue_data)-blue_thresh)
				blue_alphas=np.sqrt(np.where(blue_alphas<0,0,blue_alphas))
				ax3=plt.imshow(blue_data,origin='lower',transform=ax.get_transform(blue_wcs),cmap=bluecmap,alpha=blue_alphas,vmin=-10*np.nanmax(blue_data),vmax=10*np.nanmax(blue_data)) 

				if blue_label!=None:
					ax.plot(-100,-100,'-',c=twilightcmap(0.2),label=blue_label,linewidth=3)

			#add EMU data
			if EMUSB!='NaN':
				# calculate transparency array for map (alphas) by normalising the data to be between 0 and 1
				red_alphas=np.divide(red_data-red_thresh,vmax(red_data,99.9)-red_thresh)
				red_alphas=np.where(red_alphas<0,0,red_alphas)
				red_alphas=np.where(red_alphas>1,1,red_alphas)
				red_alphas=0.8*np.sqrt(red_alphas)
				red_contours = [red_thresh * i for i in contourmults]
				ax.contour(red_data, transform=ax.get_transform(red_wcs), colors=[twilightcmap(0.7)],levels=red_contours,linewidths=0.8)

				ax4=plt.imshow(red_data,origin='lower',transform=ax.get_transform(red_wcs),cmap=cmap2,alpha=red_alphas,vmin=-10*np.nanmax(red_data),vmax=10*np.nanmax(red_data))

				if red_label!=None:
					ax.plot(-100,-100,'-',c=twilightcmap(0.75),label=red_label,linewidth=3)
				
				ellipse= Ellipse(xy=(greyxmin+5,greyymin+5),width=r_bmaj_pix,height=r_bmin_pix,angle=r_bpa+90,edgecolor=twilightcmap(0.7),fc=twilightcmap(0.7),lw=1)
				ax.add_patch(ellipse)		

			#add possum data
			if POSSUMSB!='NaN' and POSSUMSB!='10035':
				#ax.contour(yellow_data, transform=ax.get_transform(yellow_wcs), colors=['gold'],levels=yellow_contours,linewidths=0.2)

				# calculate transparency array for map (alphas) by normalising the data to be between 0 and 1
				yellow_alphas=np.divide(yellow_data-yellow_thresh,vmax(yellow_data,99.9)-yellow_thresh)
				yellow_alphas=np.where(yellow_alphas<0,0,yellow_alphas)
				yellow_alphas=np.where(yellow_alphas>1,1,yellow_alphas)
				yellow_alphas=0.5*np.sqrt(yellow_alphas)
				ax5=plt.imshow(yellow_data,origin='lower',transform=ax.get_transform(yellow_wcs),cmap=cmap3,alpha=yellow_alphas,vmin=-10*np.nanmax(yellow_data),vmax=10*np.nanmax(yellow_data))
				yellow_contours = [yellow_thresh * i for i in contourmults]
				ax.contour(yellow_data, transform=ax.get_transform(yellow_wcs), colors=[twilightcmap(0.85)],levels=yellow_contours,linewidths=0.8)

				ellipse= Ellipse(xy=(greyxmin+5,greyymin+5),width=y_bmaj_pix,height=y_bmin_pix,angle=y_bpa+90,edgecolor=twilightcmap(0.85),fc=twilightcmap(0.85),lw=1)
				ax.add_patch(ellipse)	

			if z!=None: 
				sb_length_pix= sb_length_arcsec/pix_scale
				scalebar=Rectangle(xy=(greyxmax-sb_length_pix-10,greyymin+5),width=sb_length_pix,height=sb_length_pix/40.,edgecolor='none',fc='black',alpha=1)
				ax.add_patch(scalebar)
				try:
					scaletext='{} kpc'.format(int(scalebarkpc))
					plt.annotate(xy=(greyxmax-(sb_length_pix/2.)-10,greyymin+5+(sb_length_pix/12.)),text=scaletext,c='black',ha="center")  
				except:
					print("shrug")

			plt.gca().set_aspect("equal")
			plt.xlabel('Right Ascension [J2000]')
			plt.ylabel('Declination [J2000]')
			ax.set_xlim(greyxmin,greyxmax)
			ax.set_ylim(greyymin,greyymax)
			#plt.legend(loc=2)
			plt.savefig('/Volumes/NARNIA/pilot_cutouts/'+src+'_DES_gs.png',dpi=dpi,transparent=True)
			plt.savefig('/Users/emma/OneDrive/PhD/thesis/Figures/new/'+src+'_DES_gs.png',dpi=dpi,transparent=True)
		
		#except:
			#print("something went wrong")

if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("-f","--filename",help='Name of file with source list',
		default=defaultfile)
	args = ap.parse_args()
	main(args)


