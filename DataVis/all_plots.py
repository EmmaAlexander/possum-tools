from functions import *
import argparse
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs import utils
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle
import warnings
warnings.filterwarnings('ignore')
import os
from os import path
from matplotlib import colors
import cmasher as cmr
from astropy import units as u
from astropy.cosmology import Planck18
import math
import img_scale
import argparse
from astropy.wcs import utils
import glob
from astropy.convolution import convolve
from astropy.convolution import Gaussian2DKernel
from astropy.coordinates import SkyCoord

rootdirectory='/Volumes/NARNIA/pilot_cutouts/'
#rootdirectory='/Volumes/TARDIS/Work/askap/'
defaultfile='/Users/emma/GitHub/possum-tools/DataProcess/pilot_sources.txt'
fig_directory='/Users/emma/OneDrive/PhD/thesis/Figures/new/'

galacicRMfile='/Volumes/TARDIS/Work/askap/Faraday_cutout_pilot.fits'
galacticRMerrs='/Volumes/TARDIS/Work/askap/Faraday_error_pilot.fits'
#galacicRMfile='/Volumes/TARDIS/Work/askap/Faraday_cutout_Oph.fits'
#galacticRMerrs='/Volumes/TARDIS/Work/askap/Faraday_error_Oph.fits'

sig_crop_pix=20
xraysrcs=['2129-5053','2144-5637','2052-5246']
WISEfolder='/Volumes/NARNIA/pilot_cutouts/WISEfiles/'
dpi=300

#import colourmaps
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
twilightcmap=plt.cm.twilight
twilightcmap.set_bad('gray',0.5)
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

#update plot parameters
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

beampad=5

def main(args):
	###########################################################################
	#read in sources
	print('Reading in source list')
	sources=np.loadtxt(args.filename,dtype='str')
	nsrc=sources.shape[0]
	###########################################################################
	for i in range(0,nsrc):
		#try:
		if 1==1:
			src=sources[i,0]
			print('Working on source {} ({}/{})'.format(src,i+1,nsrc))
			directory=rootdirectory+src+'/'
			imagefolder=rootdirectory+src+'/'
			#Read in values form source list
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
			possum_rms=POSSUMrms #because I'm silly 
			rmscale=float(sources[i,11]) 
			rm_scale=rmscale #because I'm silly 
			xmin=float(sources[i,12])
			ymin=float(sources[i,13])
			xmax=float(sources[i,14])
			ymax=float(sources[i,15])
			scalebarkpc=float(sources[i,16])

			blue_im=None

			if np.isfinite(rm_scale)==False:
				rm_scale_calc=True
			else:
				rm_scale_calc=False

			if z!=None:
				dist = Planck18.comoving_distance(z) #distance in MPc
				scale=Planck18.arcsec_per_kpc_comoving(z) # scale in arcsec per kpc
				sb_length_arcsec=scale*scalebarkpc*u.kpc
			
			#Read in data
			print("Reading in the fits files")

			#print(directory+src+'_POSSUM.fits')
			#continuum image 	
			try:
				cont_im,header=fitsopen(directory+src+'_POSSUM.fits')
			except:
				cont_im,header=fitsopen(directory+src+'_widePOSSUM.fits')
			ithresh=3*possum_rms
			ithresh=3*EMUrms
			#load the cubes 
			try:
				i_cube_data,cube_head=fitsopen(directory+src+'_POSSUM_i.fits')
				q_cube_data,qcube_head=fitsopen(directory+src+'_POSSUM_q.fits')
				u_cube_data,ucube_head=fitsopen(directory+src+'_POSSUM_u.fits')
			except:
				i_cube_data,cube_head=fitsopen(directory+src+'_i.fits')
				q_cube_data,qcube_head=fitsopen(directory+src+'_q.fits')
				u_cube_data,ucube_head=fitsopen(directory+src+'_u.fits')

			cont_wcs=WCS(header)
		    
			#peak polarised intensity from RMsynth
			pkPI_im,pkPI_header=fitsopen(directory+'FDF_maxPI_dropdeg.fits')

		    #peak rm from RMsynth
			pkrm_im,pkrm_header=fitsopen(directory+'FDF_peakRM_dropdeg.fits')
			wcs=WCS(pkrm_header)
			#fitted peak rm 
			pkrm_fit,pkrm_fit_header=fitsopen(directory+'FDF_peakRM_fitted.fits')

			#fitted peak rm 
			chi0_fit,chi0_fit_header=fitsopen(directory+'FDF_chi0_fitted.fits')

			#fdf
			fdfdata,fdfheader=fitsopen(directory+'FDF_tot_dirty.fits')
			fdf_im,fdfhead=fitsopen(directory+'FDF_im_dirty.fits')
			fdf_real,fdfhead=fitsopen(directory+'FDF_real_dirty.fits')


			fdf_clean,fdf_clean_head=fitsopen(directory+'FDF_clean_tot.fits')
			clean_components,cc_head=fitsopen(directory+'FDF_CC_tot.fits')
			fdf_width,fdfwidthhead=fitsopen(directory+'FDF_peakRM_fitted_width.fits')


			rmsf_FHWMarr,rmsfhead=fitsopen(directory+'RMSF_FWHM.fits')
			rmsf_FWHM=np.nanmean(rmsf_FHWMarr)

			#######################################################################
			#beam
			pixscale=header['CDELT2']
			bmaj_pix=header['BMAJ']/header['CDELT2']
			bmin_pix=header['BMIN']/header['CDELT2']

			b_smooth=pkrm_header['BMAJ']/header['CDELT2']

			beam_pix_area=np.pi*(bmaj_pix/2.)*(bmin_pix/2.)
			bpa=header['BPA']
			pad=1.
			posx=0.5+pad+bmaj_pix/2.
			posy=0.5+pad+bmin_pix/2.

			#frequencies
			ref_freq=cube_head['CRVAL3']
			ref_freq_pos=cube_head['CRPIX3']
			delta_freq=cube_head['CDELT3']
			nfreqchans=cube_head['NAXIS3']
			freqs=np.divide(np.linspace(ref_freq,(ref_freq+(delta_freq*nfreqchans)),
				nfreqchans),1000000000)
			c=299792458 #m/s
			lamda=np.divide(c,(freqs*1000000000))
			lamdasq=np.square(lamda)

			#fdf
			fdf_ref_rm=fdfheader['CRVAL3']
			fdf_ref_chan=fdfheader['CRPIX3']
			fdfres=fdfheader['CDELT3']
			nrmchans=fdfheader['NAXIS3']
			fdfmin=fdf_ref_rm-(fdf_ref_chan*fdfres) + 1
			rmarray=np.arange(fdfmin,fdfmin+(fdfres*nrmchans),fdfres)

			#other plot params
			#contours
			contourexps=np.arange(start=0,stop=32,step=0.5)
			contourmults=np.power(2,contourexps)

			contourexps_1=np.arange(start=0,stop=32,step=1)
			contourmults_1=np.power(2,contourexps_1)

			contourexps_2=np.arange(start=0,stop=32,step=2)
			contourmults_2=np.power(2,contourexps_2)


			i_contours = [ithresh * i for i in contourmults]
			neg_contours=[-1*ithresh,np.sqrt(2)*ithresh,2*ithresh]
			#pol intensity

			#x/y min/max
			#pkrm_index=np.where(np.isfinite(pkrm_fit))

			#######################################################################
			#calculated maps and parameters
			print("Calculating some things")
			# sample the edges of the image to get a per channel noise level
			I_cube_rms=cube_rms(i_cube_data,sig_crop_pix)
			U_cube_rms=cube_rms(u_cube_data,sig_crop_pix)
			Q_cube_rms=cube_rms(q_cube_data,sig_crop_pix)
			print("rms done")
			QUrms=np.nanmean([U_cube_rms,Q_cube_rms],axis=0)
			QU_rms_arr=np.multiply(QUrms[:,np.newaxis,np.newaxis],np.ones(
				shape=u_cube_data.shape))

			I_rms_arr=np.multiply(I_cube_rms[:,np.newaxis,np.newaxis],np.ones(
				shape=i_cube_data.shape))

			p_cube_uncorr=np.sqrt(np.add(np.square(u_cube_data),np.square(
				q_cube_data)))
			# below is the correction of p rms noise
			p_cube=np.sqrt(np.square(p_cube_uncorr)-np.square(QU_rms_arr))
			print("P CUBE DONE done")
			max_p=np.nanmax(p_cube,axis=0)
			polfrac_cube=np.divide(p_cube,i_cube_data)

			temp1=np.square(np.divide(QU_rms_arr,p_cube))
			temp2=np.square(np.divide(I_rms_arr,i_cube_data))
			polfrac_cube_errs=np.multiply(np.sqrt(temp1+temp2),polfrac_cube)
			print("P CUBE ERRS DONE done")

			polfrac_map=np.nanmean(polfrac_cube,axis=0)

			hdu = fits.PrimaryHDU(polfrac_map,header=pkrm_header)
			hdu.writeto(directory+'polfrac.fits',overwrite=True)

			polfrac_max=np.nanmax(polfrac_cube,axis=0)

			#work out how the polarisation fraction changes accross the band
			#polfrac(low frequency)/polfrac(high frequency)
			polfrac_change=np.divide(np.nanmean(polfrac_cube[0:10,:,:],axis=0),np.nanmean(polfrac_cube[-10:,:,:],axis=0))
			#>1 means polarisation fraction is higher at low frequencies
			#

			pa_cube=0.5*np.arctan2(u_cube_data,q_cube_data)
			pa_cube_errs=0.5*(np.divide(QU_rms_arr,p_cube))

			average_p_map=np.nanmean(p_cube,axis=0)
			average_p_map_rms=image_rms(average_p_map,sig_crop_pix)	

			sigmaP_map=np.divide(max_p,average_p_map_rms)

			p_contours = [3. * average_p_map_rms * i for i in contourmults]

			average_i_map=np.nanmean(i_cube_data,axis=0)
			average_i_rms=image_rms(average_i_map,sig_crop_pix)	

			avithresh=3.*average_i_rms
			#if avithresh>ithresh:
			#	avithresh=ithresh

			average_chan_int=np.nanmean(i_cube_data,axis=(1,2))
			bad_chans=np.where(average_chan_int<=0.00000001)
			#mask those channels with nans
			freqs=np.where(average_chan_int<=0.00000001,np.nan,freqs)
			lamda=np.where(average_chan_int<=0.00000001,np.nan,lamda)
			lamdasq=np.where(average_chan_int<=0.00000001,np.nan,lamdasq)


			cont_mask=np.ma.masked_where(average_i_map<avithresh,average_i_map)
			cont_mask=np.ma.masked_where(np.isfinite(pkrm_im)==False,cont_mask)
			cont_pol_mask=np.ma.masked_where(sigmaP_map<3,cont_mask)

			#average_p_map_masked=np.ma.masked_where(polfrac_change<0,average_p_map_masked)
			#average_p_map_masked=np.ma.masked_where(average_p_map_masked<ithresh,average_p_map_masked)

			chi0_fit_masked=np.ma.masked_where(cont_pol_mask.mask==1,chi0_fit)
			pkrm_fit_masked=np.ma.masked_where(cont_pol_mask.mask==1,pkrm_fit)
			fdf_width_masked=np.ma.masked_where(cont_pol_mask.mask==1,fdf_width)

			polfrac_change_masked=np.ma.masked_where(cont_mask.mask==1,polfrac_change)

			polfrac_map_masked=np.ma.masked_where(cont_mask.mask==1,polfrac_map)
			polfrac_map_masked=np.ma.masked_where(polfrac_map_masked>=2,polfrac_map_masked)
			polfrac_map_masked=np.ma.masked_where(polfrac_map_masked<=0,polfrac_map_masked)


			polfrac_max_masked=np.ma.masked_where(cont_pol_mask.mask==1,polfrac_max)
			#polfrac_change_masked=np.ma.masked_where(np.isfinite(pkrm_fit_masked)==False,polfrac_change_masked)
			average_p_map_masked=np.ma.masked_where(cont_pol_mask.mask==1,average_p_map)
			#burn k law fit
			#emptymask=0*np.ones(shape=average_p_map_masked.shape)
			if os.path.isfile(directory+'burn_k.fits')==True and os.path.isfile(directory+'burn_k_chi2.fits')==True and os.path.isfile(directory+'burn_k_p0.fits')==True:
				k_arr,khead=fitsopen(directory+'burn_k.fits')
				burnredchi,bhead=fitsopen(directory+'burn_k_chi2.fits')
				p0_arr,phead=fitsopen(directory+'burn_k_p0.fits')
			else:
				print("fitting burn k law")
				k_arr,p0_arr,burnredchi=burnkfit(lamdasq,polfrac_cube,polfrac_cube_errs,cont_mask.mask)
				hdu = fits.PrimaryHDU(k_arr,header=pkrm_header)
				hdu.writeto(directory+'burn_k.fits',overwrite=True)
				hdu2 = fits.PrimaryHDU(burnredchi,header=pkrm_header)
				hdu2.writeto(directory+'burn_k_chi2.fits',overwrite=True)
				hdu3 = fits.PrimaryHDU(p0_arr,header=pkrm_header)
				hdu3.writeto(directory+'burn_k_p0.fits',overwrite=True)

			'''
			try:

				if os.path.isfile(directory+'inband_alpha.fits'==True):
					alpha_arr,ahead=fitsopen(directory+'inband_alpha.fits')
					alpha_redchi,achead=fitsopen(directory+'inband_alpha_redchi2.fits')
				else:
					print("fitting in-band spec ind")
					alpha_arr,alpha_redchi=specind(freqs,i_cube_data,I_rms_arr,cont_mask.mask)
					print("saving")
					hdu4=fits.PrimaryHDU(alpha_arr,header=pkrm_header)
					hdu4.writeto(directory+'inband_alpha.fits',overwrite=True)
					hdu5=fits.PrimaryHDU(alpha_redchi,header=pkrm_header)
					hdu5.writeto(directory+'inband_alpha_redchi2.fits',overwrite=True)
			except:
				print("debug time!!")
			'''
			
			print("got here")
			avi_5sig=np.ma.masked_where(average_i_map>=5*ithresh,average_i_map)

			#find the limits of where there is finite data
			if np.isfinite(xmin) ==False:
				xmin,xmax,ymin,ymax=finitelims(avi_5sig)
				print(xmin,xmax,ymin,ymax)
			posx=posx+xmin
			posy=posy+ymin

			if np.nanstd(rmsf_FHWMarr) <= 1:
				RMerr_arr=np.divide(rmsf_FWHM,(2*sigmaP_map))
			else:
				print("You need to check the RMSF FWHM")

			#SAVE THE RM ERROR MAP
			print("Saving the RM error map")
			hdu = fits.PrimaryHDU(RMerr_arr,header=pkrm_header)
			hdu.writeto(directory+'RM_error.fits',overwrite=True)

			RMerr_arr_masked=np.ma.masked_where(cont_mask.mask==1,RMerr_arr)

			#CC complexity 

			rm_broadcast=np.multiply(rmarray[:,np.newaxis,np.newaxis],np.ones(
				shape=clean_components.shape))

			clean_locs=np.where(np.abs(clean_components)>0,rm_broadcast,np.nan)
			cc_stdv=np.nanstd(clean_locs,axis=0)
			cc_stdv=np.ma.masked_where(average_p_map_masked==True,cc_stdv)

			#######################################################################
			#find the most polarised pixel and do example plots:
			exampley,examplex=np.unravel_index(np.nanargmax(average_p_map), average_p_map.shape)
			print(exampley,examplex)
			#examplex=228
			#exampley=372
			#plt.imshow(pkrm_fit_masked,cmap=rmcmap,vmin=-40,vmax=40,origin='lower')
			#plt.show()

			#galactic rm of example pixel
			examplecoords=utils.pixel_to_skycoord(xp=examplex,yp=exampley,wcs=cont_wcs)
			galacticRM,gRMerr=galactic_rm(examplecoords,galacicRMfile,galacticRMerrs)

			slice_i=i_cube_data[:,exampley,examplex]
			slice_q=q_cube_data[:,exampley,examplex]
			slice_u=u_cube_data[:,exampley,examplex]
			slice_p=p_cube[:,exampley,examplex]
			slice_pa=pa_cube[:,exampley,examplex]
			slice_pa_errs=pa_cube_errs[:,exampley,examplex]
			slice_polfrac=polfrac_cube[:,exampley,examplex]
			slice_polfrac_err=polfrac_cube_errs[:,exampley,examplex]

			#slice_i_fracerrs=np.divide(slice_i,I_cube_rms)
			#slice_p_fracerrs=np.divide(slice_p,QUrms)

			#quadratureerrs=np.square(slice_i_fracerrs) + np.square(slice_p_fracerrs)
			
			#slice_polfrac_errs=np.multiply(slice_polfrac,np.sqrt(quadratureerrs))

			#print(slice_polfrac_errs.shape)

			unwrapped = 2*np.copy(slice_pa) 
			unwrapped[~np.isnan(unwrapped)] = 0.5*np.unwrap(unwrapped[~np.isnan(unwrapped)],discont=-np.pi/2.)
			
			#unwrapped=np.copy(slice_pa) 
			try:
				p_rad,pError_rad,redchisqrd_rad=fit(lamdasq,unwrapped,slice_pa_errs+0.001)
			except:
				print("uh oh")


			examplerm=pkrm_fit_masked[exampley,examplex]
			examplermerr=RMerr_arr[exampley,examplex]
			#example_chi0=chi0_fit_masked[exampley,examplex]

			example_burn=k_arr[exampley,examplex]
			example_p0=p0_arr[exampley,examplex]

			burnlawfit=example_p0*np.exp(-1*example_burn*np.square(lamdasq))


			print('Example RM={:+.1f}'.format(examplerm))
			print("Making plots")

			# put things together
			
			clean_lines=np.where(np.abs(clean_components[:,exampley,examplex])>0,rmarray,np.nan)
			clean_lines_max=np.where(np.abs(clean_components[:,exampley,examplex])>0,clean_components[:,exampley,examplex],np.nan)
			clean_lines_min=np.where(np.abs(clean_components[:,exampley,examplex])>0,0,np.nan)
			fig=plt.figure(dpi=400,figsize=(6,10))

			examplesquare=10
			ax1=plt.subplot(421,projection=wcs,adjustable='box')
			ax1im=ax1.imshow(average_p_map, origin='lower', cmap=magmacmap,vmin=0,vmax=np.nanmax(average_p_map[exampley-examplesquare:exampley+examplesquare,examplex-examplesquare:examplex+examplesquare]))
			ax1.contour(cont_im,levels=i_contours,colors='white',transform=ax1.get_transform(cont_wcs))
			ax1.scatter(examplex,exampley,marker='o',edgecolor='black',facecolor='none')
			ax1.set_xlim(examplex-examplesquare,examplex+examplesquare)
			ax1.set_ylim(exampley-examplesquare,exampley+examplesquare)
			ax1.set_xlabel('Right Ascension [J2000]')
			ax1.set_ylabel('Declination [J2000]')
			cbar1=plt.colorbar(ax1im,fraction=0.0475, pad=0.005)
			cbar1.ax.set_ylabel('Polarisation [Jy/bm]')
			ax1.set_title('a) Pixel location')	

			rm_cut_a=pkrm_fit_masked[exampley-examplesquare:exampley+examplesquare,examplex-examplesquare:examplex+examplesquare]
			absrm=np.abs(rm_cut_a)
			rm_cut=rm_cut_a[absrm<=200]
			absrm_cut=absrm[absrm<=200]

			if rm_scale_calc == True:
				rm_scale=absrm_cut
			try:
				rm_min=np.nanmin(rm_cut)
				rm_max=np.nanmax(rm_cut)
			except:
				print("uh oh")

			ax2=plt.subplot(423,projection=wcs,adjustable='box')
			ax2im=ax2.imshow(pkrm_fit_masked, origin='lower', cmap=cmap_redblue,vmin=-1*rm_scale,vmax=rm_scale)
			ax2.contour(cont_im,levels=i_contours,colors='black',transform=ax2.get_transform(cont_wcs))
			ax2.scatter(examplex,exampley,marker='o',edgecolor='white',facecolor='none')
			ax2.set_xlim(examplex-examplesquare,examplex+examplesquare)
			ax2.set_ylim(exampley-examplesquare,exampley+examplesquare)
			ax2.set_xlabel('Right Ascension [J2000]')
			ax2.set_ylabel('Declination [J2000]')
			#ax2.set_ylabel(' ')
			cbar2=plt.colorbar(ax2im,fraction=0.0475, pad=0.005)
			cbar2.ax.set_ylabel(r'Peak RM [rad m$^{-2}$]')
			ax2.set_title('c) Pixel location')#,\n RM value = {:+.1f}'.format(examplerm)+r'$\,\pm\,$'+'{:.1f}'.format(examplermerr)+r' rad m$^{-2}$')
			print(r'RM value (synthesis) = {:+.1f}'.format(examplerm)+r'$\,\pm\,$'+'{:.1f}'.format(examplermerr)+r' rad m$^{-2}$')
			ax3=plt.subplot(426,adjustable='box')
			ax3.errorbar(lamdasq,slice_p,QUrms,linestyle="None",color=twilight(0.5),alpha=0.5)
			ax3.errorbar(lamdasq,slice_q,QUrms,linestyle="None",color=twilight(0.2),alpha=0.5)
			ax3.errorbar(lamdasq,slice_u,QUrms,linestyle="None",color=twilight(0.7),alpha=0.5)
			ax3.scatter(lamdasq,slice_p,label='P',marker='.',color=twilight(0.5),s=10)
			ax3.scatter(lamdasq,slice_q,label='Q',marker='+',linewidth=1,color=twilight(0.2),s=10)
			ax3.scatter(lamdasq,slice_u,label='U',marker='x',linewidth=1,color=twilight(0.7),s=10)
			ax3.axhline(y=0,color='k')
			ax3.set_xlim(np.nanmin(lamdasq),np.nanmax(lamdasq))
			ax3.set_xlabel(r'$\lambda^2$ [m$^2$]')
			ax3.set_ylabel('Intensity [Jy/bm]')
			ax3.set_title('f) Pixel Stokes values')
			ax3.legend()


			if np.nanmax(slice_polfrac)>=1:
				ymax_polfrac=1
			else:
				ymax_polfrac=1.05*np.nanmax(slice_polfrac)

			if np.nanmin(slice_polfrac)<=0:
				ymin_polfrac=0
			else:
				ymin_polfrac=0.95*np.nanmin(slice_polfrac)

			ax4=plt.subplot(422,adjustable='box')
			ax4.scatter(lamdasq,slice_polfrac,color=twilight(0.5),s=10)
			ax4.errorbar(lamdasq,slice_polfrac,slice_polfrac_err,color=twilight(0.5),linestyle='none',alpha=0.5)
			ax4.set_xlim(np.nanmin(lamdasq),np.nanmax(lamdasq))
			ax4.set_ylim(ymin_polfrac,ymax_polfrac)
			ax4.set_xlabel(r'$\lambda^2$ [m$^2$]')
			ax4.set_ylabel('Polarisation fraction')
			ax4.set_title('b) Pixel polarisation fraction')

			ax4.plot(lamdasq,burnlawfit,linewidth=1,color=twilight(0.25))

			urange=np.nanmax(slice_u)-np.nanmin(slice_u)
			qrange=np.nanmax(slice_q)-np.nanmin(slice_q)
			qurange=0.5*np.nanmax([urange,qrange])
			umid=0.5*(np.nanmax(slice_u)+np.nanmin(slice_u))
			qmid=0.5*(np.nanmax(slice_q)+np.nanmin(slice_q))

			umin=umid-qurange
			umax=umid+qurange
			qmin=qmid-qurange
			qmax=qmid+qurange

			ax5=plt.subplot(425,adjustable='box')
			ax5im=ax5.scatter(slice_u,slice_q,c=lamdasq,cmap=plt.cm.plasma,s=10)
			ax5.errorbar(slice_u,slice_q,xerr=U_cube_rms,yerr=Q_cube_rms,linestyle='none',color='k')
			ax5.set_xlabel('U [Jy/bm]')
			ax5.set_ylabel('Q [Jy/bm]')
			ax5.set_xlim(umin,umax)
			ax5.set_ylim(qmin,qmax)
			ax5.axhline(y=0,color='k')
			ax5.axvline(x=0,color='k')
			#ax5.set_aspect('equal')
			cbar5=plt.colorbar(ax5im,fraction=0.048, pad=0.005)
			cbar5.ax.set_ylabel(r'$\lambda^2$ [m$^2$]')
			ax5.set_title('e) Pixel Q-U values')

			ax6=plt.subplot(424,adjustable='box')
			#ax6.scatter(lamdasq,slice_pa,color='b')
			#ax6.errorbar(lamdasq,slice_pa,slice_pa_errs,linestyle='None',color='b')
			ax6.scatter(lamdasq,unwrapped,color=twilight(0.5),s=10)
			ax6.errorbar(lamdasq,unwrapped,slice_pa_errs,linestyle='None',color=twilight(0.5),alpha=0.5)
			ax6.plot(lamdasq,(p_rad[0]*lamdasq+p_rad[1]),linewidth=1,color=twilight(0.25))
			#ax6.plot(lamdasq,examplerm*lamdasq+example_chi0,linewidth=1,color=twilightcmap(0.75))
			#ax6.plot(lamdasq,examplerm*lamdasq+example_chi0-np.pi,linewidth=1,color=twilightcmap(0.75))
			#ax6.plot(lamdasq,examplerm*lamdasq+example_chi0-2*np.pi,linewidth=1,color=twilightcmap(0.75))
			ax6.set_xlabel(r'$\lambda^2$ [m$^2$]')
			ax6.set_ylabel('Polarisation angle [radians]')
			ax6.set_xlim(np.nanmin(lamdasq),np.nanmax(lamdasq))
			ax6.set_ylim(np.nanmin(unwrapped)-0.1,np.nanmax(unwrapped)+0.1)
			#plt.annotate(text='{:+.1f}'.format(p_rad[0])+r'$\,\pm\,$'+'{:.1f}'.format(pError_rad[0])+r' rad m$^{-2}$',xy=(0.64,0.45),xycoords='figure fraction',color='k',fontsize=SMALL_SIZE)
			ax6.set_title('d) Pixel polarisation angles')#,\n RM value = {:+.1f}'.format(p_rad[0])+r'$\,\pm\,$'+'{:.1f}'.format(pError_rad[0])+r' rad m$^{-2}$')
			print(r'RM value (angle fit) = {:+.1f}'.format(p_rad[0])+r'$\,\pm\,$'+'{:.1f}'.format(pError_rad[0])+r' rad m$^{-2}$')
			ax7=plt.subplot(414)

			ax7.plot(rmarray,fdf_clean[:,exampley,examplex],label='Cleaned',linewidth=1,color=twilight(0.5))
			ax7.plot(rmarray,fdfdata[:,exampley,examplex],label='Dirty',linestyle='dashed',linewidth=0.7,color=twilight(0.15))
			#ax7.axhline(y=0,color='k')
			ax7.vlines(clean_lines,clean_lines_min,clean_lines_max,color=twilight(0.8),linewidth=1,label='Clean components')
			ax7.set_ylabel(r'|P| [mJy beam$^{-1}$ rmsf$^{-1}$]')
			ax7.set_xlabel(r'$\phi$  [rad m$^{-2}]$')
			ax7.set_xlim(-1000,1000)
			ax7.axvline(x=galacticRM,color=twilight(0.3),linestyle=':',linewidth=1,label='Galactic RM')

			try:
				ymaxfdf=1.05*np.nanmax(fdfdata[:,exampley,examplex])
				ax7.set_ylim(0,ymaxfdf)
				ax7.fill_betweenx([0,ymaxfdf],x1=galacticRM-gRMerr, x2=galacticRM+gRMerr,color=twilight(0.3),alpha=0.5)
			except:
				print("uh oh")
			ax7.legend(loc='upper left')
			ax7.set_title('g) Pixel FDF')

			plt.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.95, wspace=0.6, hspace=0.4)
			#plt.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.95, wspace=0.6, hspace=0.35)
			plt.savefig(fig_directory+src+'_example_pix{}_{}'.format(examplex,exampley)+'.png',dpi=400,bbox_inches='tight',transparent=True)
			plt.savefig(directory+src+'_example_pix{}_{}'.format(examplex,exampley)+'.png',dpi=400,bbox_inches='tight',transparent=True)
			

			###################################################################################
			#Maps 

			#widfth of fdf
			fdf_width_median=np.nanmedian(fdf_width)
			fdf_vmin=np.percentile(fdf_width[np.isfinite(fdf_width)],1)
			fdf_vmax=np.percentile(fdf_width[np.isfinite(fdf_width)],90)

			
			#zoomed in pol plot
			chi_0_fdf_deg=np.rad2deg(chi0_fit_masked)
			x=np.cos(np.deg2rad(chi_0_fdf_deg+90))
			y=np.sin(np.deg2rad(chi_0_fdf_deg+90))

			shortest_axis=np.nanmin([xmax-xmin,ymax-ymin])
			step=int(np.floor(shortest_axis/100.))
			if step <= 2:
				step=2

			#step=2
			x_reduced=x[0::step,0::step]
			y_reduced=y[0::step,0::step]
			x_pos=np.arange(start=0,stop=x.shape[1],step=step)
			y_pos=np.arange(start=0,stop=x.shape[0],step=step)

			pol_angle_mask=average_p_map_masked.mask[0::step,0::step]

			x_reduced=np.where(pol_angle_mask==1,np.nan,x_reduced)
			y_reduced=np.where(pol_angle_mask==1,np.nan,y_reduced)

			absrm=np.abs(pkrm_fit_masked)
			rm_cut=pkrm_fit_masked[absrm<=200]
			absrm_cut=absrm[absrm<=200]

			plt.figure(dpi= 600)
			ax=plt.subplot(projection=wcs)
			c=ax.imshow(average_p_map, origin='lower', cmap=magmacmap,vmin=0,vmax=vmax(average_p_map_masked,99.99))#norm=colors.LogNorm(vmin=0.001,vmax=0.01))#,vmin=0,vmax=vmax(average_p_map_masked,99.99))
			q = ax.quiver(x_pos,y_pos,x_reduced,y_reduced,color='white',scale=50,width=0.002,headlength=0,headwidth=0)
			#ax.contour(cont_im,levels=i_contours,colors='black')
			ax.contour(cont_im,levels=i_contours,colors='gray',transform=ax.get_transform(cont_wcs))
			plt.xlabel('Right Ascension [J2000]')
			plt.ylabel('Declination [J2000]')
			plt.gca().set_aspect("equal")
			ellipse= Ellipse(xy=(posx,posy),width=bmaj_pix,height=bmin_pix,angle=bpa+90,edgecolor='gray',fc='gray',lw=0.5)
			ax.add_patch(ellipse)
			#smoothed beam
			ellipse= Ellipse(xy=(posx,posy),width=b_smooth,height=b_smooth,angle=0,edgecolor='black',fc='none',lw=0.5)
			ax.add_patch(ellipse)
			cbar=plt.colorbar(c,fraction=0.046, pad=0.04)
			cbar.ax.set_ylabel('Band-averaged polarisation [Jy/bm]')
			plt.xlim(xmin,xmax)
			plt.ylim(ymin,ymax)
			#plt.show()
			plt.savefig(fig_directory+src+'_avpol_quiver.png',dpi=400,bbox_inches='tight',transparent=True)
			plt.savefig(directory+src+'_avpol_quiver.png',dpi=400,bbox_inches='tight',transparent=True)

			plt.figure(dpi= 600)
			ax=plt.subplot(projection=wcs)
			c=ax.imshow(average_p_map, origin='lower', cmap=magmacmap,vmin=0,vmax=vmax(average_p_map_masked,99.99))
			q = ax.quiver(x_pos,y_pos,x_reduced,y_reduced,color='gray',scale=50,width=0.002,headlength=0,headwidth=0)
			ax.contour(cont_im,levels=i_contours,colors='white',transform=ax.get_transform(cont_wcs))
			plt.xlabel('Right Ascension [J2000]')
			plt.ylabel('Declination [J2000]')
			plt.gca().set_aspect("equal")
			ellipse= Ellipse(xy=(posx,posy),width=bmaj_pix,height=bmin_pix,angle=bpa+90,edgecolor='gray',fc='gray',lw=0.5)
			ax.add_patch(ellipse)
			#smoothed beam
			ellipse= Ellipse(xy=(posx,posy),width=b_smooth,height=b_smooth,angle=0,edgecolor='black',fc='none',lw=0.5)
			ax.add_patch(ellipse)
			cbar=plt.colorbar(c,fraction=0.046, pad=0.04)
			cbar.ax.set_ylabel('Band-averaged polarisation [Jy/bm]')
			#plt.show()
			plt.savefig(fig_directory+src+'_avpol_quiver_full.png',dpi=400,bbox_inches='tight',transparent=True)
			plt.savefig(directory+src+'_avpol_quiver_full.png',dpi=400,bbox_inches='tight',transparent=True)

			#polfrac
			plt.figure(dpi= 600)
			ax=plt.subplot(projection=wcs)
			c=ax.imshow(polfrac_map_masked, origin='lower', cmap=magmacmap,norm=colors.LogNorm(vmin=0.01, vmax=1))#,vmin=0,vmax=1)
			q = ax.quiver(x_pos,y_pos,x_reduced,y_reduced,color='white',scale=50,width=0.002,headlength=0,headwidth=0)
			ax.contour(cont_im,levels=i_contours,colors='gray',transform=ax.get_transform(cont_wcs))
			plt.xlabel('Right Ascension [J2000]')
			plt.ylabel('Declination [J2000]')
			plt.gca().set_aspect("equal")
			ellipse= Ellipse(xy=(posx,posy),width=bmaj_pix,height=bmin_pix,angle=bpa+90,edgecolor='gray',fc='gray',lw=0.5)
			ax.add_patch(ellipse)
			#smoothed beam
			ellipse= Ellipse(xy=(posx,posy),width=b_smooth,height=b_smooth,angle=0,edgecolor='black',fc='none',lw=0.5)
			ax.add_patch(ellipse)
			cbar=plt.colorbar(c,fraction=0.046, pad=0.04)
			cbar.ax.set_ylabel('Average polarisation fraction')
			plt.xlim(xmin,xmax)
			plt.ylim(ymin,ymax)
			#plt.show()
			plt.savefig(fig_directory+src+'_avpolfrac_quiver.png',dpi=400,bbox_inches='tight',transparent=True)
			plt.savefig(directory+src+'_avpolfrac_quiver.png',dpi=400,bbox_inches='tight',transparent=True)

			#chi 0 plot
			chi_0_fdf_deg=np.where(pkrm_fit_masked.mask==True,np.nan,chi_0_fdf_deg)
			plt.figure(dpi= 600)
			ax=plt.subplot(projection=wcs)
			c=ax.imshow(chi_0_fdf_deg, origin='lower', cmap=twilight,vmin=-90,vmax=90)
			ax.contour(cont_im,levels=i_contours,colors='black',transform=ax.get_transform(cont_wcs))
			plt.xlabel('Right Ascension [J2000]')
			plt.ylabel('Declination [J2000]')
			ellipse= Ellipse(xy=(posx,posy),width=bmaj_pix,height=bmin_pix,angle=bpa+90,edgecolor='black',fc='black',lw=0.5)
			ax.add_patch(ellipse)
			#smoothed beam
			ellipse= Ellipse(xy=(posx,posy),width=b_smooth,height=b_smooth,angle=0,edgecolor='black',fc='none',lw=0.5)
			ax.add_patch(ellipse)
			cbar=plt.colorbar(c,fraction=0.046, pad=0.04)
			cbar.ax.set_ylabel(r'$\chi_0$ from FDF [degrees]')
			plt.gca().set_aspect("equal")
			plt.xlim(xmin,xmax)
			plt.ylim(ymin,ymax)
			#plt.show()
			plt.savefig(fig_directory+src+'_chi0',dpi=400,bbox_inches='tight',transparent=True)
			plt.savefig(directory+src+'_chi0',dpi=400,bbox_inches='tight',transparent=True)

			#I don't know why I need to do the followling line but here we are
			pkrm_fit_masked=np.where(pkrm_fit_masked.mask==1,np.nan,pkrm_fit_masked)
			plt.figure(dpi= 600)
			ax=plt.subplot(projection=wcs)
			#print(np.nanmedian(pkrm_fit_masked))
			#print(np.nanmean(pkrm_fit_masked))
			#print(np.nanstd(pkrm_fit_masked))

			c=ax.imshow(pkrm_fit_masked, origin='lower', cmap=cmap_redblue,vmin=-1*rm_scale,vmax=rm_scale)
			ax.contour(cont_im,levels=i_contours,colors='black',transform=ax.get_transform(cont_wcs))
			plt.xlabel('Right Ascension [J2000]')
			plt.ylabel('Declination [J2000]')
			ellipse= Ellipse(xy=(posx,posy),width=bmaj_pix,height=bmin_pix,angle=bpa+90,edgecolor='black',fc='black',lw=0.5)
			ax.add_patch(ellipse)
			#smoothed beam
			ellipse= Ellipse(xy=(posx,posy),width=b_smooth,height=b_smooth,angle=0,edgecolor='black',fc='none',lw=0.5)
			ax.add_patch(ellipse)
			ax.set_facecolor('grey')
			cbar=plt.colorbar(c,fraction=0.046, pad=0.04)
			cbar.ax.set_ylabel(r'RM from synthesis [rad m$^{-2}$]')
			plt.gca().set_aspect("equal")
			plt.xlim(xmin,xmax)
			plt.ylim(ymin,ymax)
			#plt.show()
			plt.savefig(fig_directory+src+'RM_peak_fitted.png',dpi=400,bbox_inches='tight',transparent=True)
			plt.savefig(directory+src+'RM_peak_fitted.png',dpi=400,bbox_inches='tight',transparent=True)
			'''
			plt.figure(dpi= 600)
			ax=plt.subplot(projection=wcs)
			c=ax.imshow(alpha_arr, origin='lower',cmap=cividis)
			ax.contour(cont_im,levels=i_contours,colors='black',transform=ax.get_transform(cont_wcs))
			plt.xlabel('Right Ascension [J2000]')
			plt.ylabel('Declination [J2000]')
			ellipse= Ellipse(xy=(posx,posy),width=bmaj_pix,height=bmin_pix,angle=bpa+90,edgecolor='black',fc='black',lw=0.5)
			ax.add_patch(ellipse)
			#smoothed beam
			#ellipse= Ellipse(xy=(posx,posy),width=b_smooth,height=b_smooth,angle=0,edgecolor='black',fc='none',lw=0.5)
			#ax.add_patch(ellipse)
			ax.set_facecolor('grey')
			cbar=plt.colorbar(c,fraction=0.046, pad=0.04)
			cbar.ax.set_ylabel(r'In-band spectral index')
			plt.gca().set_aspect("equal")
			plt.xlim(xmin,xmax)
			plt.ylim(ymin,ymax)
			#plt.show()
			plt.savefig(fig_directory+src+'alpha.png',dpi=400,bbox_inches='tight',transparent=True)
			plt.savefig(directory+src+'alpha.png',dpi=400,bbox_inches='tight',transparent=True)
			'''


			#rmplot with pol contours
			plt.figure(dpi= 600)
			ax=plt.subplot(projection=wcs)
			c=ax.imshow(pkrm_fit_masked, origin='lower', cmap=cmap_redblue,vmin=-1*rm_scale,vmax=rm_scale)
			ax.contour(average_p_map,levels=i_contours,colors='black')
			plt.xlabel('Right Ascension [J2000]')
			plt.ylabel('Declination [J2000]')
			ellipse= Ellipse(xy=(posx,posy),width=bmaj_pix,height=bmin_pix,angle=bpa+90,edgecolor='black',fc='black',lw=0.5)
			ax.add_patch(ellipse)
			#smoothed beam
			ellipse= Ellipse(xy=(posx,posy),width=b_smooth,height=b_smooth,angle=0,edgecolor='black',fc='none',lw=0.5)
			ax.add_patch(ellipse)
			cbar=plt.colorbar(c,fraction=0.046, pad=0.04)
			cbar.ax.set_ylabel(r'RM from synthesis [rad m$^{-2}$]')
			plt.gca().set_aspect("equal")
			plt.xlim(xmin,xmax)
			plt.ylim(ymin,ymax)
			#plt.show()
			plt.savefig(fig_directory+src+'RM_peak_fitted_polconts.png',dpi=400,bbox_inches='tight',transparent=True)
			plt.savefig(directory+src+'RM_peak_fitted_polconts.png',dpi=400,bbox_inches='tight',transparent=True)

			#rmplot with pol contours
			plt.figure(dpi= 600)
			ax=plt.subplot(projection=wcs)
			c=ax.imshow(pkrm_fit_masked, origin='lower', cmap=plt.cm.viridis,vmin=np.nanmin(rm_cut),vmax=np.nanmax(rm_cut))
			ax.contour(average_p_map,levels=i_contours,colors='black')
			plt.xlabel('Right Ascension [J2000]')
			plt.ylabel('Declination [J2000]')
			ellipse= Ellipse(xy=(posx,posy),width=bmaj_pix,height=bmin_pix,angle=bpa+90,edgecolor='black',fc='black',lw=0.5)
			ax.add_patch(ellipse)
			#smoothed beam
			ellipse= Ellipse(xy=(posx,posy),width=b_smooth,height=b_smooth,angle=0,edgecolor='black',fc='none',lw=0.5)
			ax.add_patch(ellipse)
			cbar=plt.colorbar(c,fraction=0.046, pad=0.04)
			cbar.ax.set_ylabel(r'RM from synthesis [rad m$^{-2}$]')
			plt.gca().set_aspect("equal")
			plt.xlim(xmin,xmax)
			plt.ylim(ymin,ymax)
			#plt.show()
			plt.savefig(fig_directory+src+'RM_peak_fitted_polconts_viridis.png',dpi=400,bbox_inches='tight',transparent=True)
			plt.savefig(directory+src+'RM_peak_fitted_polconts_viridis.png',dpi=400,bbox_inches='tight',transparent=True)

			#rm serrs
			RMerr_arr_masked=np.where(cont_pol_mask.mask==True,np.nan,RMerr_arr)
			plt.figure(dpi= 600)
			ax=plt.subplot(projection=wcs)
			c=ax.imshow(RMerr_arr_masked, origin='lower', cmap=cividis)
			ax.contour(cont_im,levels=i_contours,colors='black',transform=ax.get_transform(cont_wcs))
			plt.xlabel('Right Ascension [J2000]')
			plt.ylabel('Declination [J2000]')
			ellipse= Ellipse(xy=(posx,posy),width=bmaj_pix,height=bmin_pix,angle=bpa+90,edgecolor='black',fc='black',lw=0.5)
			ax.add_patch(ellipse)
			#smoothed beam
			ellipse= Ellipse(xy=(posx,posy),width=b_smooth,height=b_smooth,angle=0,edgecolor='black',fc='none',lw=0.5)
			ax.add_patch(ellipse)
			cbar=plt.colorbar(c,fraction=0.046, pad=0.04)
			cbar.ax.set_ylabel(r'Error in RM [rad m$^{-2}$]')
			plt.gca().set_aspect("equal")
			plt.xlim(xmin,xmax)
			plt.ylim(ymin,ymax)
			#plt.show()
			plt.savefig(fig_directory+src+'RM_peak_fitted_errs.png',dpi=400,bbox_inches='tight',transparent=True)
			plt.savefig(directory+src+'RM_peak_fitted_errs.png',dpi=400,bbox_inches='tight',transparent=True)

			#polfrac change map
			plt.figure(dpi= 400)
			ax=plt.subplot(projection=wcs)
			c=ax.imshow(polfrac_change_masked, origin='lower', cmap=viridis,vmin=0,vmax=2)
			ax.contour(cont_im,levels=i_contours,colors='black',transform=ax.get_transform(cont_wcs))
			plt.xlabel('Right Ascension [J2000]')
			plt.ylabel('Declination [J2000]')
			ellipse= Ellipse(xy=(posx,posy),width=bmaj_pix,height=bmin_pix,angle=bpa+90,edgecolor='black',fc='gray',lw=0.5)
			ax.add_patch(ellipse)
			#smoothed beam
			ellipse= Ellipse(xy=(posx,posy),width=b_smooth,height=b_smooth,angle=0,edgecolor='black',fc='none',lw=0.5)
			ax.add_patch(ellipse)
			cbar=plt.colorbar(c,fraction=0.046, pad=0.04)
			cbar.ax.set_ylabel(r'Band depolarisation factor')
			plt.gca().set_aspect("equal")
			plt.xlim(xmin,xmax)
			plt.ylim(ymin,ymax)
			#plt.show()
			plt.savefig(fig_directory+src+'polfrac_change.png',dpi=400,bbox_inches='tight',transparent=True)
			plt.savefig(directory+src+'polfrac_change.png',dpi=400,bbox_inches='tight',transparent=True)

			plt.figure(dpi= 400)
			ax=plt.subplot(projection=wcs)
			c=ax.imshow(k_arr, origin='lower', cmap=viridis,vmin=vmax(k_arr,1),vmax=vmax(k_arr,99))
			ax.contour(cont_im,levels=i_contours,colors='black',transform=ax.get_transform(cont_wcs))
			plt.xlabel('Right Ascension [J2000]')
			plt.ylabel('Declination [J2000]')
			ellipse= Ellipse(xy=(posx,posy),width=bmaj_pix,height=bmin_pix,angle=bpa+90,edgecolor='black',fc='gray',lw=0.5)
			ax.add_patch(ellipse)
			#smoothed beam
			ellipse= Ellipse(xy=(posx,posy),width=b_smooth,height=b_smooth,angle=0,edgecolor='black',fc='none',lw=0.5)
			ax.add_patch(ellipse)
			cbar=plt.colorbar(c,fraction=0.046, pad=0.04)
			cbar.ax.set_ylabel(r'Burn law k [rad$^{-4}$ m$^{-4}$]')
			plt.gca().set_aspect("equal")
			plt.xlim(xmin,xmax)
			plt.ylim(ymin,ymax)
			#plt.show()
			plt.savefig(fig_directory+src+'burnk.png',dpi=400,bbox_inches='tight',transparent=True)
			plt.savefig(directory+src+'burnk.png',dpi=400,bbox_inches='tight',transparent=True)

			plt.figure(dpi= 400)
			ax=plt.subplot(projection=wcs)
			c=ax.imshow(burnredchi, origin='lower', cmap=viridis)
			ax.contour(cont_im,levels=i_contours,colors='black',transform=ax.get_transform(cont_wcs))
			plt.xlabel('Right Ascension [J2000]')
			plt.ylabel('Declination [J2000]')
			ellipse= Ellipse(xy=(posx,posy),width=bmaj_pix,height=bmin_pix,angle=bpa+90,edgecolor='black',fc='gray',lw=0.5)
			ax.add_patch(ellipse)
			#smoothed beam
			ellipse= Ellipse(xy=(posx,posy),width=b_smooth,height=b_smooth,angle=0,edgecolor='black',fc='none',lw=0.5)
			ax.add_patch(ellipse)
			cbar=plt.colorbar(c,fraction=0.046, pad=0.04)
			cbar.ax.set_ylabel(r'Reduced $\chi^{2}$ for Burn law fit')
			plt.gca().set_aspect("equal")
			plt.xlim(xmin,xmax)
			plt.ylim(ymin,ymax)
			#plt.show()
			plt.savefig(fig_directory+src+'burnk_chi2.png',dpi=400,bbox_inches='tight',transparent=True)
			plt.savefig(directory+src+'burnk_chi2.png',dpi=400,bbox_inches='tight',transparent=True)

			plt.figure(dpi= 400)
			ax=plt.subplot(projection=wcs)
			c=ax.imshow(cc_stdv, origin='lower', cmap=inferno,vmin=0,vmax=100)
			ax.contour(cont_im,levels=i_contours,colors='black',transform=ax.get_transform(cont_wcs))
			plt.xlabel('Right Ascension [J2000]')
			plt.ylabel('Declination [J2000]')
			ellipse= Ellipse(xy=(posx,posy),width=bmaj_pix,height=bmin_pix,angle=bpa+90,edgecolor='black',fc='gray',lw=0.5)
			ax.add_patch(ellipse)
			#smoothed beam
			ellipse= Ellipse(xy=(posx,posy),width=b_smooth,height=b_smooth,angle=0,edgecolor='black',fc='none',lw=0.5)
			ax.add_patch(ellipse)
			cbar=plt.colorbar(c,fraction=0.046, pad=0.04)
			cbar.ax.set_ylabel(r'Clean componant spread [rad m$^{-2}$]')
			plt.gca().set_aspect("equal")
			plt.xlim(xmin,xmax)
			plt.ylim(ymin,ymax)
			#plt.show()
			plt.savefig(fig_directory+src+'cc_stdv.png',dpi=400,bbox_inches='tight',transparent=True)
			plt.savefig(directory+src+'cc_stdv.png',dpi=400,bbox_inches='tight',transparent=True)		

			fdf_width=np.where(cont_pol_mask.mask==1,np.nan,fdf_width)
			plt.figure(dpi= 400)
			ax=plt.subplot(projection=wcs)
			c=ax.imshow(fdf_width, origin='lower', cmap=inferno,vmin=0,vmax=vmax(fdf_width,99))
			ax.contour(cont_im,levels=i_contours,colors='black',transform=ax.get_transform(cont_wcs))
			plt.xlabel('Right Ascension [J2000]')
			plt.ylabel('Declination [J2000]')
			ellipse= Ellipse(xy=(posx,posy),width=bmaj_pix,height=bmin_pix,angle=bpa+90,edgecolor='black',fc='gray',lw=0.5)
			ax.add_patch(ellipse)
			#smoothed beam
			ellipse= Ellipse(xy=(posx,posy),width=b_smooth,height=b_smooth,angle=0,edgecolor='black',fc='none',lw=0.5)
			ax.add_patch(ellipse)
			cbar=plt.colorbar(c,fraction=0.046, pad=0.04)
			cbar.ax.set_ylabel(r'FDF width [rad m$^{-2}$]')
			plt.gca().set_aspect("equal")
			plt.xlim(xmin,xmax)
			plt.ylim(ymin,ymax)
			#plt.show()
			plt.savefig(fig_directory+src+'fdf_width.png',dpi=400,bbox_inches='tight',transparent=True)
			plt.savefig(directory+src+'fdf_width.png',dpi=400,bbox_inches='tight',transparent=True)

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
				WISEim=glob.glob(WISEfolder+'*-w1-*ra{}*dec{}*'.format(ra_str,dec_str))
				blue_im=WISEim[0]
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
			
			R,Rhead=fitsopen(glob.glob(imagefolder+'DES*_i*')[0])
			G,Ghead=fitsopen(glob.glob(imagefolder+'DES*_r*')[0])
			B,Bhead=fitsopen(glob.glob(imagefolder+'DES*_g*')[0])
			DES_wcs=WCS(Rhead)

			img = np.zeros((R.shape[0], R.shape[1], 3), dtype=float)
			img[:,:,0] = img_scale.linear(R, scale_min=0.1, scale_max=200)
			img[:,:,1] = img_scale.linear(G, scale_min=0.1, scale_max=200)
			img[:,:,2] = img_scale.linear(B, scale_min=0.1, scale_max=200)
			
			plt.figure(dpi=dpi)
			# load the data for the background greyscale 
			greyscale_data,greyscale_header=fitsopen(greyscale_im)
			greyscale_wcs=WCS(greyscale_header)
			pix_scale=greyscale_header['CDELT2']*3600*u.arcsec
			# plot the greyscale
			ax=plt.subplot(projection=greyscale_wcs)
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
				ax.contour(smoothed, transform=ax.get_transform(purple_wcs), colors=[twilightcmap(0.4)],levels=purple_contours,linewidths=0.3)

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
				ax.contour(blue_data, transform=ax.get_transform(blue_wcs), colors=[twilightcmap(0.2)],levels=blue_contours,linewidths=0.3)

				# calculate transparency array for map (alphas) by normalising the data to be between 0 and 1
				#blue_alphas=np.divide(blue_data-blue_thresh,np.nanmax(blue_data)-blue_thresh)
				#blue_alphas=np.where(blue_alphas<0,0,blue_alphas)
				#ax3=plt.imshow(blue_data,origin='lower',transform=ax.get_transform(blue_wcs),cmap=bluecmap,alpha=blue_alphas) 

				if blue_label!=None:
					ax.plot(-100,-100,'-',c=twilightcmap(0.2),label=blue_label,linewidth=3)


			#add EMU data
			if EMUSB!='NaN':
				red_data,red_header=fitsopen(imagefolder+src+'_EMU.fits')
				red_label='ASKAP 944 MHz'
				red_wcs=WCS(red_header)
				red_thresh=3.*EMUrms
				red_contours = [red_thresh * i for i in contourmults]
				ax.contour(red_data, transform=ax.get_transform(red_wcs), colors=[twilightcmap(0.75)],levels=red_contours,linewidths=0.3)
				# calculate transparency array for map (alphas) by normalising the data to be between 0 and 1
				red_alphas=np.divide(red_data-red_thresh,np.nanmax(red_data)-red_thresh)
				red_alphas=np.where(red_alphas<0,0,red_alphas)
				ax5=plt.imshow(red_data,origin='lower',transform=ax.get_transform(red_wcs),cmap=redcmap,alpha=red_alphas)

				if red_label!=None:
					ax.plot(-100,-100,'-',c=twilightcmap(0.75),label=red_label,linewidth=3)

				r_bmaj=red_header['BMAJ']*u.degree
				r_bmin=red_header['BMIN']*u.degree
				r_bpa=red_header['BPA']
				r_bmaj_pix=r_bmaj/pix_scale
				r_bmin_pix=r_bmin/pix_scale

				lowerleft=utils.skycoord_to_pixel(utils.pixel_to_skycoord(xmin,ymin,red_wcs),greyscale_wcs)
				upperright=utils.skycoord_to_pixel(utils.pixel_to_skycoord(xmax,ymax,red_wcs),greyscale_wcs)

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
				yellow_contours = [yellow_thresh * i for i in contourmults]
				ax.contour(yellow_data, transform=ax.get_transform(yellow_wcs), colors=['gold'],levels=yellow_contours,linewidths=0.3)

				# calculate transparency array for map (alphas) by normalising the data to be between 0 and 1
				#yellow_alphas=np.divide(yellow_data-yellow_thresh,np.nanmax(yellow_data)-yellow_thresh)
				#yellow_alphas=np.where(yellow_alphas<0,0,yellow_alphas)
				#ax5=plt.imshow(yellow_data,origin='lower',transform=ax.get_transform(yellow_wcs),cmap=yellowcmap,alpha=yellow_alphas)

				if yellow_label!=None:
					ax.plot(-100,-100,'-',c='gold',label=yellow_label,linewidth=3)

				y_bmaj=yellow_header['BMAJ']*u.degree
				y_bmin=yellow_header['BMIN']*u.degree
				y_bpa=yellow_header['BPA']
				y_bmaj_pix=y_bmaj/pix_scale
				y_bmin_pix=y_bmin/pix_scale	

				lowerleft=utils.skycoord_to_pixel(utils.pixel_to_skycoord(xmin,ymin,yellow_wcs),greyscale_wcs)
				upperright=utils.skycoord_to_pixel(utils.pixel_to_skycoord(xmax,ymax,yellow_wcs),greyscale_wcs)

				greyxmin=lowerleft[0]
				greyymin=lowerleft[1]
				greyxmax=upperright[0]
				greyymax=upperright[1]

				rectangle=Rectangle(xy=(greyxmin+beampad-(0.6*r_bmaj_pix),greyymin+beampad-(0.6*r_bmaj_pix)),width=r_bmaj_pix*1.2,height=r_bmaj_pix*1.2,fc='white',alpha=0.5)
				ax.add_patch(rectangle)
			if 	EMUSB!='NaN':
				ellipse= Ellipse(xy=(greyxmin+beampad,greyymin+beampad),width=r_bmaj_pix,height=r_bmin_pix,angle=r_bpa+90,edgecolor=twilightcmap(0.75),fc='None',lw=1)
				ax.add_patch(ellipse)	
				print('Adding EMU beam')
			if POSSUMSB!='NaN' and POSSUMSB!='10035':
				ellipse= Ellipse(xy=(greyxmin+beampad,greyymin+beampad),width=y_bmaj_pix,height=y_bmin_pix,angle=y_bpa+90,edgecolor='gold',fc='None',lw=1)
				ax.add_patch(ellipse)			
				print('Adding POSSUM beam')

			if z!=None: 
				sb_length_pix= sb_length_arcsec/pix_scale
				scalebar=Rectangle(xy=(greyxmax-sb_length_pix-10,greyymin+3),width=sb_length_pix,height=sb_length_pix/20.,edgecolor='none',fc='white',alpha=1)
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
			#ax.set_xlim(0,greyscale_data.shape[0])
			#ax.set_ylim(0,greyscale_data.shape[1])
			#plt.legend(loc=3)
			plt.savefig('/Volumes/NARNIA/pilot_cutouts/'+src+'_DES.png',dpi=dpi,transparent=True,bbox_inches='tight')
			plt.savefig('/Users/emma/OneDrive/PhD/thesis/Figures/new/'+src+'_mwl.png',dpi=dpi,transparent=True,bbox_inches='tight')


			#figure for Paddy
			print('doing this plot')
			plt.figure(figsize=(13,6))

			SMALL_SIZE = 10
			MEDIUM_SIZE = 11
			BIGGER_SIZE = 12
			plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
			plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
			plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
			plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
			plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
			plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

			'''
			ax1=plt.subplot(131,projection=greyscale_wcs,adjustable='box')
			ax1.imshow(greyscale_data,origin='lower',cmap=greyscale,vmin=vmax(greyscale_data,5),vmax=vmax(greyscale_data,99.9))
			ax1.contour(red_data,levels=i_contours,colors=[twilightcmap(0.75)],transform=ax1.get_transform(cont_wcs))
			ax1.set_xlabel('RA [J2000]')
			ax1.set_ylabel('Dec [J2000]')
			ax1.imshow(red_data,origin='lower',transform=ax1.get_transform(cont_wcs),cmap=redcmap,alpha=red_alphas)

			x,y=utils.skycoord_to_pixel(utils.pixel_to_skycoord(examplex,exampley,wcs),greyscale_wcs)
			ax1.scatter(x,y,marker='s',facecolor='none',edgecolor='yellow')
			ax1.set_xlim(0.28*greyscale_data.shape[0],0.62*greyscale_data.shape[0])
			ax1.set_ylim(0,0.6*greyscale_data.shape[1])
			'''
			ax1=plt.subplot(131,projection=wcs,adjustable='box')
			c=ax1.imshow(average_p_map, origin='lower', cmap=magmacmap,vmin=0,vmax=vmax(average_p_map_masked,99.99))#norm=colors.LogNorm(vmin=0.001,vmax=0.01))#,vmin=0,vmax=vmax(average_p_map_masked,99.99))
			q=ax1.quiver(x_pos,y_pos,x_reduced,y_reduced,color='white',scale=50,width=0.002,headlength=0,headwidth=0)
			ax1.contour(cont_im,levels=i_contours,colors='grey',transform=ax1.get_transform(cont_wcs))
			ax1.set_xlabel('Right Ascension [J2000]')
			ax1.set_ylabel('Declination [J2000]')
			ellipse= Ellipse(xy=(xmin+25+bmaj_pix,ymin+60+bmaj_pix),width=bmaj_pix,height=bmin_pix,angle=bpa+90,edgecolor='gray',fc='gray',lw=0.5)
			ax1.add_patch(ellipse)
			cbar=plt.colorbar(c,fraction=0.089, pad=0.005)
			cbar.ax.set_ylabel('Polarised intensity [Jy/bm]')
			ax1.scatter(examplex,exampley,marker='s',s=6,facecolor='none',edgecolor='black')

			ax1.set_xlim(xmin,xmax)
			ax1.set_ylim(ymin,ymax)


			ax2=plt.subplot(232,projection=wcs,adjustable='box')
			ax2im=ax2.imshow(pkrm_fit_masked,origin='lower',cmap=cmap_redblue,vmin=-1*rm_scale,vmax=rm_scale)
			ax2.contour(cont_im,levels=i_contours,colors='black',transform=ax2.get_transform(cont_wcs))
			ax2.set_xlim(examplex-examplesquare,examplex+examplesquare)
			ax2.set_ylim(exampley-examplesquare,exampley+examplesquare)
			ax2.scatter(examplex,exampley,marker='s',s=10,facecolor='none',edgecolor='black')

			ax2.set_xlabel(' ')
			ax2.set_ylabel(' ')
			#ax2.set_ylabel(' ')
			cbar2=plt.colorbar(ax2im,fraction=0.0475, pad=0.005)
			cbar2.ax.set_ylabel(r'Peak RM [rad m$^{-2}$]')

			ax3=plt.subplot(233,adjustable='box')
			ax3.errorbar(lamdasq,slice_p,QUrms,linestyle="None",color=twilight(0.5),alpha=0.5)
			ax3.errorbar(lamdasq,slice_q,QUrms,linestyle="None",color=twilight(0.2),alpha=0.5)
			ax3.errorbar(lamdasq,slice_u,QUrms,linestyle="None",color=twilight(0.7),alpha=0.5)
			ax3.scatter(lamdasq,slice_p,label='P',marker='.',color=twilight(0.5),s=10)
			ax3.scatter(lamdasq,slice_q,label='Q',marker='+',linewidth=1,color=twilight(0.2),s=10)
			ax3.scatter(lamdasq,slice_u,label='U',marker='x',linewidth=1,color=twilight(0.7),s=10)
			ax3.axhline(y=0,color='k')
			ax3.set_xlim(np.nanmin(lamdasq),np.nanmax(lamdasq))
			ax3.set_xlabel(r'$\lambda^2$ [m$^2$]')
			ax3.set_ylabel('Intensity [Jy/bm]')
			ax3.legend()

			ax4=plt.subplot(235,adjustable='box')
			ax4im=ax4.scatter(slice_u,slice_q,c=lamdasq,cmap=plt.cm.plasma,s=10)
			ax4.errorbar(slice_u,slice_q,xerr=U_cube_rms,yerr=Q_cube_rms,linestyle='none',color='k')
			ax4.set_xlabel('U [Jy/bm]')
			ax4.set_ylabel(' ')
			ax4.set_xlim(umin,umax)
			ax4.set_ylim(qmin,qmax)
			ax4.axhline(y=0,color='k')
			ax4.axvline(x=0,color='k')
			#ax5.set_aspect('equal')
			cbar4=plt.colorbar(ax5im,fraction=0.048, pad=0.005)
			cbar4.ax.set_ylabel(r'$\lambda^2$ [m$^2$]')


			ax5=plt.subplot(236,adjustable='box')
			ax5.plot(rmarray,fdf_clean[:,exampley,examplex],label='Cleaned',linewidth=1,color=twilight(0.5))
			ax5.plot(rmarray,fdfdata[:,exampley,examplex],label='Dirty',linestyle='dashed',linewidth=0.7,color=twilight(0.15))
			ax5.vlines(clean_lines,clean_lines_min,clean_lines_max,color=twilight(0.8),linewidth=1,label='CC')
			ax5.set_ylabel(r'|P| [mJy beam$^{-1}$ rmsf$^{-1}$]')
			ax5.set_xlabel(r'$\phi$  [rad m$^{-2}]$')
			ax5.set_xlim(-500,500)
			ax5.axvline(x=galacticRM,color=twilight(0.3),linestyle=':',linewidth=1,label='Galactic')
			ax5.legend()
			try:
				ymaxfdf=1.05*np.nanmax(fdfdata[:,exampley,examplex])
				ax5.set_ylim(0,ymaxfdf)
				ax5.fill_betweenx([0,ymaxfdf],x1=galacticRM-gRMerr, x2=galacticRM+gRMerr,color=twilight(0.3),alpha=0.5)
			except:
				print("uh oh")

			plt.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.95, wspace=0.4, hspace=0.2)
			plt.savefig('/Users/emma/OneDrive/PhD/thesis/Figures/new/'+src+'_compactplts.png',dpi=dpi,transparent=True,bbox_inches='tight')

			break
		#except:
			#print("Something went wrong with this source!")
###############################################################################
if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("-f","--filename",help='Name of file with source list',
		default=defaultfile)
	args = ap.parse_args()
	main(args)






