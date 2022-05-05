from functions import *
import argparse
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np 
import warnings
warnings.filterwarnings('ignore')
import os
import matplotlib.pyplot as plt

rootdirectory='/Volumes/TARDIS/Work/askap/'
defaultfile='/Users/emma/GitHub/possum-tools/DataProcess/Oph_sources.txt'

def main(args):
	###########################################################################
	#read in sources
	print('Reading in source list')
	sources=np.loadtxt(args.filename,dtype='str')
	nsrc=sources.shape[0]
	###########################################################################
	for i in range (0,nsrc):
		src=sources[i,0]
		print('Working on source {} ({}/{})'.format(src,i,nsrc))
		directory=rootdirectory+src+'/'
		#Read in values form source list
		POSSUMSB=sources[i,1]
		EMUSB=sources[i,2]
		emu_rms=float(sources[i,9])
		possum_rms=float(sources[i,10])

		ithresh=3*possum_rms
		#check that RMsynthesis files are in the directory

		#if os.path.isfile(directory+'FDF_chi0_fitted.fits'):
			#already run for this source
		#	print("Skipping: FDF fitted files already exist for source {}".format(src)) 
		if os.path.isfile(directory+'FDF_maxPI.fits') == False:
			#no FDF for soure
			print("Skipping: No FDF files for source {}".format(src))
		else:
			#Read in data
			print("Reading in the fits files")
			#continuum image
			#cont_im,cont_im_header=fitsopen(directory+src+'_POSSUM.fits')

			#wcs=WCS(cont_im_header)
			#peak polarised intensity from RMsynth
			PI_im,PI_header=fitsopen(directory+'FDF_maxPI.fits')
			#peak rm from RMsynth
			pkrm_im,pkrm_header=fitsopen(directory+'FDF_peakRM.fits')
			pkrm_im=pkrm_im[0,:,:]
			#FDF cubes
			fdfdata,fdfheader=fitsopen(directory+'FDF_tot_dirty.fits')
			fdf_im,fdfimhead=fitsopen(directory+'FDF_im_dirty.fits')
			fdf_real,fdfrealhead=fitsopen(directory+'FDF_real_dirty.fits')

			#if(cont_im.shape[0])!=fdfdata.shape[1] or (cont_im.shape[1])!=fdfdata.shape[2]:
			#	print("Skipping: file dimensions do not match. Check files.") 
			#	print(cont_im.shape)
			#	print(fdfdata.shape)
			#else:
				###########################################################################
			fdf_ref_rm=fdfheader['CRVAL3']
			fdf_ref_chan=fdfheader['CRPIX3']
			fdfres=fdfheader['CDELT3']
			nrmchans=fdfheader['NAXIS3']
			fdfmin=fdf_ref_rm-(fdf_ref_chan*fdfres) +fdfres
			rmarray=np.arange(fdfmin,fdfmin+(fdfres*nrmchans),fdfres)

			print(np.nanmin(rmarray))
			print(np.nanmax(rmarray))
			#do the FDF fitting
			print("Fitting the FDF per pixel")

			rmsf_rm_fit,rmsf_chi0_fit,pkPI_fit=fdf_fit(pkrm_im,fdfdata,fdf_real,fdf_im,rmarray,rmthresh=1995)
			#save to file
			print("Saving the fits")
			hdu = fits.PrimaryHDU(rmsf_rm_fit,header=pkrm_header)
			hdu.writeto(directory+'FDF_peakRM_fitted_CORRECT.fits',overwrite=True)

			hdu2 = fits.PrimaryHDU(rmsf_chi0_fit,header=pkrm_header)
			hdu2.writeto(directory+'FDF_chi0_fitted_CORRECT.fits',overwrite=True)
			#rmsf_rm_fit_width=fdf_fit_gauss(pkrm_im,fdfdata,fdf_real,fdf_im,rmarray,rmthresh=1995)

			#hdu3 = fits.PrimaryHDU(rmsf_rm_fit_width,header=pkrm_header)
			#hdu3.writeto(directory+'FDF_peakRM_fitted_width.fits',overwrite=True)

			hdu4 = fits.PrimaryHDU(pkPI_fit,header=PI_header)
			hdu4.writeto(directory+'FDF_pkPI_fitted_CORRECT.fits',overwrite=True)
				
###############################################################################
if __name__ == "__main__":

	ap = argparse.ArgumentParser()
	ap.add_argument("-f","--filename",help='Name of file with source list',
		default=defaultfile)
	args = ap.parse_args()

	main(args)





















