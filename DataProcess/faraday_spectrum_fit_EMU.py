from functions import *
import argparse
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np 
import warnings
warnings.filterwarnings('ignore')
import os
import matplotlib.pyplot as plt

rootdirectory='/Volumes/NARNIA/pilot_cutouts/'
defaultfile='/Users/emma/GitHub/possum-tools/DataProcess/emu_cube_sources.txt'

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

		ithresh=3*emu_rms
		#check that RMsynthesis files are in the directory

		#if os.path.isfile(directory+'FDF_chi0_fitted.fits'):
			#already run for this source
		#	print("Skipping: FDF fitted files already exist for source {}".format(src)) 
		if os.path.isfile(directory+'EMU_FDF_maxPI.fits') == False:
			#no FDF for soure
			print("Skipping: No FDF files for source {}".format(src))
		else:
			#Read in data
			print("Reading in the fits files")



			#EMU
			PI_im,PI_header=fitsopen(directory+'EMU_FDF_maxPI.fits')
			#peak rm from RMsynth
			pkrm_im,pkrm_header=fitsopen(directory+'EMU_FDF_peakRM.fits')
			pkrm_im=pkrm_im[0,:,:]
			#FDF cubes
			fdfdata,fdfheader=fitsopen(directory+'EMU_FDF_tot_dirty.fits')
			fdf_im,fdfimhead=fitsopen(directory+'EMU_FDF_im_dirty.fits')
			fdf_real,fdfrealhead=fitsopen(directory+'EMU_FDF_real_dirty.fits')

			fdf_ref_rm=fdfheader['CRVAL3']
			fdf_ref_chan=fdfheader['CRPIX3']
			fdfres=fdfheader['CDELT3']
			nrmchans=fdfheader['NAXIS3']
			fdfmin=fdf_ref_rm-(fdf_ref_chan*fdfres) + fdfres
			rmarray=np.arange(fdfmin,fdfmin+(fdfres*nrmchans),fdfres)

			#do the FDF fitting
			print("Fitting the FDF per pixel")
			rmsf_rm_fit,rmsf_chi0_fit,pkPI_fit=fdf_fit(pkrm_im,fdfdata,fdf_real,fdf_im,rmarray,rmthresh=1999)
			#save to file
			print("Saving the fits")
			hdu = fits.PrimaryHDU(rmsf_rm_fit,header=pkrm_header)
			hdu.writeto(directory+'EMU_FDF_peakRM_fitted.fits',overwrite=True)
			hdu2 = fits.PrimaryHDU(rmsf_chi0_fit,header=pkrm_header)
			hdu2.writeto(directory+'EMU_FDF_chi0_fitted.fits',overwrite=True)

			hdu4 = fits.PrimaryHDU(pkPI_fit,header=PI_header)
			hdu4.writeto(directory+'FDF_pkPI_fitted_CORRECT.fits',overwrite=True)


			#rmsf_rm_fit_width=fdf_fit_gauss(pkrm_im,fdfdata,fdf_real,fdf_im,rmarray,rmthresh=1995)

			#hdu3 = fits.PrimaryHDU(rmsf_rm_fit_width,header=pkrm_header)
			#hdu3.writeto(directory+'FDF_peakRM_fitted_width.fits',overwrite=True)
				
###############################################################################
if __name__ == "__main__":

	ap = argparse.ArgumentParser()
	ap.add_argument("-f","--filename",help='Name of file with source list',
		default=defaultfile)
	args = ap.parse_args()

	main(args)





















