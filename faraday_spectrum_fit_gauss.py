from functions import fit, peakfit, quadfit, image_rms, cube_rms, fitsopen, fdf_fit_gauss,gauss ,gaussfit
import argparse
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np 
import warnings
warnings.filterwarnings('ignore')
import os
import matplotlib.pyplot as plt

rootdirectory='/Volumes/TARDIS/Work/askap/'
defaultfile=rootdirectory+'final_sourcelist2.txt'

def main(args):
	###########################################################################
	#read in sources
	print('Reading in source list')
	sources=np.loadtxt(args.filename,dtype='str')
	nsrc=sources.shape[0]
	###########################################################################
	sourcecount=0
	for src in sources[:,0]:
		sourcecount+=1
		print('Working on source {} ({}/{})'.format(src,sourcecount,nsrc))
		directory=rootdirectory+src+'/'
		#Read in values form source list
		listpos=sourcecount-1
		POSSUMSB=sources[listpos,1]
		EMUSB=sources[listpos,2]
		emu_rms=float(sources[listpos,3])
		possum_rms=float(sources[listpos,4])

		if src[-1] =='U':
			print("This is the EMU source")
			ithresh=3*emu_rms
		else:
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
			try:
				cont_im,cont_im_header=fitsopen(directory+src+'_cont.fits')
			except:
				if src[-1] =='U':
					cont_im,cont_im_header=fitsopen(directory+src+'_EMU.fits')
				else:
					cont_im,cont_im_header=fitsopen(directory+src+'_POSSUM.fits')
			wcs=WCS(cont_im_header)
			#peak polarised intensity from RMsynth
			PI_im,PI_header=fitsopen(directory+'FDF_maxPI.fits')
			#peak rm from RMsynth
			pkrm_im,pkrm_header=fitsopen(directory+'FDF_peakRM.fits')
			pkrm_im=pkrm_im[0,:,:]
			#FDF cubes
			fdfdata,fdfheader=fitsopen(directory+'FDF_clean_tot.fits')
			fdf_im,fdfimhead=fitsopen(directory+'FDF_clean_im.fits')
			fdf_real,fdfrealhead=fitsopen(directory+'FDF_clean_real.fits')

			if(cont_im.shape[0])!=fdfdata.shape[1] or (cont_im.shape[1])!=fdfdata.shape[2]:
				print("Skipping: file dimensions do not match. Check files.") 
				print(cont_im.shape)
				print(fdfdata.shape)
			else:
				###########################################################################
				fdf_ref_rm=fdfheader['CRVAL3']
				fdf_ref_chan=fdfheader['CRPIX3']
				fdfres=fdfheader['CDELT3']
				nrmchans=fdfheader['NAXIS3']
				fdfmin=fdf_ref_rm-(fdf_ref_chan*fdfres) + 1
				rmarray=np.arange(fdfmin,fdfmin+(fdfres*nrmchans),fdfres)
				#do the FDF fitting
				print("Fitting the FDF per pixel")
				rmsf_rm_fit_width=fdf_fit_gauss(pkrm_im,fdfdata,fdf_real,fdf_im,rmarray,cont_im,rmthresh=495,contthresh=ithresh)
	

				#save to file
				print("Saving the fits")
				hdu = fits.PrimaryHDU(rmsf_rm_fit_width,header=pkrm_header)
				hdu.writeto(directory+'FDF_peakRM_fitted_width.fits',overwrite=True)

				#hdu2 = fits.PrimaryHDU(rmsf_chi0_fit,header=pkrm_header)
				#hdu2.writeto(directory+'FDF_chi0_fitted.fits',overwrite=True)


###############################################################################
if __name__ == "__main__":

	ap = argparse.ArgumentParser()
	ap.add_argument("-f","--filename",help='Name of file with source list',
		default=defaultfile)
	args = ap.parse_args()

	main(args)





















