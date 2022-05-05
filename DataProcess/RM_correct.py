#because I am stupid 
#and made a mistake when making some other files
#from functions import fit, peakfit, quadfit, image_rms, cube_rms,fitsopen, finitelims,vmax,get_rm_scale,burnkfit,weighted_sample_avg_std
from functions import *
import argparse
from astropy.io import fits
import numpy as np 
import warnings
warnings.filterwarnings('ignore')
import os

rootdirectory='/Volumes/TARDIS/Work/askap/leakage_corrected/'
defaultfile='/Users/emma/GitHub/possum-tools/DataProcess/pilot_sources.txt'

def main(args):
	###########################################################################
	#read in sources
	print('Reading in source list')
	sources=np.loadtxt(args.filename,dtype='str')
	nsrc=sources.shape[0]
	###########################################################################
	for i in range(0,1):
		src=sources[i,0]

		print('Working on source {} ({}/{})'.format(src,i+1,nsrc))
		directory=rootdirectory+src+'/'
		try:
			pkrm_im,pkrm_header=fitsopen(directory+src+'_FDF_peakRM.fits')
			pkrm_fit,pkrm_fit_header=fitsopen(directory+src+'_FDF_peakRM_fitted.fits')

			diff=np.subtract(pkrm_fit,pkrm_im)
			print(np.nanmean(np.abs(diff)))

			pkrm_fit_corrected=np.add(pkrm_fit,67.5)#67.5 for 2018, 9 otherwise
			hdu = fits.PrimaryHDU(pkrm_fit_corrected,header=pkrm_header)
			hdu.writeto(directory+'FDF_peakRM_fitted_corrected.fits',overwrite=True)
		except:
			print("oops")


if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("-f","--filename",help='Name of file with source list',
		default=defaultfile)
	args = ap.parse_args()
	main(args)


