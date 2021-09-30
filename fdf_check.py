from functions import fit, peakfit, quadfit, image_rms, cube_rms,fitsopen, finitelims,vmax,get_rm_scale
import argparse
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.patches import Ellipse
import warnings
warnings.filterwarnings('ignore')
import os
from os import path

rootdirectory='/Volumes/TARDIS/Work/askap/'
defaultfile=rootdirectory+'final_sourcelist.txt'

def main():
	###########################################################################
	#read in sources
	print('Reading in source list')
	sources=np.loadtxt(defaultfile,dtype='str')
	nsrc=sources.shape[0]
	###########################################################################
	sourcecount=0
	for src in sources[:,0]:
		sourcecount+=1
		print('Looking at source {} ({}/{})'.format(src,sourcecount,nsrc))
		directory=rootdirectory+src+'/'
		#Read in data
		print("Reading in the fits files")

		#fitted peak rm 
		try:
			pkrm_fit,pkrm_fit_header=fitsopen(directory+'FDF_peakRM_fitted.fits')
			median=np.nanmedian(pkrm_fit)
			
			if median>=1000:
				#somethign has gone wrong
				print("uh oh median is {}".format(median))
				newfolder=directory+'old_fdf_fits'
				os.system("mkdir {}".format(newfolder))
				os.system("mv {} {}".format(directory+'FDF_peakRM_fitted.fits',newfolder+'/'))
				os.system("mv {} {}".format(directory+'FDF_chi0_fitted.fits',newfolder+'/'))
				print("moved")
		except:
			print("file not found")



if __name__ == "__main__":
	main()