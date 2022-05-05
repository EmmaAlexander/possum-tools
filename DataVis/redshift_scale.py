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

def main(args):
	###########################################################################
	#read in sources
	print('Reading in source list')
	sources=np.loadtxt(args.filename,dtype='str')
	nsrc=sources.shape[0]
	###########################################################################
	sourcecount=0
	for i in range(0,nsrc):
		src=sources[i,0]
		LAS=float(sources[i,5]) #arcmin
		z=float(sources[i,8]) 

		if np.isfinite(z):
			dist = Planck18.comoving_distance(z) #distance in MPc
			scale=Planck18.arcsec_per_kpc_comoving(z) # scale in arcsec per kpc
			source_length=(LAS*60*u.arcsec)/scale
			print(src,dist,source_length)

if __name__ == "__main__":
	ap = argparse.ArgumentParser()
	ap.add_argument("-f","--filename",help='Name of file with source list',
		default=defaultfile)
	args = ap.parse_args()
	main(args)




