import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import cmasher as cmr
import astropy.units as u
import astropy.coordinates as coord
from astropy.io import ascii
from astropy.io import fits
from astropy.wcs import WCS
from functions import *


plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'cmu serif'

SMALL_SIZE = 8
MEDIUM_SIZE = 8
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=MEDIUM_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize

fig_directory='/Users/emma/OneDrive/PhD/thesis/Figures/'
cmap_blue = cmr.get_sub_cmap('twilight_shifted', 0, 0.5)
cmap_red = cmr.get_sub_cmap('twilight', 0.5, 1)
cmap_redblue=cmr.get_sub_cmap('twilight_shifted', 0.1, 0.9)
cmap=plt.cm.twilight_shifted

def main():


	faradaysky,header=fitsopen('/Volumes/TARDIS/Work/askap/Faraday_cutout_pilot.fits')
	faradayuncertainty,header2=fitsopen('/Volumes/TARDIS/Work/askap/Faraday_error_pilot.fits')

	print(faradaysky.shape)
	wcs=WCS(header)

	sources=np.loadtxt('source_coords.txt',dtype='str')

	plt.figure()
	ax=plt.subplot(projection=wcs)
	c=ax.imshow(faradaysky, origin='lower', cmap=cmap_redblue,vmin=-50,vmax=50)
	cbar=plt.colorbar(c,fraction=0.046, pad=0.04)

	for i in range(0,sources.shape[0]):
		ra_ha=coord.Angle(sources[i,0],unit=u.hourangle)
		ra = coord.Angle(ra_ha,unit=u.degree)
		dec = coord.Angle(sources[i,1],unit=u.degree)
		coords=coord.SkyCoord(ra=ra,dec=dec)
		pixcoords=wcs.world_to_pixel(coords)
		x=int(round(float(pixcoords[0])))
		y=int(round(float(pixcoords[1])))

		plt.scatter(pixcoords[0],pixcoords[1],marker='.',color='k')

		RM=faradaysky[y,x]
		RMerr=faradayuncertainty[y,x]

		print(sources[i,0],sources[i,1])
		print('{} +/- {}'.format(RM,RMerr))

	plt.show()
if __name__ == "__main__":
	main()