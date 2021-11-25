import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import healpy as hp
import cmasher as cmr
import astropy.units as u
import astropy.coordinates as coord
from astropy.io import ascii
from astropy.wcs import WCS


plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'cmu serif'
plt.rc('legend', fontsize=14)    

fig_directory='/Users/emma/OneDrive/PhD/thesis/Figures/'
cmap_blue = cmr.get_sub_cmap('twilight_shifted', 0, 0.5)
cmap_red = cmr.get_sub_cmap('twilight', 0.5, 1)
cmap=plt.cm.twilight

def main():


	fullsky='/Volumes/TARDIS/Work/askap/spass_dr1_1902_healpix_Tb.i.fits'
	data = hp.read_map(fullsky)
	#plt.subplot(111, projection="mollweide")
	hp.visufunc.mollview(data,min=0.1,max=10,cmap=cmap_blue,title=' ',norm='log')
	hp.visufunc.graticule()


	fields = ascii.read("pilotcoords.txt")
	ra_fields=coord.Angle(fields['RA'],unit=u.hourangle)
	ra = coord.Angle(ra_fields,unit=u.degree)
	dec = coord.Angle(fields['Dec'],unit=u.degree)
	boxsize=5.5*u.degree

	#add on pilot fields
	pitchangle=0.75
	d=pitchangle
	
	for k in range(0,len(ra)):
		#these are the centre coordinates of beam zero for closepack
		x0=(ra[k]/u.degree) + (3*d - d/(2*np.sqrt(2)))
		y0=(dec[k]/u.degree) + (-1*d/(2*np.sqrt(2)) - np.sqrt(3)*d)

		for j in range(0,6):
			for i in range(0,6):
				y=y0 + j*d*(np.sqrt(3)/2)
				if j==0 or j==2 or j==4:
					x=x0 - 2*i*d 
				else: 
					x=x0 - 2*(i+0.5)*d 
				
				hp.visufunc.projscatter(x,y,lonlat='True',marker='o',s=2,color=cmap(0.65),coord='EG')
	hp.visufunc.projscatter(np.nanmean(ra),np.nanmean(dec),lonlat='True',marker='o',s=14,color=cmap(0.65),coord='EG',label='POSSUM Pilot')		

	#now do DRAGN field
	DRAGN_ra_h=coord.Angle('17h04m58s',unit=u.hourangle)
	DRAGN_ra=coord.Angle(DRAGN_ra_h,unit=u.degree)/u.degree
	DRAGN_dec=coord.Angle('-24:40:48',unit=u.degree)/u.degree

	x0=DRAGN_ra+2.5*d
	y0=DRAGN_dec-2.5*d


	for j in range(0,6):
		for i in range(0,6):
			x=x0-(i*d)
			y=y0+(j*d)
			hp.visufunc.projscatter(x,y,lonlat='True',marker='o',s=2,color=cmap(0.85),coord='EG')	
	hp.visufunc.projscatter(DRAGN_ra,DRAGN_dec,lonlat='True',marker='o',s=14,color=cmap(0.85),coord='EG',label='Ophiuchus field')	
			
	plt.legend()
	plt.savefig(fig_directory+'footprints.png',transparent='True',dpi=800)
	plt.show()	


if __name__ == "__main__":
	main()