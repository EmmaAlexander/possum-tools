import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits
import math as mth
import warnings
warnings.filterwarnings('ignore')
from scipy import spatial
from functions import fitsopen


rmcmap=plt.cm.seismic
rmcmap.set_bad('gray',1)

rootdirectory='/Volumes/TARDIS/Work/askap/'
defaultfile=rootdirectory+'possum_plots_sourcelist.txt'
sources=np.loadtxt(defaultfile,dtype='str')
nsrc=sources.shape[0]

rmerr_thresh=20


sourcecount=0
for src in sources[:,0]:
	sourcecount+=1
	print('Working on source {} ({}/{})'.format(src,sourcecount,nsrc))
	directory=rootdirectory+src+'/'

	#saved fitted RM array
	hdu=fits.open(directory+'FDF_peakRM_fitted.fits')
	data=hdu[0].data
	hdu.close()

	#load the original fits file to get some details out automatically
	hdu=fits.open(directory+'FDF_peakRM.fits')
	header=hdu[0].header
	hdu.close()

	#load the RM errors file
	RMerrs,RMerrshead=fitsopen(directory+'RM_error.fits')

	#beam
	pixsize=header['CDELT2']
	bmaj=header['BMAJ']
	bmin=header['BMIN']
	bpa=header['BPA']
	bmaj_pix=bmaj/pixsize
	bmin_pix=bmin/pixsize

	samplepix=bmaj_pix
	npix_beam=np.pi*bmaj_pix*bmin_pix/4.

	y_size=data.shape[0]
	x_size=data.shape[1]

	max_samples=np.int(np.floor(np.divide(np.max([y_size,x_size]),samplepix)))

	SF=np.zeros(shape=(5,max_samples),dtype=float)
	#SF_randomised=np.zeros(shape=(2,max_samples),dtype=float)

	#make high error RMs into NaNs
	data=np.where(RMerrs<rmerr_thresh,data,np.nan)

	#plt.imshow(data,origin='lower')
	#plt.show()
	#break

	pix=np.argwhere(np.isfinite(data))
	print(pix.shape)

	#randomise the RM data
	smooshed_data=np.reshape(data[np.isfinite(data)],-1)
	shuffled = sorted(smooshed_data, key=lambda k: np.random.random())
	shuffled_map=np.nan*np.ones(shape=data.shape)
	for j in range(0,len(pix)):
		m=pix[j]
		shuffled_map[m[0],m[1]]=shuffled[j]

	d = spatial.distance.cdist(pix,pix,'euclidean')

	print(max_samples)

	for i in range(1,max_samples):
		print(i)
		dtheta_min=(i-1)*samplepix
		dtheta_max=i*samplepix
		SF[0,i]=(dtheta_min+dtheta_max)/2.
		#SF_randomised[0,i]=(dtheta_min+dtheta_max)/2.
		#SF[3,i]=(dtheta_min+dtheta_max)/2.
		numpoints=0

		values=[]
		value_errs=[]
		#values_rand=[]
		biases=[]
		for j in range(0,len(pix)):
			#get the index of this pixel 
			m=pix[j]
			pointA=data[m[0],m[1]]
			pointA_err=RMerrs[m[0],m[1]]
			#pointA_rand=shuffled_map[m[0],m[1]]
			q=d[:,j]
			roi=np.where((q>=dtheta_min),1,0)
			roi=np.where((q<=dtheta_max),roi,0)
			n_list=pix[np.where(roi==1)]

			for n in n_list:
				numpoints+=1
				pointB=data[n[0],n[1]]
				pointB_err=RMerrs[n[0],n[1]]
				#pointB_rand=shuffled_map[n[0],n[1]]

				value=np.square(np.subtract(pointA,pointB))
				bias=np.square(np.subtract(pointA_err,pointB_err))
				value_err=np.sqrt(2)*np.sqrt(np.square(pointA_err)+np.square(pointB_err))
				#value_rand=np.square(np.subtract(pointA_rand,pointB_rand))

				values.append(value)
				value_errs.append(value_err)
				biases.append(bias)
				#values_rand.append(values_rand)

		if numpoints ==0:
			SF[1,i]=np.nan
			SF[2,i]=0
			SF[3,i]=np.nan
			SF[4,i]=np.nan
		else:
			values=np.asarray(values)
			value_errs=np.asarray(value_errs)
			biases=np.asarray(biases)
			#values_rand=np.asarray(values_rand)
			unbiased_values=np.subtract(values,biases)
			w=1./np.square(value_errs)
			weightedmean=np.divide(np.sum(np.multiply(values,w)),np.sum(w))
			errorweightedmean=1./np.sum(w)
			weighted_errors=np.divide(np.sum(np.multiply(value_errs,w)),np.sum(w))

			weightedunbiased=np.divide(np.sum(np.multiply(unbiased_values,w)),np.sum(w))

			print(weightedmean)
			print(weighted_errors)
			print(weightedunbiased)

			#randommean=np.divide(np.sum(np.multiply(values_rand,w)),np.sum(w))
			
			SF[1,i]=weightedmean
			SF[2,i]=weighted_errors
			SF[3,i]=numpoints
			SF[4,i]=weightedunbiased


	filename=directory+src+'SF_witherrs_new.txt'
	np.savetxt(filename,SF)


