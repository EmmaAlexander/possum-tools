import numpy as np 
import matplotlib.pyplot as plt 
from astropy.io import fits
import math as mth
import warnings
warnings.filterwarnings('ignore')
from scipy import spatial
from functions import fitsopen

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
rmcmap=plt.cm.seismic
rmcmap.set_bad('gray',1)

rootdirectory='/Volumes/TARDIS/Work/askap/leakage_corrected/'
defaultfile='/Users/emma/GitHub/possum-tools/DataProcess/pilot_sources.txt'
sources=np.loadtxt(defaultfile,dtype='str')
nsrc=sources.shape[0]

rmerr_thresh=100
rmval_thresh=500



sourcecount=0
for j in range(0,1):
	if 1+1==2:
		src=sources[j,0]
		xmin=50
		ymin=50
		xmax=450
		ymax=700
		LAS=float(sources[j,5])*60. #get into arcsec

		print('Working on source {} ({}/{})'.format(src,j+1,nsrc))
		directory=rootdirectory+src+'/'

		E_backflow=np.loadtxt(directory+'2018-5540_region_E_backflow.txt')
		W_backflow=np.loadtxt(directory+'2018-5540_region_W_backflow.txt')
		core=np.loadtxt(directory+'2018-5540_region_core.txt')
		N_lobe=np.loadtxt(directory+'2018-5540_region_N_lobe.txt')
		S_lobe=np.loadtxt(directory+'2018-5540_region_S_lobe.txt')

		regions=sum([E_backflow,W_backflow,core,N_lobe,S_lobe])
		regions=regions[79:817,181:682]


		#saved fitted RM array
		hdu=fits.open(directory+'2018-5540_pkrm_forSF.fits')
		data=hdu[0].data
		hdu.close()

		#load the original fits file to get some details out automatically
		hdu=fits.open(directory+'2018-5540_pkrm_dropdeg_forSF.fits')
		header=hdu[0].header
		hdu.close()

		#load the RM errors file
		RMerrs,RMerrshead=fitsopen(directory+src+'_RM_error_forSF.fits')

		#beam
		pixsize=header['CDELT2']
		bmaj=header['BMAJ']
		bmin=header['BMIN']
		bpa=header['BPA']
		bmaj_pix=bmaj/pixsize
		bmin_pix=bmin/pixsize

		samplepix=bmaj_pix
		npix_beam=np.pi*bmaj_pix*bmin_pix/4.

		y_size=ymax-ymin
		x_size=xmax-xmin

		LAS_pix=LAS/pixsize

		max_samples=np.int(np.floor(np.divide(np.max([y_size,x_size]),samplepix)))

		SF=np.zeros(shape=(4,max_samples),dtype=float)
		#SF_randomised=np.zeros(shape=(2,max_samples),dtype=float)

		#make high error RMs into NaNs
		data=np.where(RMerrs<rmerr_thresh,data,np.nan)
		data=np.where(np.abs(data)<rmval_thresh,data,np.nan)
		data=np.where(regions>=1,data,np.nan)

		data[0:ymin,:]=np.nan
		data[ymax:,:]=np.nan
		data[:,0:xmin]=np.nan
		data[:,xmax:]=np.nan
		#plt.imshow(data,origin='lower')
		#plt.show()
		#break

		'''
		print(np.nanmean(data))
		print(np.nanstd(data))
		print(np.nanmin(data))
		print(np.nanmax(data))
		bins=np.arange(start=-2000,stop=2000,step=50)
		plt.hist(np.ndarray.flatten(data),bins=bins)
		plt.show()
		break
		'''
		print("got here")

		pix=np.argwhere(np.isfinite(data))

		print("got here")

		#randomise the RM data
		'''
		smooshed_data=np.reshape(data[np.isfinite(data)],-1)
		shuffled = sorted(smooshed_data, key=lambda k: np.random.random())
		shuffled_map=np.nan*np.ones(shape=data.shape)
		for j in range(0,len(pix)):
			m=pix[j]
			shuffled_map[m[0],m[1]]=shuffled[j]
		'''
		d = spatial.distance.cdist(pix,pix,'euclidean')
		np.savetxt(directory+'2018_d.txt',d)
		break

		masked_RM_errs=np.where(np.isfinite(data),RMerrs,np.nan)

		print("got here")


		for i in range(0,max_samples):
			print(i)
			dtheta_min=(i)*samplepix
			dtheta_max=(i+1)*samplepix
			SF[0,i]=3600*pixsize*(dtheta_min+dtheta_max)/(2.)
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
				pointA_err=masked_RM_errs[m[0],m[1]]
				#pointA_rand=shuffled_map[m[0],m[1]]
				q=d[:,j]
				roi=np.where((q>=dtheta_min),1,0)
				roi=np.where((q<=dtheta_max),roi,0)
				n_list=pix[np.where(roi==1)]

				for n in n_list:
					numpoints+=1
					pointB=data[n[0],n[1]]
					pointB_err=masked_RM_errs[n[0],n[1]]
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
				SF[2,i]=np.nan
				SF[3,i]=0
			else:
				values=np.asarray(values)
				value_errs=np.asarray(value_errs)
				biases=np.asarray(biases)

				SFval=np.nanmean(values)-np.nanmean(biases)
				err_SFval=np.sqrt(np.nanmean(np.square(value_errs)))

				SF[1,i]=SFval
				SF[2,i]=err_SFval
				SF[3,i]=numpoints
				print(np.sqrt(numpoints))

		filename=directory+src+'_RM_SF_witherrs.txt'
		np.savetxt(filename,SF,header='#scale (arcsec), SVval, err, npoints')

		markersizes=50*SF[3,:]/np.nanmax(SF[3,:])
		plt.figure()
		plt.scatter(SF[0,:],SF[1,:],s=markersizes)
		plt.errorbar(SF[0,:],SF[1,:],SF[2,:])
		plt.xlabel('Scale [arcsec]')
		plt.ylabel('SF')
		plt.yscale('log')
		plt.xlim(0,np.nanmax(SF[0,:]))
		plt.savefig(directory+src+'_SF_plot.png',dpi=300)

	#except:
		#print("something went wrong with this source")


