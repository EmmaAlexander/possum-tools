from functions import fit, peakfit, quadfit, image_rms, cube_rms, fitsopen,fdf_fit
import argparse
from astropy.io import fits
from astropy.wcs import WCS
import numpy as np 
import warnings
warnings.filterwarnings('ignore')
import os
c=299792458 #m/s

rootdirectory='/Volumes/TARDIS/Work/askap/'
defaultfile=rootdirectory+'final_sourcelist.txt'
sig_crop_pix=20

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
		ithresh=3*possum_rms

		#Read in data
		print("Reading in the fits files")

		#continuum image 
		cont_im,header=fitsopen(directory+src+'_cont.fits')
		wcs= WCS(header)

		#load the cubes 
		i_cube_data,cube_head=fitsopen(directory+src+'_i.fits')
		q_cube_data,qcube_head=fitsopen(directory+src+'_q.fits')
		u_cube_data,ucube_head=fitsopen(directory+src+'_u.fits')

		#frequencies
		ref_freq=cube_head['CRVAL3']
		ref_freq_pos=cube_head['CRPIX3']
		delta_freq=cube_head['CDELT3']
		nfreqchans=cube_head['NAXIS3']
		freqs=np.divide(np.linspace(ref_freq,(ref_freq+(delta_freq*nfreqchans)),
			nfreqchans),1000000000)
		lamda=np.divide(c,(freqs*1000000000))
		lamdasq=np.square(lamda)

		I_cube_rms=cube_rms(i_cube_data,sig_crop_pix)
		U_cube_rms=cube_rms(u_cube_data,sig_crop_pix)
		Q_cube_rms=cube_rms(q_cube_data,sig_crop_pix)

		QUrms=np.nanmean([U_cube_rms,Q_cube_rms],axis=0)
		QU_rms_arr=np.multiply(QUrms[:,np.newaxis,np.newaxis],np.ones(
			shape=u_cube_data.shape))

		p_cube_uncorr=np.sqrt(np.add(np.square(u_cube_data),np.square(
			q_cube_data)))
		p_cube=np.sqrt(np.square(p_cube_uncorr)-np.square(QU_rms_arr))

		polfrac_cube=np.divide(p_cube,i_cube_data)
		polfrac_map=np.nanmean(polfrac_cube,axis=0)

		pa_cube=0.5*np.arctan2(u_cube_data,q_cube_data)
		pa_cube_errs=0.5*(np.divide(QU_rms_arr,p_cube))

		average_p_map=np.nanmean(p_cube,axis=0)
		average_p_map_rms=image_rms(average_p_map,sig_crop_pix)


		for i in range (0,average_p_map.shape[1]):
			for j in range (0,average_p_map.shape[0]):
				if cont_im[j,i]>=contthresh:
					slice_pa_errs=pa_cube_errs[:,j,i]
					slice_pa_errs=np.nan_to_num(slice_pa_errs,nan=np.nanmax(
						slice_pa_errs))
            
					unwrapped = 2*slice_pa
					unwrapped[~np.isnan(unwrapped)] = 0.5*np.unwrap(
						unwrapped[~np.isnan(unwrapped)],discont=np.pi/2.)
            
					p_rad,pError_rad,redchisqrd_rad=fit(lamdasq,unwrapped,
						slice_pa_errs)
					RM=p_rad[0]
					RM_err=pError_rad[0]
					chi_0=p_rad[1]
					chi_0_err=pError_rad[1]
            
					fitted_RMs[j,i]=RM
					fitted_RM_errs[j,i]=RM_err
					redchisqrds[j,i]=redchisqrd_rad
					chi_0s[j,i]=chi_0
					chi_0_errs[j,i]=chi_0_err


###############################################################################
if __name__ == "__main__":

	ap = argparse.ArgumentParser()
	ap.add_argument("-f","--filename",help='Name of file with source list',
		default=defaultfile)
	args = ap.parse_args()

	main(args)