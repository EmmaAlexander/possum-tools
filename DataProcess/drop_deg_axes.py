#CASA script to create cutouts of fits cubes

directoryA = '/Volumes/TARDIS/Work/askap/'
directoryB = '/Volumes/NARNIA/pilot_cutouts/'

import numpy as np
sources=np.loadtxt('/Users/emma/GitHub/possum-tools/DataProcess/Oph_sources.txt',dtype='str')

for i in range(0,sources.shape[0]):
#if 1+1==2:
	src=sources[i,0]
	print(src)

	maxPI=directoryA+src+'/FDF_maxPI.fits'
	peakRM=directoryA+src+'/FDF_peakRM.fits'
	imsubimage(imagename=maxPI,outfile='maxPI_temp',dropdeg=True,overwrite=True)
	imsubimage(imagename=peakRM,outfile='peakRM_temp',dropdeg=True,overwrite=True)
	exportfits(imagename='maxPI_temp',fitsimage='FDF_maxPI_dropdeg.fits',overwrite=True)
	exportfits(imagename='peakRM_temp',fitsimage='FDF_peakRM_dropdeg.fits',overwrite=True)

	if sources[i,4]=='10083' or sources[i,4]=='10635':
		maxPIEMU=directoryB+src+'/EMU_FDF_maxPI.fits'
		peakRMEMU=directoryB+src+'/EMU_FDF_peakRM.fits'
		imsubimage(imagename=maxPIEMU,outfile='EMU_maxPI_temp',dropdeg=True,overwrite=True)
		imsubimage(imagename=peakRMEMU,outfile='EMU_peakRM_temp',dropdeg=True,overwrite=True)
		exportfits(imagename='EMU_maxPI_temp',fitsimage='EMU_FDF_maxPI_dropdeg.fits',overwrite=True)
		exportfits(imagename='EMU_peakRM_temp',fitsimage='EMU_FDF_peakRM_dropdeg.fits',overwrite=True)

	os.system("rm -r *_temp")
	os.system("mv *FDF_maxPI_dropdeg.fits {}/".format(src))
	os.system("mv *FDF_peakRM_dropdeg.fits {}/".format(src))

