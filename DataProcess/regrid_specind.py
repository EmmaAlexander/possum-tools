#CASA script to create cutouts of fits cubes

directoryA = '/Volumes/TARDIS/Work/askap/leakage_corrected/'
#directoryB = '/Volumes/NARNIA/pilot_cutouts/'

import numpy as np
sources=np.loadtxt('/Users/emma/GitHub/possum-tools/DataProcess/pilot_sources.txt',dtype='str')

#for i in range(0,sources.shape[0]):
for i in range(1,sources.shape[0]):
	objectname=sources[i,0]
	POSSUMSB=sources[i,3]
	EMUSB=sources[i,4]
	ra=sources[i,1]
	dec=sources[i,2]
	sourcecentre=ra+','+dec
	fov=sources[i,6]#arcsec
	print(objectname)

	srcpath=directoryA+objectname+'/'

	possum_cont_filename = srcpath+objectname+'_POSSUM.fits'
	emu_cont_filename=srcpath+objectname+'_EMU.fits'

	print(emu_cont_filename)
	print(possum_cont_filename)

	if POSSUMSB == '10035':
		print('Skipping: bad SB10035')
	else:

		importfits(fitsimage=emu_cont_filename,imagename='EMU_temp')
		importfits(fitsimage=possum_cont_filename,imagename='POSSUM_temp')
		imregrid(imagename='EMU_temp',template='POSSUM_temp',output='EMU_regrid_temp')
		imsmooth(imagename='POSSUM_temp',major="13arcsec",minor="13arcsec",pa="0deg",targetres=True,outfile='POSSUM_smooth13_temp')
		imsmooth(imagename='EMU_regrid_temp',major="13arcsec",minor="13arcsec",pa="0deg",targetres=True,outfile='EMU_regrid_smooth13_temp')

		exportfits(imagename='POSSUM_smooth13_temp',fitsimage=srcpath+objectname+'_POSSUM_smooth13.fits')
		exportfits(imagename='EMU_regrid_smooth13_temp',fitsimage=srcpath+objectname+'_EMU_smooth13.fits')


	#tidy up
	os.system("rm -r *_temp")
	#os.system("mv *{}* {}/".format(objectname,objectname))
