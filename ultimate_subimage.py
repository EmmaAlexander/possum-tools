#trying to tie all my bloody subimage scripts together

directory = '/Volumes/TARDIS/Work/askap/'
altdirectory = '/Volumes/NARNIA/'

import numpy as np
sources=np.loadtxt('subimage_sources.txt',dtype='str')

for i in range(0,sources.shape[0]):
	objectname=sources[i,0]
	POSSUMSB=sources[i,1]
	EMUSB=sources[i,2]
	sourcecentre=sources[i,3]
	print(objectname)

	region='centerbox[['+sourcecentre+'], [10arcmin, 10arcmin]]'
	possum_outfile=directory+objectname+'/'+objectname+'_POSSUM.fits'
	emu_outfile=directory+objectname+'/'+objectname+'_EMU.fits'

	#POSSUM
	if POSSUMSB == '5038':
		#this is the Early Science data 
		possum_cont_filename = '/Volumes/NARNIA/PawseySync/DRAGN_1_0p8_A/DRAGN_1_0p8_A/image.i.SB5038.cont.restored.fits'

	else:
		possum_cont_filename = directory +'fullfields/image.i.SB'+POSSUMSB+'.cont.taylor.0.restored.fits'

	if POSSUMSB == '10035':
		print('Skipping POSSUM: bad SB10035')
	else:
		imsubimage(imagename=possum_cont_filename,outfile='possum_cont_temp',region=region,overwrite=True,dropdeg=True)
		exportfits(imagename='possum_cont_temp',fitsimage=possum_outfile,overwrite=True)

		#cubes 
		i_filename = directory + 'fullfields/image.restored.i.SB'+POSSUMSB+'.contcube.total.fits'
		q_filename = directory + 'fullfields/image.restored.q.SB'+POSSUMSB+'.contcube.total.fits'
		u_filename = directory + 'fullfields/image.restored.u.SB'+POSSUMSB+'.contcube.total.fits'

		imsubimage(imagename=i_filename,outfile='i_im_temp',region=region,overwrite=True,dropdeg=True)
		imsubimage(imagename=q_filename,outfile='q_im_temp',region=region,overwrite=True,dropdeg=True)
		imsubimage(imagename=u_filename,outfile='u_im_temp',region=region,overwrite=True,dropdeg=True)
		exportfits(imagename='i_im_temp',fitsimage=objectname+'_POSSUM_i.fits',overwrite=True)
		exportfits(imagename='q_im_temp',fitsimage=objectname+'_POSSUM_q.fits',overwrite=True)
		exportfits(imagename='u_im_temp',fitsimage=objectname+'_POSSUM_u.fits',overwrite=True)
	#EMU
	if EMUSB != 'NaN':

		if EMUSB=='10083':
			i_EMU_filename = '/Volumes/NARNIA/fullfields/image.restored.i.SB10083.contcube.conv.fits'
			q_EMU_filename = '/Volumes/NARNIA/fullfields/image.restored.q.SB10083.contcube.conv.fits'
			u_EMU_filename = '/Volumes/NARNIA/fullfields/image.restored.u.SB10083.contcube.conv.fits'
			cont_EMU_filename= '/Volumes/NARNIA/fullfields/image.i.SB10083.cont.taylor.0.restored.conv.fits'

			imsubimage(imagename=i_EMU_filename,outfile='i_EMU_im_temp',region=region,overwrite=True,dropdeg=True)
			imsubimage(imagename=q_EMU_filename,outfile='q_EMU_im_temp',region=region,overwrite=True,dropdeg=True)
			imsubimage(imagename=u_EMU_filename,outfile='u_EMU_im_temp',region=region,overwrite=True,dropdeg=True)
			imsubimage(imagename=cont_EMU_filename,outfile='EMU_cont_im_temp',region=region,overwrite=True,dropdeg=True)
			exportfits(imagename='i_EMU_im_temp',fitsimage=objectname+'_EMU_i.fits',overwrite=True)
			exportfits(imagename='q_EMU_im_temp',fitsimage=objectname+'_EMU_q.fits',overwrite=True)
			exportfits(imagename='u_EMU_im_temp',fitsimage=objectname+'_EMU_u.fits',overwrite=True)
			exportfits(imagename='EMU_cont_im_temp',fitsimage=emu_outfile,overwrite=True)
		elif EMUSB=='10635':
			i_EMU_filename = '/Volumes/NARNIA/fullfields/image.restored.i.SB10635.contcube.v2.conv.fits'
			q_EMU_filename = '/Volumes/NARNIA/fullfields/image.restored.q.SB10635.contcube.v2.conv.fits'
			u_EMU_filename = '/Volumes/NARNIA/fullfields/image.restored.u.SB10635.contcube.v2.conv.fits'
			cont_EMU_filename= '/Volumes/NARNIA/fullfields/image.i.SB10635.cont.taylor.0.restored.fits'
			imsubimage(imagename=i_EMU_filename,outfile='i_EMU_im_temp',region=region,overwrite=True,dropdeg=True)
			imsubimage(imagename=q_EMU_filename,outfile='q_EMU_im_temp',region=region,overwrite=True,dropdeg=True)
			imsubimage(imagename=u_EMU_filename,outfile='u_EMU_im_temp',region=region,overwrite=True,dropdeg=True)
			imsubimage(imagename=cont_EMU_filename,outfile='EMU_cont_im_temp',region=region,overwrite=True,dropdeg=True)
			exportfits(imagename='i_EMU_im_temp',fitsimage=objectname+'_EMU_i.fits',overwrite=True)
			exportfits(imagename='q_EMU_im_temp',fitsimage=objectname+'_EMU_q.fits',overwrite=True)
			exportfits(imagename='u_EMU_im_temp',fitsimage=objectname+'_EMU_u.fits',overwrite=True)
			exportfits(imagename='EMU_cont_im_temp',fitsimage=emu_outfile,overwrite=True)
		else:
			#no cubes 
			emu_filename= directory +'fullfields/image.i.SB'+EMUSB+'.cont.taylor.0.restored.fits'
			imsubimage(imagename=emu_filename,outfile='emu_cont_temp',region=region,overwrite=True,dropdeg=True)
			exportfits(imagename='emu_cont_temp',fitsimage=emu_outfile,overwrite=True)
			os.system("rm -r emu_cont_temp")

	#tidy up
	os.system("rm -r *_temp")
	os.system("mv *{}* {}/".format(objectname,objectname))
