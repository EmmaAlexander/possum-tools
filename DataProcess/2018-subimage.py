src='2018-5540'

region='box[[181pix,79pix], [681pix,816pix]]'

directory='/Volumes/NARNIA/pilot_cutouts/leakage_corrected/2018-5540/'

imsubimage(imagename='FDF_peakRM_fitted_corrected.fits',outfile='pkrm_smol_temp',region=region,overwrite=True,dropdeg=True)
exportfits(imagename='pkrm_smol_temp',fitsimage='2018-5540_pkrm_forSF.fits',overwrite=True)

imsubimage(imagename='2018-5540_FDF_peakRM_dropdeg.fits',outfile='dd_smol_temp',region=region,overwrite=True,dropdeg=True)
exportfits(imagename='dd_smol_temp',fitsimage='2018-5540_pkrm_dropdeg_forSF.fits',overwrite=True)

imsubimage(imagename='2018-5540_RM_error.fits',outfile='err_smol_temp',region=region,overwrite=True,dropdeg=True)
exportfits(imagename='err_smol_temp',fitsimage='2018-5540_RM_error_forSF.fits',overwrite=True)


os.system("rm -r *_temp")