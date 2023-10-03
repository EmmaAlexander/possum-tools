import matplotlib.pyplot as plt
import matplotlib.colors as colors
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs import utils
import numpy as np 
#from astroquery.skyview import SkyView
import astropy.units as u
from astropy.coordinates import SkyCoord
import glob, os, sys
from astropy.nddata import Cutout2D
import math
import cmasher as cmr
from astropy.visualization import lupton_rgb
#from reproject.mosaicking import find_optimal_celestial_wcs
#from reproject import reproject_interp
#from reproject.mosaicking import reproject_and_coadd
from matplotlib.patches import Ellipse
import warnings
from spectral_cube import SpectralCube
import regions 
from astropy.visualization import lupton_rgb
from reproject.mosaicking import find_optimal_celestial_wcs
from reproject import reproject_interp
from reproject.mosaicking import reproject_and_coadd
from astroquery.ipac.ned import Ned
import pandas as pd

dataloc='/share/nas2/ela/ASKAP/data/'
SB='10083'
cat=dataloc+'AS103_Continuum_Island_Catalogue_10083_1641.csv'

imagefile=dataloc+'image.i.SB'+SB+'.cont.taylor.0.restored.fits'
ifile=dataloc+'image.restored.i.SB'+SB+'.contcube.fits'
ufile=dataloc+'image.restored.u.SB'+SB+'.contcube.fits'
qfile=dataloc+'image.restored.q.SB'+SB+'.contcube.fits'


sbdir='/share/nas2/ela/ASKAP/data/'+SB
if os.path.isdir(sbdir) == False:
    os.system('mkdir '+ sbdir)

contourexps=np.arange(start=0,stop=32,step=1)
contourmults=np.power(2,contourexps)

DEStiles=np.genfromtxt('DES_tiles_dims.txt',dtype='str')

warnings.filterwarnings('ignore')
hdu=fits.open(imagefile)
image=hdu[0].data.squeeze()
header=hdu[0].header
wcs= WCS(hdu[0].header).celestial

cont=SpectralCube.read(hdu)
hdu.close()

pixscale=header['CDELT2']
bmaj_pix=header['BMAJ']/header['CDELT2']
bmin_pix=header['BMIN']/header['CDELT2']
bpa=header['BPA']

ihdu=fits.open(ifile)
#icube=ihdu[0].data.squeeze()
cubeheader=ihdu[0].header
cubewcs= WCS(ihdu[0].header).celestial
icube=SpectralCube.read(ihdu)
ihdu.close()

qhdu=fits.open(qfile)
#qcube=qhdu[0].data.squeeze()
qcube=SpectralCube.read(qhdu)
qhdu.close()

uhdu=fits.open(ufile)
#ucube=uhdu[0].data.squeeze()
ucube=SpectralCube.read(uhdu)
uhdu.close()

df=pd.read_csv(cat)

min_size=0.5
extended=df.loc[(df['solid_angle'] >= min_size)]
extended.to_csv(sbdir+'/'+SB+'_extended.csv')

bands=['g','i','r']
des_data_loc=dataloc+'DES/'
resume=True
for ind in extended.index:
    name=str(extended['island_name'][ind])
    #if name=='J220045-492301':
        #resume=True
    if resume==True:
        print(name)
        obdir=sbdir+'/'+name
        if os.path.isdir(obdir) == False:
            os.system('mkdir '+ obdir)
        ra=extended['ra_deg_cont'][ind]
        dec=extended['dec_deg_cont'][ind]
        coords=SkyCoord(ra,dec,unit=u.deg)
    
        background_noise=extended['background_noise'][ind] 
        flux_peak=extended['flux_peak'][ind]/1000 
        
        basecont=3.*background_noise/1000.
        radio_contours = [basecont * i for i in contourmults]
    
        #cutoutsize=10*u.arcmin #arcmins
        #radio_cutout = Cutout2D(image, wcs=wcs, mode='trim',position=coords, size=cutoutsize)
    
        x_min=int(extended['x_min'][ind])-10
        x_max=int(extended['x_max'][ind])+10
        y_min=int(extended['y_min'][ind])-10
        y_max=int(extended['y_max'][ind])+10
    
        c1=utils.pixel_to_skycoord(x_min,y_min,wcs)
        c2=utils.pixel_to_skycoord(x_max,y_max,wcs)
    
        cubemin=utils.skycoord_to_pixel(c1,cubewcs)
        cubemax=utils.skycoord_to_pixel(c2,cubewcs)
    
        ra_min=c2.ra.deg
        ra_max=c1.ra.deg
        dec_min=c1.dec.deg
        dec_max=c2.dec.deg
    
        #maximum dimensions of cutout in arcmins
        arcmins=np.nanmax([ra_max-ra_min,dec_max-dec_min])*60
    
        #results table from NED
        result_table = Ned.query_region(coords, radius=arcmins*u.arcmin, equinox='J2000')
        result_table.write(obdir+'/'+name+'_ned.txt', format='ascii', overwrite=True)  
    
        contsubcube = cont[:,y_min:y_max, x_min:x_max]
        contsubcube.write(obdir+"/"+name+"_cont.fits",overwrite=True)
        #alternative below
        #fits.PrimaryHDU(contsubcube,header=contsubcube.header).writeto(obdir+"/"+name+"_cont.fits",overwrite=True)
    
        isubcube=icube[:, int(cubemin[1]):int(cubemax[1]), int(cubemin[0]):int(cubemax[0])]
        isubcube.write(obdir+"/"+name+"_icube.fits",overwrite=True)
    
        qsubcube=qcube[:, int(cubemin[1]):int(cubemax[1]), int(cubemin[0]):int(cubemax[0])]
        qsubcube.write(obdir+"/"+name+"_qcube.fits",overwrite=True)
    
        usubcube=ucube[:, int(cubemin[1]):int(cubemax[1]), int(cubemin[0]):int(cubemax[0])]
        usubcube.write(obdir+"/"+name+"_ucube.fits",overwrite=True)
        
        plt.figure(dpi=300)
        #plt.imshow(radio_cutout.data,vmin=basecont/3.,vmax=flux_peak,origin='lower')
        #plt.contour(radio_cutout.data,levels=radio_contours,colors='white',linewidths=0.2)       
        
        plt.imshow(contsubcube[0,:,:].data,vmin=basecont/3.,vmax=flux_peak,origin='lower')
        plt.contour(contsubcube[0,:,:].data,levels=radio_contours,colors='white',linewidths=0.2)       
        
        plt.axis('off')
        plt.annotate(name,(0.01,0.96),color='white',xycoords='axes fraction')
        plt.annotate('SB '+SB,(0.01,0.92),color='white',xycoords='axes fraction')
        plt.savefig(obdir+'/'+name+'_radio.png',dpi=300,bbox_inches='tight')
    
        indices=np.where((DEStiles[:,2].astype(float)>=ra_min) & (DEStiles[:,1].astype(float)<=ra_max) & (DEStiles[:,3].astype(float)>=dec_min) & (DEStiles[:,4].astype(float)<=dec_max))
        DES_tiles_to_use=DEStiles[indices,0][0]
        DES_tiles_to_use_coords=DEStiles[indices,:]
    
        R_list=[]
        G_list=[]
        B_list=[]
    
        for j in range(0,len(DES_tiles_to_use)):
            tile=DES_tiles_to_use[j]
            if len(glob.glob(dataloc+'DES/'+tile+'*_i.fits*'))==0:
                print('Downloading DES tile'.format(tile))
                for b in bands:
                    if os.path.isfile(des_data_loc+'{}_r4575p01_{}.fits.fz'.format(tile,b))==False and os.path.isfile(des_data_loc+'{}_r4575p02_{}.fits.fz'.format(tile,b))==False:
                        #print('{}_r4575p0X_{}.fits.fz tileset {}/{}'.format(tile,b))
                        print('Downloading {}_r4575p0X_{}.fits.fz'.format(tile,b))
                        os.system('curl -f https://desdr-server.ncsa.illinois.edu/despublic/dr2_tiles/{}/{}_r4575p02_{}.fits.fz -L -o '.format(tile,tile,b)+des_data_loc+'{}_r4575p02_{}.fits.fz'.format(tile,b))
                        os.system('curl -f https://desdr-server.ncsa.illinois.edu/despublic/dr2_tiles/{}/{}_r4575p01_{}.fits.fz -L -o '.format(tile,tile,b)+des_data_loc+'{}_r4575p01_{}.fits.fz'.format(tile,b))
                        os.system('curl -f https://desdr-server.ncsa.illinois.edu/despublic/dr2_tiles/{}/{}_r5137p01_{}.fits.fz -L -o '.format(tile,tile,b)+des_data_loc+'{}_r4575p01_{}.fits.fz'.format(tile,b))
                        #if anything fails check for the dr1 tile
                        #os.system('curl -f https://desdr-server.ncsa.illinois.edu/despublic/dr1_tiles/{}/{}_r2597p01_{}.fits.fz -L -o '.format(tile,tile,b)+des_data_loc+'{}_r4575p01_{}.fits.fz'.format(tile,b))
                    else:
                        print("Tile already exists for band {}".format(b))
        
        for j in range(0,len(DES_tiles_to_use)):
            #for j in ind:
            Rhdu=fits.open(glob.glob(dataloc+'DES/'+DES_tiles_to_use[j]+'*_i.fits*')[0])
            R=Rhdu[1].data
            des_wcs=WCS(Rhdu[1].header)
            Rhdu.close()
    
            Ghdu=fits.open(glob.glob(dataloc+'DES/'+DES_tiles_to_use[j]+'*_r.fits*')[0])
            G=Ghdu[1].data
            Ghdu.close()
    
            Bhdu=fits.open(glob.glob(dataloc+'DES/'+DES_tiles_to_use[j]+'*_g.fits*')[0])
            B=Bhdu[1].data
            Bhdu.close()
    
            R_cutout=Cutout2D(R,position=coords,size=1.05*arcmins*u.arcmin,wcs=des_wcs,mode='trim')
            G_cutout=Cutout2D(G,position=coords,size=1.05*arcmins*u.arcmin,wcs=des_wcs,mode='trim')
            B_cutout=Cutout2D(B,position=coords,size=1.05*arcmins*u.arcmin,wcs=des_wcs,mode='trim')
    
            R_hdu=fits.PrimaryHDU(data=R_cutout.data, header=R_cutout.wcs.to_header())
            G_hdu=fits.PrimaryHDU(data=G_cutout.data, header=G_cutout.wcs.to_header())
            B_hdu=fits.PrimaryHDU(data=B_cutout.data, header=B_cutout.wcs.to_header())
    
            R_list.append(R_hdu)
            G_list.append(G_hdu)
            B_list.append(B_hdu)  
    
        if len(R_list)==0:
            print("something wrong")
        elif len(R_list)==1:
            #only one image so no need to mosaic
            R=R_list[0].data
            G=G_list[0].data
            B=B_list[0].data
            des_wcs=WCS(R_list[0].header)
        else:
            #need to combine them
            des_wcs, shape_out = find_optimal_celestial_wcs(R_list)
            R, footprint_R = reproject_and_coadd(R_list,des_wcs,shape_out=shape_out,reproject_function=reproject_interp)
            G, footprint_G = reproject_and_coadd(G_list,des_wcs,shape_out=shape_out,reproject_function=reproject_interp)
            B, footprint_B = reproject_and_coadd(B_list,des_wcs,shape_out=shape_out,reproject_function=reproject_interp)
            
        img=lupton_rgb.make_lupton_rgb(R,G,B,Q=10,stretch=50,minimum=1)
        print("plotting")
        
        plt.figure(dpi=300)
        ax=plt.subplot(111,projection=contsubcube.wcs.celestial)
        ax.imshow(img,origin='lower',transform=ax.get_transform(des_wcs))
        ax.contour(contsubcube[0,:,:].data,levels=radio_contours,colors='white',linewidths=0.2) 
        ax.set_xlim(0,contsubcube[0,:,:].data.shape[1])
        ax.set_ylim(0,contsubcube[0,:,:].data.shape[0])
        plt.axis('off')
        ax.annotate(name,(0.01,0.96),color='white',xycoords='axes fraction')
        ax.annotate('SB '+SB,(0.01,0.92),color='white',xycoords='axes fraction')
        for result in result_table:
            if result['Type']=='RadioS':
                ax.scatter(result['RA'],result['DEC'],transform=ax.get_transform('fk5'),marker='+',color='green',linewidths=0.5)
            else:
                if result['Redshift']!='--':
                    ax.scatter(result['RA'],result['DEC'],transform=ax.get_transform('fk5'),marker='o',edgecolor='green',facecolor='None')
                else:
                    ax.scatter(result['RA'],result['DEC'],transform=ax.get_transform('fk5'),marker='o',edgecolor='green',facecolor='None',linestyle=':')
        plt.savefig(obdir+'/'+name+'_optical.png',dpi=300,bbox_inches='tight')

