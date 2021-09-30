import numpy as np 
import math as mth
import matplotlib.pyplot as plt
from mpl_toolkits import axes_grid1
from astropy.io import fits
import math as m
import sys
from scipy.optimize import curve_fit


def fit(x,y,yerr):
    #fit a straight line 
    w=np.divide(1,yerr*yerr)
    df= len(x)-2 #degrees of freedom
    #filter out NaNs by indexing
    idx = np.isfinite(x) & np.isfinite(y) & np.isfinite(w)
    if len(x[idx]) >=2:
        try:
            p,cov=np.polyfit(x[idx],y[idx],1,cov=True,w=w[idx])
            # calculate chisquared
            chisqrd = 0
            for i, j, k in zip(x[idx], y[idx], yerr[idx]):
                c = pow(((j - np.polyval(p, i))/k), 2)
                chisqrd += c
            if df !=0:
                redchisqrd = chisqrd/df
            #cov = fit[1] * (len(x) - 2 - 2)/chisqrd
            pError = np.sqrt(np.diag(cov))
        except: 
            print("Something went wrong with a fit")
            p=[np.nan,np.nan]
            pError=[np.nan,np.nan]
            redchisqrd=np.nan
    else:
        p=[np.nan,np.nan]
        pError=[np.nan,np.nan]
        redchisqrd=np.nan
    return p,pError,redchisqrd

def peakfit(x,y):
    #find location of peak datapoint
    peak=np.nanmax(y)
    peakloc=np.where(y==peak)[0]

    if len(peakloc)>=2 or peakloc==0:
        #something is wrong
        fitted_peak=np.nan
        fitted_peak_loc=np.nan   
    else: 
    
        if len(peakloc)==2:
            #two adjacent points with same value, take the first
            peakloc=peakloc[0]

        #get a selection of points either side of this (5 in total)
        x_fit=x[int(peakloc-1):int(peakloc+2)]
        y_fit=y[int(peakloc-1):int(peakloc+2)]
        
        #do a quadratic fit
        peakfit=np.polyfit(x_fit,y_fit,2)
        #y = ax^2 + bx +c
        #dy/dx = 2ax+b = 0 at peak
        a=peakfit[0]
        b=peakfit[1]
        c=peakfit[2]
        fitted_peak_loc=-0.5*(b/a)
        fitted_peak=a*np.square(fitted_peak_loc) + b*fitted_peak_loc + c

    return fitted_peak_loc,fitted_peak

def gauss(x, H, A, x0, sigma):
    return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

def gaussfit(x,y):
    peak=np.nanmax(y)
    peakloc=np.where(y==peak)[0]

    sigma=50. #just for a placeholder width

    if len(peakloc)>=2 or peakloc==0:
        #something is wrong
        fitted_peak=np.nan
        fitted_peak_loc=np.nan   
    else: 
    
        if len(peakloc)==2:
            #two adjacent points with same value, take the first
            peakloc=peakloc[0]

        mean=float(x[peakloc])

        #do a guassian fit python
        popt, pcov = curve_fit(gauss, np.array(x,dtype='float64'), np.array(y,dtype='float64'), p0=[np.nanmin(y), np.nanmax(y), mean, sigma])
        return popt



def quadfit(x,y,xloc):
    #Do a quadratic fit to data
    idx = np.nanmedian(np.searchsorted(x, xloc, side="left"))
    x_fit=x[int(idx-1):int(idx+2)]
    y_fit=y[int(idx-1):int(idx+2)]
    try:
        p=np.polyfit(x_fit,y_fit,2)
    except:
        print("Something went wrong with a fit")
        p=[np.nan,np.nan,np.nan]
    return p

def image_rms(a,rms_crop_pix):
    '''Take the rms around the edge of the image then use it to do a second pass'''
    rms1=np.sqrt(np.nanmean(np.square((a[0:rms_crop_pix,0:-rms_crop_pix]))))
    rms2=np.sqrt(np.nanmean(np.square((a[rms_crop_pix:,0:rms_crop_pix]))))
    rms3=np.sqrt(np.nanmean(np.square((a[:-rms_crop_pix,rms_crop_pix:]))))
    rms4=np.sqrt(np.nanmean(np.square((a[0:-rms_crop_pix,-rms_crop_pix:]))))
    rms=np.median([rms1,rms2,rms3,rms4]) 
    return rms

def cube_rms(a,crop_pix):
    rms1=np.sqrt(np.nanmean(np.square(a[:,0:crop_pix,0:-crop_pix]),axis=(1,2)))
    rms2=np.sqrt(np.nanmean(np.square(a[:,crop_pix:,0:crop_pix]),axis=(1,2)))
    rms3=np.sqrt(np.nanmean(np.square(a[:,:-crop_pix,crop_pix:]),axis=(1,2)))
    rms4=np.sqrt(np.nanmean(np.square(a[:,0:-crop_pix,-crop_pix:]),axis=(1,2)))
    rms=np.nanmean([rms1,rms2,rms3,rms4],axis=0)
    return rms    

def add_colorbar(im, aspect=40, pad_fraction=0.5, **kwargs):
    """Add a vertical color bar to an image plot."""
    divider = axes_grid1.make_axes_locatable(im.axes)
    width = axes_grid1.axes_size.AxesY(im.axes, aspect=1./aspect)
    pad = axes_grid1.axes_size.Fraction(pad_fraction, width)
    current_ax = plt.gca()
    cax = divider.append_axes("right", size=width, pad=pad)
    plt.sca(current_ax)
    return im.axes.figure.colorbar(im, cax=cax, **kwargs)


def fdf_fit(pkrm_im,fdfdata,fdf_real,fdf_im,rmarray,cont_im,rmthresh=1000,contthresh=0):
    #fit the peak RM of a FDF
    rmsf_rm_fit=np.nan*np.ones(shape=pkrm_im.shape)
    rmsf_chi0_fit=np.nan*np.ones(shape=pkrm_im.shape)
    shape=pkrm_im.shape[1]
    print(shape)
    for i in range (0,shape):
        print(i,shape)
        for j in range (0,pkrm_im.shape[0]):
            if np.abs(pkrm_im[j,i])<=rmthresh and cont_im[j,i]>=contthresh:
                total_peak,total_peak_amp=peakfit(rmarray,fdfdata[:,j,i])
                p_real=quadfit(rmarray,fdf_real[:,j,i],total_peak)
                p_im=quadfit(rmarray,fdf_im[:,j,i],total_peak)
                q_amp=np.multiply(np.square(total_peak),p_real[0]) + np.multiply(total_peak,p_real[1]) + p_real[2]
                u_amp=np.multiply(np.square(total_peak),p_im[0]) + np.multiply(total_peak,p_im[1]) + p_im[2]
                chi_0_rmsf=0.5*np.arctan2(u_amp,q_amp)
                if np.isfinite(total_peak):
                    rmsf_rm_fit[j,i]=total_peak
                    rmsf_chi0_fit[j,i]=chi_0_rmsf

    return(rmsf_rm_fit,rmsf_chi0_fit)

def fdf_fit_gauss(pkrm_im,fdfdata,fdf_real,fdf_im,rmarray,cont_im,rmthresh=1000,contthresh=0):
    #fit the peak RM of a FDF
    rmsf_rm_fit_width=np.nan*np.ones(shape=pkrm_im.shape)
    shape=pkrm_im.shape[1]
    print(shape)
    for i in range (0,shape):
        for j in range (0,pkrm_im.shape[0]):
            if np.abs(pkrm_im[j,i])<=rmthresh and cont_im[j,i]>=contthresh:
                try:
                    H, A, x0, sigma=gaussfit(rmarray,fdfdata[:,j,i])
                except:
                    print(x0)
                FWHM=2.35482 * sigma

                if np.isfinite(x0):
                    rmsf_rm_fit_width[j,i]=FWHM

    return(rmsf_rm_fit_width)

def fitsopen(filepath):
    #open a fits file and return the data and header
    hdu=fits.open(filepath)
    header=hdu[0].header
    data=hdu[0].data
    hdu.close()
    return(data,header)

def finitelims(array):
    #finds the limits of actual data in an images
    x_collapse=np.nanmean(array,axis=0)
    y_collapse=np.nanmean(array,axis=1)
    x_coords=np.where(np.isfinite(x_collapse))
    xmin=np.nanmin(x_coords)
    xmax=np.nanmax(x_coords)
    y_coords=np.where(np.isfinite(y_collapse))
    ymin=np.nanmin(y_coords)
    ymax=np.nanmax(y_coords)

    return(xmin,xmax,ymin,ymax)

def vmax(array,percent):
    vmax_val=np.percentile(array[np.isfinite(array)],percent)
    return vmax_val

def get_rm_scale(rmarray,percent):
    rm_abs=np.abs(rmarray[np.isfinite(rmarray)])
    rm_abs=rm_abs[rm_abs<=1000]
    rm_scale_val=np.percentile(rm_abs,percent)
    print(rm_scale_val)
    return rm_scale_val




