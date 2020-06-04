"""
Code interpolates a FITS image into the grid defined by another FITS header, using ndimage.map_coordinates from scipy.
"""
import sys
import numpy as np
from scipy.ndimage import map_coordinates
from astropy.io import fits
from astropy.wcs import WCS

def loadfits(namefile):

    hdus=fits.open(namefile)
    if (isinstance(hdus,list)):
        hdu=hdus[0]
    else:
        hdu=hdus
    datacube = hdu.data
    hdr = hdu.header
    
    return datacube, hdr


def gridding(arg1, imagefile_2,fileout=False,fullWCS=True,ReturnHDU=False):
    """
    Interpolates Using ndimage and astropy.wcs for coordinate system.
    arg1 is the input data to be gridded, can be either a fits filename or an HDU
    arg2 contains the reference hdr, can either be a fits filename or a header
    """

    if (isinstance(arg1,str)):
        im1, hdr1 = loadfits(arg1)
    elif (isinstance(arg1,fits.hdu.image.PrimaryHDU)):
        im1 = arg1.data
        hdr1 = arg1.header
    elif (isinstance(arg1,fits.hdu.hdulist.HDUList)):
        im1 = arg1[0].data
        hdr1 = arg1[0].header
    else:
        sys.exit("not an recognized input format")
        
    
    if (isinstance(imagefile_2,str)):
        im2, hdr2 = loadfits(imagefile_2)
    else:
        hdr2=imagefile_2

    SetCRVAL3_1=False
    
    if ('CRVAL3' in hdr1.keys()):
        SetCRVAL3_1=True
        hdr1_CRVAL3=hdr1['CRVAL3']
    SetCRVAL3_2=False
    if ('CRVAL3' in hdr2.keys()):
        SetCRVAL3_2=True
        hdr2_CRVAL3=hdr2['CRVAL3']
        
    hdr1.pop('CRVAL3', None)  
    hdr2.pop('CRVAL3', None)  

        
    w1 = WCS(hdr1)
    w2 = WCS(hdr2)
    
    n2x = hdr2['NAXIS1']
    n2y = hdr2['NAXIS2']
    k2s=np.arange(0,n2x)
    l2s=np.arange(0,n2y)
    kk2s, ll2s = np.meshgrid(k2s, l2s)

    if (fullWCS):
        xxs2wcs, yys2wcs = w2.all_pix2world(kk2s, ll2s, 0)
        kk1s, ll1s = w1.all_world2pix(xxs2wcs,yys2wcs,0,tolerance=1e-12)
    else:
        xxs2wcs, yys2wcs = w2.wcs_pix2world(kk2s, ll2s, 0)
        kk1s, ll1s = w1.wcs_world2pix(xxs2wcs,yys2wcs,0)


    im1=np.nan_to_num(im1)
  
    resamp = map_coordinates(im1, [ll1s, kk1s],prefilter=False,order=1) #,order=1

    resamp=np.nan_to_num(resamp)

    if (fileout):
        fits.writeto(fileout,resamp, hdr2, overwrite=True)

    if (SetCRVAL3_1):
        hdr1['CRVAL3']=hdr1_CRVAL3
    if (SetCRVAL3_2):
        hdr2['CRVAL3']=hdr2_CRVAL3


    if ReturnHDU:
        rethdu=fits.PrimaryHDU()
        rethdu.data=resamp
        rethdu.header=hdr2
        return rethdu
    else:
        return resamp


    


