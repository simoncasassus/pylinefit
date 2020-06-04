import sys
import numpy as np
import astropy
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.widgets import RectangleSelector
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pprint import pprint as pp

include_path='/Users/simon/common/python/include/'
sys.path.append(include_path)
import Resamp

## WIDGET INBUILT KEY STROKES
# g: grid


# https://matplotlib.org/users/event_handling.html
# https://matplotlib.org/3.1.1/gallery/widgets/rectangle_selector.html

def pix2wcs_0CRVAL(j,i):
    hdr=hdu.header
    if ('pixel' in hdr['CTYPE1']):
        pixscale=1.
    else:
        pixscale=3600.*hdr['CDELT2']

    a=-pixscale*(i-(hdr['CRPIX1']-1.))
    d=pixscale*(j-(hdr['CRPIX2']-1.))
        
    return (d,a)

def wcs2pix_0CRVAL(d,a):
    hdr=hdu.header
    if ('pixel' in hdr['CTYPE1']):
        pixscale=1.
    else:
        pixscale=3600.*hdr['CDELT2']
        
    i = -a/pixscale+(hdr['CRPIX1']-1.)
    j=d/pixscale+(hdr['CRPIX2']-1.)
    return (j,i)

def line_select_callback(eclick, erelease):
    'eclick and erelease are the press and release events'
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata
    diagdistance=np.sqrt( (y2-y1)**2 + (x2-x1)**2)
    print("(%3.2f, %3.2f) --> (%3.2f, %3.2f): Diag %3.2f" % (x1, y1, x2, y2,diagdistance))
    #print(" The buttons you used were: %s %s" % (eclick.button, erelease.button))

    im=hdu.data
    (j1,i1)=wcs2pix_0CRVAL(y1,x1)
    (j2,i2)=wcs2pix_0CRVAL(y2,x2)
    iis=[i1,i2]
    (i1m,i2m) =map(lambda i: int(i), sorted(iis))
    jjs=[j1,j2]
    (j1m,j2m) =map(lambda j: int(j), sorted(jjs))

    #print("i1m ",i1m," j1m ",j1m,"i2m ",i2m," j2m ",j2m)
    subim=im[j1m:j2m,i1m:i2m]

    print("Flux: %.3e  Rms: %.3e  Min: %.3e Max: %.3e S/N: %.3e" % (np.sum(subim),np.std(subim),np.min(subim),np.max(subim),  np.max(subim)/np.std(subim)))

    range2=np.max(subim[np.where(np.isfinite(subim))])
    range1=np.min(subim[np.where(np.isfinite(subim))])

    theimage.set_clim(vmin=range1, vmax=range2)


    #zetransform=theimage.get_transform()
    #invzetransform=zetransform.inverted()
    ##print("zetransform ",zetransform)
    ##help(zetransform)
    #(i1,j1)=zetransform.transform_affine((x1,y1))
    #(i2,j2)=zetransform.transform_affine((x2,y2))
    #(j1,i1)=zetransform.transform_affine((y1,x1))
    #(j2,i2)=zetransform.transform_affine((y2,x2))
    #(invx1,invy1)=invzetransform.transform((i1,j1))
    #(invx2,invy2)=invzetransform.transform((i2,j2))

    if (eclick.dblclick):
        # on double click clear and deactivate RectangleSelector
        # had to place this here as a mouse event because it bugs in toggle_selector(event), where it gives a logx axis and an error.
        toggle_selector.RS.set_visible(False)
        toggle_selector.RS.set_active(False)

    fig1.canvas.draw()


def toggle_selector(event):
    #print(' Key pressed.',event.key)

    # return
    if event.key in ['C', 'c']:
        print("COLOR MAPS: Accent, Accent_r, Blues, Blues_r, BrBG, BrBG_r, BuGn, BuGn_r, BuPu, BuPu_r, CMRmap, CMRmap_r, Dark2, Dark2_r, GnBu, GnBu_r, Greens, Greens_r, Greys, Greys_r, OrRd, OrRd_r, Oranges, Oranges_r, PRGn, PRGn_r, Paired, Paired_r, Pastel1, Pastel1_r, Pastel2, Pastel2_r, PiYG, PiYG_r, PuBu, PuBuGn, PuBuGn_r, PuBu_r, PuOr, PuOr_r, PuRd, PuRd_r, Purples, Purples_r, RdBu, RdBu_r, RdGy, RdGy_r, RdPu, RdPu_r, RdYlBu, RdYlBu_r, RdYlGn, RdYlGn_r, Reds, Reds_r, Set1, Set1_r, Set2, Set2_r, Set3, Set3_r, Spectral, Spectral_r, Wistia, Wistia_r, YlGn, YlGnBu, YlGnBu_r, YlGn_r, YlOrBr, YlOrBr_r, YlOrRd, YlOrRd_r, afmhot, afmhot_r, autumn, autumn_r, binary, binary_r, bone, bone_r, brg, brg_r, bwr, bwr_r, cividis, cividis_r, cool, cool_r, coolwarm, coolwarm_r, copper, copper_r, cubehelix, cubehelix_r, flag, flag_r, gist_earth, gist_earth_r, gist_gray, gist_gray_r, gist_heat, gist_heat_r, gist_ncar, gist_ncar_r, gist_rainbow, gist_rainbow_r, gist_stern, gist_stern_r, gist_yarg, gist_yarg_r, gnuplot, gnuplot2, gnuplot2_r, gnuplot_r, gray, gray_r, hot, hot_r, hsv, hsv_r, inferno, inferno_r, jet, jet_r, magma, magma_r, nipy_spectral, nipy_spectral_r, ocean, ocean_r, pink, pink_r, plasma, plasma_r, prism, prism_r, rainbow, rainbow_r, seismic, seismic_r, spring, spring_r, summer, summer_r, tab10, tab10_r, tab20, tab20_r, tab20b, tab20b_r, tab20c, tab20c_r, terrain, terrain_r, twilight, twilight_r, twilight_shifted, twilight_shifted_r, viridis, viridis_r, winter, winter_r")
        cmap=input('Enter new color map> ')
        theimage.set_cmap(cmap)
        print('Switch to ',cmap)
        fig1.canvas.draw()
        return

    if event.key in ['A', 'a'] and not toggle_selector.RS.active:
        print('RectangleSelector activated.')
        toggle_selector.RS.set_active(True)
        toggle_selector.RS.set_visible(True)
        return

    if event.key in ['H', 'h']:
        print('key a: activate RectangleSelector (dbleclick deactivates)')
        print('key c: change colormap')
        return

def colorbar(Mappable, Orientation='vertical'):
    Ax = Mappable.axes
    fig = Ax.figure
    divider = make_axes_locatable(Ax)
    Cax = divider.append_axes("right", size="5%", pad=0.08)
    return fig.colorbar(
        mappable=Mappable,
        cax=Cax,
        use_gridspec=True,
        orientation=Orientation,
        format="%.1e"
    )

def View(indata,cmap='RdBu_r'):

    #help(indata)
    #indatatype=type(indata)

    global hdu

    AddContours=False
    if isinstance(indata,list): # astropy.io.fits.hdu.hdulist.HDUList):
        #print("This is an HDU List")
        hdu=indata[0]
        if (len(indata) > 1):
            hducontours=indata[1]           
            AddContours=True

    elif isinstance(indata,np.ndarray):
        hdu = fits.PrimaryHDU()
        hdu.data = indata
        hdr=hdu.header
        hdr['CDELT1']=1.
        hdr['CDELT2']=1.
        (nx,ny)=indata.shape
        hdr['CRPIX1']=int(nx/2.)
        hdr['CRPIX2']=int(ny/2.)
        hdr['CTYPE1']='pixel'
        hdr['CTYPE2']='pixel'
            
    else: #  isinstance(indata,astropy.io.fits.hdu.hdulist):
        #print("This is an HDU")
        hdu=indata

    im=hdu.data
    hdr=hdu.header

    if (not 'CTYPE1' in hdr.keys()):
        hdr['CTYPE1']='pixel'
        hdr['CTYPE2']='pixel'
        (nx,ny)=im.shape
        hdr['CRPIX1']=int(nx/2.)
        hdr['CRPIX2']=int(ny/2.)
        
    
    # mpl.use('TkAgg')


    global fig1
    global ax1
    fig1,ax1=plt.subplots()

    (d0,a0)=pix2wcs_0CRVAL(0.,0.)
    (d1,a1)=pix2wcs_0CRVAL(hdr['NAXIS2']-1,hdr['NAXIS1']-1)
    #(j0,i0)=wcs2pix_0CRVAL(d0,a0)
    #(j1,i1)=wcs2pix_0CRVAL(d1,a1)

    #print("a0 d0 a1 d1", a0,d0,a1,d1)
    #print("i0 j0 i1 j1", i0,j0,i1,j1)


    range2=np.max(im[np.where(np.isfinite(im))])
    range1=np.min(im[np.where(np.isfinite(im))])

    MedianvalRange=False
    if MedianvalRange:
        typicalvalue=np.median(im)
        medrms=np.sqrt(np.median( (im - typicalvalue)**2))
        print("typical value ",typicalvalue," rms ",rms,"medrms",medrms)
        range1=typicalvalue-3.*medrms
        range2=typicalvalue+3.*medrms


    global theimage

    theimage=ax1.imshow(im, origin='lower', cmap=cmap, #norm=norm,
                               extent=[a0,a1,d0,d1], vmin=range1, vmax=range2, interpolation='nearest') #'nearest'  'bicubic'

    colorbar(theimage)

    if AddContours:
        hducontours_matched=Resamp.gridding(hducontours,hdr,ReturnHDU=True)
                                            
        imcontours=hducontours_matched.data
        levels=np.array([0.2,0.4,0.6,0.8])*np.max(imcontours)
        ax1.contour(imcontours,origin='lower', cmap=cmap, extent=[a0,a1,d0,d1], levels=levels)
                     
    
    #help(theimage)
    global ClearAll
    ClearAll=False
    toggle_selector.RS = RectangleSelector(ax1, line_select_callback,
                                           drawtype='box', useblit=True,
                                           button=[1, 3],  # don't use middle button
                                           minspanx=5, minspany=5,
                                           spancoords='pixels', interactive=True)
    #help(toggle_selector.RS)
    #print("Pretty printing togle_selector")
    #pp(toggle_selector.RS)

    toggle_selector.RS.set_active(False)


    plt.connect('key_press_event', toggle_selector)


    plt.show()



def Spec(indata):

    if not isinstance(indata,list): #
        indata=[indata,]


    global fig1
    global ax1
    fig1,ax1=plt.subplots()
    
    for ispec,aspec in enumerate(indata):
        theplot=ax1.plot(aspec[:,0],aspec[:,1],label=str(ispec))

    plt.legend()
    plt.show()
        
