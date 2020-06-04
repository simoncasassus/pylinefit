"""
Package for CO rotational line fitting using uniform slab models. 
S. Casassus & F. Alarcon
"""

import MolData	 # molecular data

from multiprocessing import Pool
import numpy as np
import scipy as sp
from scipy.integrate import quad
import astropy.io.fits as pf
import os
import math
import sys
from iminuit import Minuit
from astropy.convolution import Gaussian2DKernel, convolve_fft
import astropy.units as u
import astropy.constants as const
#from astropy.modeling.blackbody import blackbody_nu
#from astropy.modeling.models import BlackBody
from copy import deepcopy

from tqdm import tqdm
import re

include_path='/Users/simon/common/python/include/'
sys.path.append(include_path)
import Vtools


def Tbrightness(I_nu,nu):
    # input I_nu in erg/s/cm2/sr/Hz
    # input nu in Hz
    if (I_nu < 0.):
        Tb=1.
    else:
        h_P=const.h.cgs.value
        k_B=const.k_B.cgs.value
        c_light=const.c.cgs.value
        Tb=h_P*nu/(k_B*np.log( 1. + (2. * h_P * nu**3  / (c_light**2 * I_nu))))
        
    return Tb

def loadfitsdata(namefile):
    hdu=pf.open(namefile)
    datacube = hdu[0].data
    hdr = hdu[0].header

    if 'BMAJ' in hdr.keys():
        bmaj=hdr['BMAJ']
        bmin=hdr['BMIN']
        bpa=hdr['BPA']
    elif (len(hdu)>1):
        print("no beam info, look for extra HDU")
        beamhdr=hdu[1].header
        beamdata=hdu[1].data
        bmaj=beamdata[0][0]
        bmin=beamdata[0][1]
        bpa=beamdata[0][2]
        hdr['BMAJ']=bmaj/3600.
        hdr['BMIN']=bmin/3600.
        hdr['BPA']=bpa
        
    if (len(datacube.shape) > 3):
        print("len(datacube)",len(datacube.shape))
        datacube=datacube[0,:,:,:]

    return datacube, hdr




#def bbody(T,nu):
#    """
#    Blackbody flux for a given temperature and frequency erg / (cm2 Hz s sr) (cgs system)
#    """
#    return blackbody_nu(nu, T).cgs.value
#


def bbody(T,nu):

    bb =  ( (2. * h_P * nu**3 ) / c_light**2) / (np.exp( h_P * nu /  (k_B * T)) - 1.)
    
    return bb





def phi(Tk,nu,nu0,vturb, molecule_mass):
    """
    Returns the normalized line profile.
    Tk: Temperature.
    nu: Array of frecuencies to sample the line profile.
    nu0: Center of line emission.
    vturb: Turbulent velocity or dispersion velocity along the line of sight (cgs system). 
    molecule_mass: Molecular mass, in g.
    """ 
    sigma_nu = (nu0/c_light)*np.sqrt(k_B*Tk/molecule_mass + vturb**2 )
    phi0 = 1./(sigma_nu*np.sqrt(2*np.pi))
    gaussprofile = phi0 * np.exp(-((nu-nu0)**2.0)/(2.*(sigma_nu**2.0)))

    #print('phi0',np.max(gaussprofile),phi0,nu0,np.mean(nu),sigma_nu,molecule_mass)

    return gaussprofile




def Kappa_line(Tk,iiso):

    levelenergies=levelenergiess[iiso]
    B_21=B_21s[iiso]
    E_lo=E_los[iiso]
    restfreq=restfreqs[iiso]
    
    Zpart = Part(levelenergies, g_Js, Tk)
    B_12 = B_21 * g_Jup/g_Jlo
    frac_lowerlevel = g_Jlo*np.exp(-(E_lo/ (k_B * Tk)))/ Zpart


    kappa_L = (h_P * restfreq / (4. *np.pi))* frac_lowerlevel * B_12 * (1. - np.exp( - (h_P * restfreq / (k_B * Tk)))) / mH2
    return kappa_L


def intensity(nu, Tk, nu0, Sigma_g, vturb,iiso):
        
    kappa_L=Kappa_line(Tk,iiso)

    molecule_mass=molecule_masses[iiso]

    f_abund=f_abunds[iiso]

    phiprof=phi(Tk,nu,nu0,vturb,molecule_mass)
    
    tau_L = kappa_L * Sigma_g * f_CO * f_abund *  phiprof

    phiprof0=phi(Tk,nu0,nu0,vturb,molecule_mass)

    tau_nu0 = kappa_L *  Sigma_g * f_CO *  f_abund *  phiprof0


    Iemerge = bbody(Tk,nu)*(1.0-np.exp(-tau_L))   # return units in CGS


    return  Iemerge, tau_nu0, tau_L


#def intensity_continuum(nu, T, nu0, alpha, Sigma_g, vturb, Icont_0):
#    cont = Icont_0*np.exp(-tau(T,nu,nu0,N_CO,vturb,angle, f_abund,  molecule_mass, sigma))*(nu/nu0)**alpha
#    opt_depth = tau(T,nu,nu0,N_CO,vturb,angle, f_abund,  molecule_mass, sigma)
#    blackbody = bbody(T,nu)*(1.0-np.exp(-opt_depth)) #*scaling
#    tau_nu0 = tau(T,nu0,nu0,N_CO,vturb,angle, f_abund,  molecule_mass, sigma)
#    return  blackbody + cont , tau_nu0, opt_depth
#


def Part(levelenergies, g_Js, Tk):
    return np.sum(g_Js*np.exp(-levelenergies/(k_B*Tk)))



def intensity_err(nu, nu0, Tk,Sigma_g, vturb, datos, rms,iiso):
    """
    returns chi2 for model vs data
    """
    model ,tau0, taus = intensity(nu, Tk, nu0, Sigma_g, vturb,iiso)
    aux = (datos-model)**2
    chi2 = np.sum(aux)/rms**2
    return chi2

def master_chi2(nuss, v0, Temp, Sigma_g, vturb, datas, rmss):

    chi2=0.
    for iiso,adata in enumerate(datas):
        nus=nuss[iiso]
        rms=rmss[iiso]
        restfreq=restfreqs[iiso]
        nu0=restfreq-(v0/c_light)*restfreq
        chi2+=intensity_err(nus, nu0, Temp, Sigma_g, vturb, adata, rms,iiso)

    return chi2

        
def parspar(n):
    j = n[0]
    i = n[1]


    T_inits=[]
    vel_peaks=[]
    I_peaks=[]
    datas=[]
    for iiso,acubo in enumerate(cubos):
        
        data = acubo[:,j,i]
        datas.append(data)
        nus=nuss[iiso]
        
        datamax = data.max()
        nu0_init= nus[np.argmax(data)] # selected_velocities[signal_a==signal_a.max()]

        aT_init = Tbrightness(datamax,nu0_init)
        
        T_inits.append(aT_init)

        velocities=velocitiess[iiso]
        vel_peak = velocities[data==data.max()][0]
        vel_peaks.append(vel_peak)
        I_peaks.append(data.max())
        
    T_init=T_inits[0] # max(T_inits)
    T_limits = (0.5*T_init,1.5*T_init)


    vel_peak=vel_peaks[0]
    v0_init=vel_peak * 1E3 * 1E2 # init centroid velocity in CGS

    vturb_init=0.

    Sigma_g_thins=[]
    Sigma_g_tauones=[]
    rmss=[]
    for iiso,acubo in enumerate(cubos):
        data = acubo[:,j,i]
        nus=nuss[iiso]

        datamax = data.max()
        nu0_init= nus[np.argmax(data)] 


        noise = data[(velocities<vel_peak-1.) | (velocities>vel_peak+1.)]
        rms = np.std(noise)

        rmss.append(rms)
        molecule_mass=molecule_masses[iiso]
        restfreq=restfreqs[iiso]
        kappa_L = Kappa_line(T_init,iiso)
        f_abund=f_abunds[iiso]
        #initialize Sigma_g so that tau_0 = 1
        Sigma_g_tauone =  (1. / (kappa_L * f_CO * f_abund *  phi(T_init,restfreq,restfreq,vturb_init,molecule_mass)))
        Sigma_g_tauones.append(Sigma_g_tauone)
        
        typicalint=datamax
        if (datamax < (3. *rms)):
            typicalint = 3.*rms

            
        Sigma_g_thin = typicalint/ (bbody(T_init,nu0_init)*kappa_L*f_CO*f_abund*phi(T_init,restfreq,restfreq,vturb_init,molecule_mass))
        Sigma_g_thins.append(Sigma_g_thin)

        if ViewIndividualFits:
            print("iiso ",iiso,"typical int ", typicalint,"Sigma_g_thin",Sigma_g_thin, "f_CO", f_CO, "f_abund", f_abund)

    Sigma_g_init=max(Sigma_g_thins)

    datamin = data.min()
    #Icont_0 = datamin
    #if Icont_0==0:
    #    Icont_0=1e-10
    #Icont_0_lim=(0.5*Icont_0,1.2*Icont_0)

    if ViewIndividualFits:
        print("Initial Conditions")
        print("T_init",T_init)
        print("Sigma_g_init",Sigma_g_init)

        for iiso,restfreq in enumerate(restfreqs):
            nus=nuss[iiso]
            
            v0=v0_init
            nu0=restfreq-(v0/c_light)*restfreq
            data=datas[iiso]
            rms=rmss[iiso]

            initfit=[T_init,vturb_init,Sigma_g_init,v0_init]
            modelij, tau0ij, taus = intensity(nus, initfit[0], nu0, initfit[2], initfit[1],iiso)

            print("iiso ",iiso)
            specobs=np.zeros((len(data),2))
            specmod=np.zeros((len(data),2))
            specobs[:,0]=velocitiess[iiso]
            specobs[:,1]=datas[iiso]
            specmod[:,0]=velocitiess[iiso]
            specmod[:,1]=modelij

            Vtools.Spec([specobs,specmod])




        

    f = lambda Temp,vturb,Sigma_g, v0: master_chi2(nuss, v0, Temp, Sigma_g, vturb, datas, rmss)
    m = Minuit(f, Temp=T_init, vturb=vturb_init, Sigma_g=Sigma_g_init, 
               v0=v0_init,
               error_Temp=1.,
               error_Sigma_g=0.0001,
               error_vturb=100.,
               error_v0=0.01*1E5, 
               limit_Temp=T_limits,
               limit_vturb=(0.0, 1E5),
               limit_Sigma_g=(0., 1.5*Sigma_g_init),
               limit_v0=(v0_init-10.*1E5, v0_init+10.*1E5),
               errordef=1,
               fix_vturb = Fix_vturb,
               #fix_Temp = True,
               #fix_nu0 = True,
    )
    m.migrad()
    errmod = f(m.values['Temp'], m.values['vturb'], m.values['Sigma_g'], m.values['v0'])
    fit = [m.values['Temp'], m.values['vturb'], m.values['Sigma_g'], m.values['v0']];

    isomodelsij=[]
    isotaus0ij=[]
    
    for iiso,restfreq in enumerate(restfreqs):
        nus=nuss[iiso]
        v0=fit[3]
        nu0=restfreq-(v0/c_light)*restfreq
        data=datas[iiso]
        rms=rmss[iiso]
        modelij, tau0ij, taus = intensity(nus, fit[0], nu0, fit[2], fit[1],iiso)
        isomodelsij.append(modelij)
        isotaus0ij.append(tau0ij)
        

    errmodelij = errmod

    if ViewIndividualFits:
        print("Best fit:")
        for aparam in m.values.keys():
            print(aparam,m.values[aparam])
        
        for iiso,restfeq in enumerate(restfreqs):
            nus=nuss[iiso]
            data=datas[iiso]
            specobs=np.zeros((len(data),2))
            specmod=np.zeros((len(data),2))
            specobs[:,0]=velocitiess[iiso]
            specobs[:,1]=datas[iiso]
            specmod[:,0]=velocitiess[iiso]
            specmod[:,1]=isomodelsij[iiso]
            print("iiso ",iiso)
            Vtools.Spec([specobs,specmod])
        
    #return [j,i,fit, model[j,i], tau0[j,i]]
    passout=[j,i,fit, errmodelij, isomodelsij, isotaus0ij]
    #pbar.update(ncores)
    return passout





def exec_optim(inputcubefiles,InputDataUnits='head',maxradius=0.5,moldatafiles=['LAMDAmoldatafiles/molecule_12c16o.inp',],J_up=2,ncores=30,outputdir='./output_iminuit_fixvturb/',ViewIndividualSpectra=False,Fix_vturbulence=False):

    
    global h_P, c_light, k_B,  mp, mH2 
    global cubos
    global nuss, dnus
    global velocitiess #, alpha, Icont_0
    global sigma, molecule_masses, B_21s, g_Jlo, g_Jup, E_los, restfreqs
    global levelenergiess, g_Js
    global f_CO, f_abunds
    global ViewIndividualFits
    global Fix_vturb


    Fix_vturb=Fix_vturbulence
    
    f_CO=1E-4

    ViewIndividualFits=ViewIndividualSpectra
    

    # constants in cgs units
    h_P = const.h.cgs.value
    c_light = const.c.cgs.value
    k_B = const.k_B.cgs.value
    mp = const.m_p.cgs.value
    meanmolweight=2.17
    mH2= meanmolweight * mp

    cubos=[]
    heads=[]
    unitfactors=[]
    print('Opening FITS images ')
    for ainputcubefile in inputcubefiles:
        cubo, head = loadfitsdata(ainputcubefile)
        pixscl = head['CDELT2'] * 3600.
        unitfactor=1.
        if re.search(r"head",InputDataUnits,re.IGNORECASE):
            InputDataUnits=head['BUNIT']
            
        if re.search(r"Jy.*beam",InputDataUnits,re.IGNORECASE):
            print("converting input data units from Jy/beam to CGS/sr, using beam", head['BMAJ'],head['BMIN'])
            omegabeam = (np.pi/(4.*np.log(2.))) * (np.pi/180.)**2 * (head['BMAJ']*head['BMIN'])
            unitfactor = 1E-26 * 1E7 * 1E-4 / omegabeam
            cubo *= unitfactor
        elif re.search(r"Jy.*pix",InputDataUnits,re.IGNORECASE):
            print("converting input data units from Jy/pix to CGS/sr, using pixel", head['CDELT2'])
            omegapix = (np.pi/180.)**2 * (head['CDELT2']**2)
            unitfactor = 1E-26 * 1E7 * 1E-4 / omegapix
            cubo *= unitfactor
        else:
            sys.exit("scale units")
        cubos.append(cubo)
        heads.append(head)
        unitfactors.append(unitfactor)

        
    head=heads[0]
    
    #alpha = 2.3
    #Icont_0 = 0.0   # continuum guess

    if (not re.search(r"\/$",outputdir)):
        outputdir+='/'
        print("added trailing back slack to outputdir")

    #os.system("rm -rf "+outputdir)
    os.system("mkdir "+outputdir)


    maskradpixels = int(maxradius / pixscl)
    print("maskradpixels ",maskradpixels)
    nx=head['NAXIS1']
    ny=head['NAXIS2']
    ii=np.arange(0,nx)
    jj=np.arange(0,ny)
    iis, jjs = np.meshgrid(ii, jj)
    


    tasks=[]
    if ViewIndividualFits:
        for apos in ViewIndividualFits:
            xoffset=apos[1]
            yoffset=apos[0]
            ioff=int(((xoffset/3600.)/head['CDELT1'])+(head['CRPIX1']-1))
            joff=int(((yoffset/3600.)/head['CDELT2'])+(head['CRPIX2']-1))
            print("ioff ",ioff," joff ",joff)
            tasks.append([joff,ioff])
            
    else:
        X0 = ((float(nx)-1.)/2.)
        Y0 = ((float(ny)-1.)/2.)
        irrs=np.sqrt( (iis-X0)**2 + (jjs-Y0)**2)
        mask=np.zeros([ny,nx])
        mask[np.where(irrs < maskradpixels)]=1
        for i in ii:
            for j in jj:
                if (mask[j,i]==1):
                    tasks.append([j,i])
            

    #pbar=tqdm(total=len(tasks))

    MasterMolDatas=[]
    levelenergiess=[]
    E_los=[]
    B_21s=[]
    restfreqs=[]
    molecule_masses=[]
    g_Jss=[] # should all be the same but store for testing 
    dnus=[]
    velocitiess=[]
    nuss=[]


    f_abunds=[]
    isonames=[]
 
    for iiso,amoldatafile in enumerate(moldatafiles):
        MasterMolData = MolData.load_moldata(amoldatafile)
        MasterMolDatas.append(MasterMolData)

        molname=MasterMolData['name']
        f_abund=MolData.molecular_fraction(molname)
        f_abunds.append(f_abund)
        isonames.append(molname)
        
        levelenergies=np.array(MasterMolData['levelenergies'])
        levelenergiess.append(levelenergies)
        
        g_Js = np.array(MasterMolData['g_Js'])
        g_Jss.append(g_Js)

        levelJs = MasterMolData['levelJs']
        levelnumber = MasterMolData['levelnumbers']
    
        iJ_up=levelJs.index(J_up)
        iJ_lo=iJ_up-1
        g_Jup=g_Js[iJ_up]
        g_Jlo=g_Js[iJ_lo]
        n_up=levelnumber[iJ_up]

        E_lo=levelenergies[iJ_lo]
        E_los.append(E_lo)
    
        alltransitions=MasterMolData['transitions'].keys()
        for itransition in alltransitions:
            thistransition=MasterMolData['transitions'][itransition]
            if (thistransition['nlevelup'] == n_up):
                Einstein_A=thistransition['Einstein_A']
                restfreq=thistransition['restfreq']
                restfreqs.append(restfreq)
                Einstein_B21=thistransition['Einstein_B21']
                B_21=Einstein_B21
                B_21s.append(B_21)
                break
        
        molecule_mass = MasterMolData['molecularmass']
        molecule_masses.append(molecule_mass)

        if ViewIndividualFits:
            print("restfreq :",restfreq)
            print("Einstein_A :",Einstein_A)
            print("molecule_mass ",molecule_mass)



        print("using header number",iiso)
        ahead=heads[iiso]
        dnu = ahead['CDELT3']
        len_nu = ahead['NAXIS3']
        nus= ahead['CRVAL3']+(np.arange(ahead['NAXIS3'])-ahead['CRPIX3']+1)*ahead['CDELT3']
        velocities = -(nus-restfreq)*c_light*1E-5/restfreq # velocities in km/s

        nuss.append(nus)
        velocitiess.append(velocities)
        dnus.append(dnu)

    print("Molecule names:",isonames)
    print("Molecule fractions:",f_abunds)
    
        
    ndim = head['NAXIS1']
    Temperature = np.zeros((ndim,ndim))
    tau0 = np.zeros((ndim,ndim))
    Sigma_g_im = np.zeros((ndim,ndim))
    Turbvel = np.zeros((ndim,ndim))
    velo_centroid = np.zeros((ndim,ndim))
    errmodel = np.zeros((ndim,ndim))
    dust = np.zeros((ndim,ndim,cubo.shape[0]))

    nisos=len(inputcubefiles)
    
    models = []
    isotau0s = []
    for iiso in list(range(nisos)):
        model = np.zeros(cubo.shape)
        isotau0 = np.zeros((ndim,ndim))
        models.append(model)
        isotau0s.append(isotau0)
        
    mom2 = np.zeros((ndim,ndim))
    
    
    #pool = Pool(ncores)
    #todo=pool.map(parspar, tasks)
    #pbar.close()


    with Pool(ncores) as pool:
        Pooloutput = list(tqdm(pool.imap(parspar, tasks), total=len(tasks)))
        pool.close()
        pool.join()



    print("Done whole pool")


    for ls in Pooloutput:
        if len(ls)!=6:
            continue
        j = ls[0]
        i = ls[1]
        fit = ls[2]
        Temperature[j,i] = fit[0]
        Sigma_g_im[j,i] = fit[2]
        Turbvel[j,i] = fit[1]
        velo_centroid[j,i]= fit[3]*1E-5

                     
        rettau0s  = ls[-1]
        retmodels = ls[-2]
        for iiso in list(range(nisos)):
            models[iiso][:,j,i]=retmodels[iiso]
            isotau0s[iiso][j,i]=rettau0s[iiso]
                  
        errmodel[j,i] = ls[-3]




    punchout=[]
    punchout.append({'data':Sigma_g_im,'BUNIT':'g/cm2','BTYPE':'MassColumn','outfile':'Sigma_g.fits'}) 
    punchout.append({'data':Turbvel,'BUNIT':'cm/s','BTYPE':'Velocity','outfile':'vturb.fits'}) 
    punchout.append({'data':Temperature,'BUNIT':'K','BTYPE':'Temperature','outfile':'temperature.fits'}) 
    for iiso in list(range(nisos)):
        punchout.append({'data':isotau0s[iiso],'BUNIT':'N/A','BTYPE':'OpticalDepth','outfile':'tau0_'+isonames[iiso]+'.fits'}) 

    punchout.append({'data':velo_centroid,'BUNIT':'km/s','BTYPE':'Velocity','outfile':'velocentroid.fits'}) 
    punchout.append({'data':errmodel,'BUNIT':'erg/s/cm2/Hz/sr','BTYPE':'Intensity','outfile':'fiterror.fits'}) 
    
    for apunchout in punchout:
        dataout=np.nan_to_num(apunchout['data'])
        rout=pf.PrimaryHDU(dataout)
        headout = deepcopy(head)
        headout['BUNIT']=apunchout['BUNIT']
        headout['BTYPE']=apunchout['BTYPE']
        rout.header=headout
        rout.writeto(outputdir+apunchout['outfile'], overwrite=True)

    StoreModels=True
    if StoreModels:
        for iiso in list(range(nisos)):
            unitfactor=unitfactors[iiso]
            amodel=models[iiso]/unitfactor
            ahead=heads[iiso]
            rout=pf.PrimaryHDU(amodel)
            rout.header=ahead
            rout.writeto(outputdir+'model_'+isonames[iiso]+'.fits', overwrite=True)
        

    return


