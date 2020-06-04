import os
import sys
import re
include_path=os.environ['HOME']+'/common/python/include/'
sys.path.append(include_path)
import  DGaussMinuit 
import  DGaussMoments 

import time
#from time import gmtime, strftime


inputcubefiles = [
    '/strelka_ssd/simon/HD135344B/red/tclean_contsubHD135344Bbriggs2.0_12CO.fits',  
    '/strelka_ssd/simon/HD135344B/red/tclean_contsubHD135344Bbriggs2.0_13CO.fits',
    '/strelka_ssd/simon/HD135344B/red/tclean_contsubHD135344Bbriggs2.0_C18O.fits']



start_time=time.time()
print("start Curve_fit:",time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()))


filein = '/strelka_ssd/simon/HD135344B/guvmem_runs/belka_results/12CO21/smooth_model_cube_lS0.003_lL0.0.fits'  
workdir='12CO_sgauss_contsub_guvmem_lS0.003_lL0.0_smooth'
os.system("rm -rf "+workdir)
os.system('mkdir '+workdir)
DGaussMoments.exec_Gfit(filein,workdir,wBaseline=False,n_cores=30,zoom_area=1.0,Noise=2.5E-3,Clip=True,DoubleGauss=False,StoreModel=False,Randomize2ndGauss=False,ShrinkCanvas=True,UseCommonSigma=False)

end_time=time.time()
print( "curve_fit done in (elapsed time):",   end_time - start_time)

