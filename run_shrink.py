import sys
import os
import re

include_path='/Users/simon/common/python/include/'
sys.path.append(include_path)
import ShrinkCanvas


ShrinkCanvas.Zoom('output_iminuit_multiso_dev/model_CO.fits',zoom_area=0.6)
ShrinkCanvas.Zoom('output_iminuit_multiso_dev/model_13C16O.fits',zoom_area=0.6)
ShrinkCanvas.Zoom('output_iminuit_multiso_dev/model_C18O.fits',zoom_area=0.6)



inputcubefiles = [
    '/strelka_ssd/simon/HD135344B/red/tclean_contsubHD135344Bbriggs2.0_12CO.fits',  
    '/strelka_ssd/simon/HD135344B/red/tclean_contsubHD135344Bbriggs2.0_13CO.fits',
    '/strelka_ssd/simon/HD135344B/red/tclean_contsubHD135344Bbriggs2.0_C18O.fits']

for afile in inputcubefiles:
    ShrinkCanvas.Zoom(afile,zoom_area=0.6)
    
