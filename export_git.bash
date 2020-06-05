

destination=/Users/simon/common/python/pylinefit/
origin=/strelka_ssd/simon/HD135344B/linefit/

rsync -va $origin\README                            $destination

rsync -va $origin\MolData.py                        $destination
rsync -va $origin\gen_Dmoments.py		    $destination
rsync -va $origin\run_shrink.py			    $destination
rsync -va $origin\gen_Dmoments_residuals.py	    $destination
rsync -va $origin\Linefit_iminuit.py		    $destination
rsync -va $origin\run_linefit.py		    $destination
rsync -va $origin\run_summary.py 		    $destination
rsync -va $origin\gen_comparspectra.py		    $destination
rsync -va $origin\LAMDAmoldatafiles                 $destination
rsync -va $origin\export_git.bash                   $destination
rsync -va $origin\SummaryLineFit.py                   $destination

rsync -va /strelka_ssd/simon/HD135344B/red/tclean_contsubHD135344Bbriggs2.0_12CO_z.fits    $destination\example_data/
rsync -va /strelka_ssd/simon/HD135344B/red/tclean_contsubHD135344Bbriggs2.0_13CO_z.fits	   $destination\example_data/
rsync -va /strelka_ssd/simon/HD135344B/red/tclean_contsubHD135344Bbriggs2.0_C18O_z.fits    $destination\example_data/

rsync -va /Users/simon/common/python/include/DGaussMoments.py    $destination
rsync -va /Users/simon/common/python/include/Cube2Im.py		 $destination
rsync -va /Users/simon/common/python/include/Resamp.py		 $destination
rsync -va /Users/simon/common/python/include/Vtools.py           $destination


tar cvfz $origin\example_output.tgz  $origin\/output_iminuit_multiso

rsync -va $origin\example_output.tgz             $destination
rsync -va $origin\fig_comparespectra*pdf $destination





										    
