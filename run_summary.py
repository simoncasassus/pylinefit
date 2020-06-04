
import SummaryLineFit

workdirs=['output_iminuit_multiso/',]
vsyst= 7.2
for aworkdir in workdirs:
        fileout = aworkdir+'fig_summary.pdf'
        fix_vturb=False
        SummaryLineFit.exec_summary(aworkdir,fileout,vsyst=vsyst, vrange=10.,fix_vturb=fix_vturb, WCont=False, Zoom=False)
