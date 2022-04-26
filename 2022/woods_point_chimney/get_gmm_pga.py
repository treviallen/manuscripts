import warnings
warnings.filterwarnings('ignore')


    
'''
start main
'''
# start of plotting routine
from numpy import array, arange, sqrt, exp, log, unique, argsort, isnan, loadtxt
from sys import argv
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
import matplotlib as mpl
mpl.style.use('classic')

from calc_oq_gmpes import inslab_gsims, scr_gsims, tang2019_cam_gsim, \
                          nga_east_mean, get_station_vs30, adjust_gmm_with_SS14, \
                          adjust_gmm_with_nga_east, gaull1990_gsim
from data_fmt_tools import return_sta_data
from mapping_tools import distance
from misc_tools import listdir_extension
from gmt_tools import cpt2colormap, remove_last_cmap_colour
from misc_tools import get_mpl2_colourlist

# set eq deets
ztor = 4. # from JLG
rake = 0. # GA CMT
dip  = 83.
mag = 5.9
dep = 10.
vs30 = 560.

# parse chimney file
csvfile = 'Chimney_RJB_GMPE.csv'
data = loadtxt(csvfile, skiprows=1, delimiter=',')
rjbs = data[:,1]
rhyps = data[:,2] # now added
n_chim = data[:,0]

txt = '#CHIMNEY,RHYP,RJB,Gea90SEA,Gea90WA,AB06,Sea09NC,Sea09YC,A12,Bea14,NGA-E\n'

for chim, rjb, rhyp in zip(n_chim, rjbs, rhyps):
    print(chim)
    rrup = sqrt(rjb**2 + ztor**2)
    #rhyp = rrup = sqrt(rjb**2 + dep**2) # assume for now as need site locations
    
    Tea02imt, C03imt, AB06imt, AB11imt, Sea09imt_SS14, Sea09YCimt_SS14, Pea11imt, A12imt, A12imt_SS14, Bea14imt \
             = scr_gsims(mag, dep, ztor, dip, rake, rrup, rjb, vs30)
    
    #Tea19 = tang2019_cam_gsim(mag, dep, rrup, vs30)
    
    G90WAimt, G90SEAimt, G90INDimt, G90WA_PGVimt, G90SEA_PGVimt, G90IND_PGVimt = gaull1990_gsim(mag, dep, rhyp)
    
    # adjust some GMMs using Seyhan & Stewart 2014
    #A12imt = adjust_gmm_with_SS14(A12imt, 820., vs30) - corrected in gsim code
    #Sea09imt = adjust_gmm_with_SS14(Sea09imt, 865., vs30)
    
    
    
    # get mean USGS NGA-E
    nga_e_imt = nga_east_mean(mag, dep, dip, rake, rrup, vs30)
    #plt.loglog(nga_e_imt['per'], exp(nga_e_imt['sa']),'--' , lw=1.5, color=cs[5]) # for testing amp factors
    nga_e_imt_hold = nga_e_imt
    
    # adjust NGA-E from 3000 -> target
    nga_e_imt = adjust_gmm_with_nga_east(nga_e_imt, vs30)
    
    
    txt += ','.join((str(chim), str(rhyp), str(rjb), str(exp(G90SEAimt['pga'][0][0])), str(exp(G90WAimt['pga'][0][0])), \
                     str(exp(AB06imt['pga'][0][0])), str(exp(Sea09imt_SS14['pga'][0][0])), str(exp(Sea09YCimt_SS14['pga'][0][0])), \
                     str(exp(A12imt_SS14['pga'][0][0])), str(exp(Bea14imt['pga'][0][0])), str(exp(nga_e_imt['sa'][0])))) + '\n'
                     

outfile = 'Chimney_RJB_GMPE_'+str('%0.0f' % vs30)+'.csv'
f = open(outfile, 'w')
f.write(txt)
f.close()
    
