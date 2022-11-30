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
                          adjust_gmm_with_nga_east, gaull1990_gsim, crustal_gsims
from data_fmt_tools import return_sta_data
from mapping_tools import distance
from misc_tools import listdir_extension
from gmt_tools import cpt2colormap, remove_last_cmap_colour
from misc_tools import get_mpl2_colourlist
from shakemap_tools import parse_faultdat, make_fault_mesh

# set eq deets
ztor = 4. # from JLG
rake = 0. # GA CMT
dip  = 83.
mag = 5.9
dep = 10.
vs30s = [270., 400., 560., 760., 1100.]

# get fault deets
faultfile = 'woods_point_fault.txt'
res = 0.1 # km
faultdat = make_fault_mesh(faultfile, res)

# parse chimney file
csvfile = 'Chimney_v2.csv'
data = loadtxt(csvfile, skiprows=1, delimiter=',')
lats = data[:,1]
lons = data[:,0] # now added
n_chim = data[:,2]
	
# get distance for each of the chimneys
epilon = 146.402
epilat = -37.506
hypdep = 12.7

rjbs = []
rrups = []
rhyps = []
for i, chim in enumerate(n_chim):
    minrjb = 9999.
    minrrup = 9999.
    for flo, fla, fde in zip(faultdat['lon'], faultdat['lat'], faultdat['dep']):
        rngkm = distance(lats[i], lons[i], fla, flo)[0]
        if rngkm < minrjb:
            minrjb = rngkm
            minrrup = sqrt(rngkm**2 + fde**2)
    rjbs.append(minrjb)
    rrups.append(minrrup)
    rhyps.append(sqrt(distance(lats[i], lons[i], fla, flo)[0]**2 + hypdep**2))

for vs30 in vs30s:
    print(vs30)
    txt = '#CHIMNEY,RHYP,RJB,RRUP,Gea90SEA,Gea90WA,AB06,CY08SWISS,Sea09NC,Sea09YC,A12,Bea14,CY14,NGA-E\n'
    
    for chim, rjb, rhyp, rrup in zip(n_chim, rjbs, rhyps, rrups):
        print(chim)
        #rrup = sqrt(rjb**2 + ztor**2)
        #rhyp = rrup = sqrt(rjb**2 + dep**2) # assume for now as need site locations
        
        # get SCR GMMs
        Tea02imt, C03imt, AB06imt, AB11imt, CY08imtSWISS, Sea09imt_SS14, Sea09YCimt_SS14, Pea11imt, A12imt, A12imt_SS14, Bea14imt \
                 = scr_gsims(mag, dep, ztor, dip, rake, rrup, rjb, vs30)
        
        # get AC GMMs
        Bea97imt, Zea06imt, CB08imt, CY08imt, Bea11imt, BA11imt, Aea14imt, Bea14imt, CB14imt, CY14imt \
                 = crustal_gsims(mag, dep, ztor, dip, rake, rrup, rjb, vs30)
        #Tea19 = tang2019_cam_gsim(mag, dep, rrup, vs30)
        
        G90WAimt, G90SEAimt, G90INDimt, G90WA_PGVimt, G90SEA_PGVimt, G90IND_PGVimt = gaull1990_gsim(mag, dep, rhyp)
        
        # get mean USGS NGA-E
        nga_e_imt = nga_east_mean(mag, dep, dip, rake, rrup, vs30, ztor=4.0)[0]
        
        #plt.loglog(nga_e_imt['per'], exp(nga_e_imt['sa']),'--' , lw=1.5, color=cs[5]) # for testing amp factors
        nga_e_imt_hold = nga_e_imt
        
        # adjust NGA-E from 3000 -> target
        nga_e_imt = adjust_gmm_with_nga_east(nga_e_imt, vs30)
        
        txt += ','.join((str('%.0f' % chim), str(rhyp), str(rjb), str(rrup), str(exp(G90SEAimt['pga'][0][0])), str(exp(G90WAimt['pga'][0][0])), \
                         str(exp(AB06imt['pga'][0][0])), str(exp(CY08imtSWISS['pga'][0][0])), str(exp(Sea09imt_SS14['pga'][0][0])), str(exp(Sea09YCimt_SS14['pga'][0][0])), \
                         str(exp(A12imt_SS14['pga'][0][0])), str(exp(Bea14imt['pga'][0][0])), str(exp(CY14imt['pga'][0][0])), str(exp(nga_e_imt['sa'][0])))) + '\n'
                         
    
    outfile = 'Chimney_RJB_GMPE_'+str('%0.0f' % vs30)+'.csv'
    f = open(outfile, 'w')
    f.write(txt)
    f.close()
    
