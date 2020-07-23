from __future__ import unicode_literals
from numpy import array, exp, log, interp, ceil, mean, vstack, around, arange, loadtxt
from tools.oq_tools import return_annualised_haz_curves
import matplotlib.pylab as plt
from os import path
import warnings, sys
from gmt_tools import cpt2colormap 
from misc_tools import get_mpl2_colourlist
import matplotlib as mpl
mpl.style.use('classic')

#reload(sys) # for unicode chars
#sys.setdefaultencoding("latin-1")
#warnings.filterwarnings("ignore")

##############################################################################
# parse param file
##############################################################################
paramfile = 'multi_fault_hazcurves.conf'
lines = open(paramfile).readlines()
hazpaths = []
modnames = []
outfile = lines[0].strip()
for line in lines:
    modnames.append(line.strip().split(';')[0])
    hazpaths.append(line.strip().split(';')[1])

period = '1.0'

###############################################################################
# parse site file
###############################################################################
sitelistfile = '/Users/trev/Documents/NRCan/2015_National_Hazard/2015_gsc_nshm/shared/swcan_sites_cities.csv'
lines = open(sitelistfile).readlines()
places = []
place_lat = []
place_lon = []

for line in lines:
    dat = line.strip().split(',')
    place_lon.append(float(dat[0]))
    place_lat.append(float(dat[1]))
    places.append(dat[2])
    
###############################################################################
# make colormap
###############################################################################
'''
cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/GMT_no_green.cpt'
#cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/temperature.cpt'
cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/gay-flag-1978.cpt'
#cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/qual-dark-06.cpt'
'''
ncolours = 5
#cmap, zvals = cpt2colormap(cptfile, ncolours, rev=False)
cmap = plt.cm.get_cmap('hsv', ncolours)
cs = (cmap(arange(ncolours)))
cs = get_mpl2_colourlist()

###############################################################################
# read frisk data
###############################################################################

tmpper = period.replace('.','')
friskfile = 'OQ_CISbGMPE_1.0_000thperc.edit.sol'
friskhaz = loadtxt(friskfile, delimiter=',', skiprows=1, usecols = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11))
friskprobs = [0.02, 0.01375, 0.01, 0.00445, 0.0021, 0.001, 0.0005, 0.000404, 0.0002, 0.0001]

# get place names
places = []
place_lon = []
place_lat = []
lines = open(friskfile, encoding = 'latin-1').readlines()[1:]
for line in lines:
    dat = line.strip().split(',')
    places.append(dat[-1])
    place_lon.append(float(dat[0]))
    place_lat.append(float(dat[1]))

###############################################################################
# plot frisk data
###############################################################################

i = 1
fig = plt.figure(i, figsize=(14, 10))
yhaz = 1./2475.

###############################################################################
# read OQ data
###############################################################################
jj = 0
for j, hazcurvefile in enumerate(hazpaths):   
    print(hazcurvefile)
    ii = 0
    siteDict, imls, investigation_time = return_annualised_haz_curves(hazcurvefile)
    
    # loop thru sites in hazcurvefile
    for place, plon, plat, fh in zip(places, place_lon, place_lat, friskhaz):
        for sd in siteDict:
            curlon = sd['lon']
            curlat = sd['lat']
            curve = sd['poe_probs_annual']
            
            ###############################################################################
            # plt OQ & Frisk data
            ###############################################################################
        
            # loop thru OQ curves
        
            if around(plon, decimals=2) == around(sd['lon'], decimals=2) \
               and around(plat, decimals=2) == around(sd['lat'], decimals=2):
                print(place)
                ii += 1 # subplot index
                plt.subplot(2,3,ii)
                
                # plt FRISK curve
                if j == 0:
                    h2 = plt.semilogy(fh[6:], friskprobs[4:], '--', c='k', lw=3, label='GSCFRISK')
                
                # plt OQ curves    
                h1 = plt.semilogy(imls, curve, '-', c=cs[j], lw=1.5, label=modnames[j])
                
                # now make pretty
                if j == 3:
                    # get x lims from frisk haz
                    thaz = exp(interp(log(1e-4), log(friskprobs[::-1]), log(fh[2:][::-1])))
                    
                    # round to nearest 0.1
                    if thaz < 0.1:
                        xmax = ceil(thaz / 0.01) * 0.01
                    else:
                        xmax = ceil(thaz / 0.1) * 0.1
                    plt.xlim([0, xmax])
                    plt.ylim([1e-4, 3e-3])
                    
                    plt.title(place)#.encode('utf8'))
                    plt.grid(which='both')
                    plt.semilogy([0, 2.5], [yhaz, yhaz], 'k--')
                    plt.xticks(rotation=45)

                    if ii == 3:
                        plt.legend(fontsize=9)
                    
                    if ii == 1 or ii == 4:
                        plt.ylabel('Annual Probabability of Exceedance', fontsize=14)
                        
                    if ii == 4 or ii == 5 or ii == 6:
                        plt.xlabel(' '.join(('Mean','SA['+period+']', 'Hazard (g)')), fontsize=14)
                        
                    if ii == 6:
                      plt.savefig('scalrel_hazard.png', format='png',bbox_inches='tight', dpi=800)
                     
                      i += 1
                      ii = 0
                      #fig = plt.figure(i, figsize=(14, 10)
                  
plt.subplots_adjust(bottom=0.1)
plt.savefig('scalrel_hazard.png', format='png',bbox_inches='tight', dpi=300)
plt.show()
