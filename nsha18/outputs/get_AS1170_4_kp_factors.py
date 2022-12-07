from __future__ import unicode_literals
#from openquake.nrmllib.hazard.parsers import HazardCurveXMLParser
from numpy import array, arange, exp, log, interp, vstack, mean, around
from tools.oq_tools import return_annualised_haz_curves
from gmt_tools import cpt2colormap
from os import path, mkdir
import matplotlib.pyplot as plt
import warnings, sys
#reload(sys) # for unicode chars
#sys.setdefaultencoding("latin-1")
warnings.filterwarnings("ignore")
plt.rcParams['pdf.fonttype'] = 42

import matplotlib as mpl
mpl.style.use('classic')

###############################################################################
# read config file
###############################################################################

period = 'PGA'
'''
conf_file = sys.argv[1]

# get paths for input files
lines = open(conf_file).readlines()
hazcurvefile1 = lines[0].split('=')[-1].strip()
job1 = hazcurvefile1.split('/')[0]
hazcurvefile2 = lines[1].split('=')[-1].strip()
job2 = hazcurvefile2.split('/')[0]
outputdir     = lines[2].split('=')[-1].strip()
sitelistfile  = lines[3].split('=')[-1].strip()
period        = lines[4].split('=')[-1].strip()

# check to see if exists
if path.isdir(outputdir) == False:
    mkdir(outputdir)
'''
###############################################################################
# read OQ data
###############################################################################


###############################################################################
# parse site file
###############################################################################

sitelistfile = path.join('X:\\NSHA2018\\shared\\nsha_cities.csv')
#sitelistfile = 'nsha_cities.csv'
lines = open(sitelistfile).readlines() 
places = []
place_lat = []
place_lon = []

for line in lines:
    dat = line.strip().split(',')
    place_lon.append(float(dat[0]))
    place_lat.append(float(dat[1]))
    places.append(dat[2])
"""
###############################################################################
# parse 2007 AS1170.4 PGA
###############################################################################

as1170file = path.join('2007_AS1170.4', 'AS1170.4_pga.csv')
place2007 = []
pga2007 = []

lines = open(as1170file).readlines()

for line in lines:
    dat = line.strip().split(',')
    place2007.append(dat[0])
    pga2007.append(float(dat[1]))

###############################################################################
# parse 2012 NSHA PGA
###############################################################################

nsha12file = path.join('2012_NSHA', 'nsha13_localities_pga.csv')
lon2012 = []
lat2012 = []
place2012 = []
pga2012_10 = []
pga2012_02 = []

lines = open(nsha12file).readlines()[1:]

for line in lines:
    dat = line.strip().split(',')
    lon2012.append(float(dat[0]))
    lat2012.append(float(dat[1]))
    pga2012_10.append(float(dat[2]))
    pga2012_02.append(float(dat[3]))
    place2012.append(dat[-1])
"""    
###############################################################################
# parse 2018 NSHA PGA from OQ
###############################################################################
#hazcurvefile = path.join('Complete_Model', 'hazard_curve-mean_1.csv')
print('\n!!!! Final NSHA18 hazard curves !!!!\n')
#print('tallen/NSHA2018/source_models/complete_model/final/results_fractilesUHS/hazard_curve-mean-PGA_1.csv')

hazcurvefile = 'X:\\NSHA2018\\source_models\\complete_model\\final\\results_fractilesUHS\\hazard_curve-mean-PGA_1.csv'

# get annualize the curves.
siteDict, imls, investigation_time = return_annualised_haz_curves(hazcurvefile)

###############################################################################
# write OQ & Frisk data to file
###############################################################################
i = 1
ii = 0
yhaz = 1./475.

# set NBCC probabilities
def get_prob_from_pc(pc, yr):
    p0 = 1 - (pc / 100.)
    n = -log(p0)
    r = n / yr
    
    return 1/r, r
    
#kpprobs = array([0.046052, 0.0277258, 0.01386, 0.01, 0.00446, 0.00211, 0.0.00103, 0.000506, 0.000404, 0.000201, 0.0001002])
# get kp probs
#num_yrs = 50.
#kpprobs = get_prob_from_pc(percentages, 50.)[1]

'''
In [9]: 1./kpprobs
array([   21.7147241 ,    26.35573944,    36.06737602,    54.5678334 ,
          72.13475204,    97.88075945,   140.1836626 ,   224.07100589,
         474.56107905,   974.78628731,  1224.82991308,  1641.53975526,
        2474.91582263,  4974.95812367,  9974.97911442])
'''
from hazard_tools import get_percent_chance_from_return_period, get_probability_from_percent_chance

return_periods = array([10000, 7500, 5000, 4000, 3000, 2500, 2475, 2000, 1500, 1000, 800, 500, 475, 250, 200, 100, 50, 25, 20])
investigation_time = 50.
percent_chance = get_percent_chance_from_return_period(return_periods, investigation_time)

return_period, probability = get_probability_from_percent_chance(percent_chance, investigation_time)
kpprobs = probability
# loop thru AS1170.4 cities
i = 0

kpcurves = []

capitals = ['Adelaide', 'Brisbane', 'Canberra', 'Darwin', 'Hobart', 'Melbourne', 'Perth', 'Sydney']

target_imt = 'poe'
imlkey = target_imt+'_imls'
probkey = target_imt+'_probs_annual'
     
# now loop through OQ runs and match lon/lats
for sd2018 in siteDict:
    matchPlace = False
    # get city name
    for place, pllo, plla in zip(places, place_lon, place_lat):
        if round(sd2018['lon'],2) == round(pllo,2) \
            and round(sd2018['lat'],2) == round(plla,2):
            matchPlace = True
            curve_city = place.strip()
            print(place)
    
    # now check if city in list of capitals
    if matchPlace == True:
        
        for capital in capitals:
            #print(capital
            if capital.strip() == curve_city:
        
                # interp to as1170.4 probs
                haz2018 = exp(interp(log(kpprobs), log(sd2018[probkey][::-1]), log(imls[::-1]))) #[::-1]
                #haz2018 = exp(interp(log(kpprobs), log(curve2018[::-1]), log(imls[::-1])))#[::-1]
                kpfacts = haz2018 / haz2018[9] # assumes 5000-year max
                kpfacts = haz2018 / haz2018[11] # assumes 10,000-year max
                
                tmp = {'place': curve_city, 'pga2018': haz2018, 'kpfact':kpfacts}
                       
                # add loc dict
                kpcurves.append(tmp)
                
                # export city kp factor
                kptxt = ''.join(('KP FACTORS FOR ', curve_city.upper(),'\nRETURN_PERIOD,KP_FACTOR\n'))
                for rp, kp in zip(return_periods, kpfacts):
                    kptxt += ','.join((str(rp), str('%0.2f' % kp))) + '\n'
                    
                f = open(path.join('kp_factors', curve_city+'_kpfacts.csv'), 'w')
                f.write(kptxt)
                f.close()
                  
# begin plotting
ncols = len(capitals)+1
#cptfile = '\\prod.lan\\active\\ops\\community_safety\\ehp\\georisk_earthquake\\hazard\\DATA\\cpt\\gay-flag-1978.cpt'
cptfile = 'Z:\\DATA\\cpt\\gay-flag-1978.cpt'
cmap, zvals = cpt2colormap(cptfile, ncols)
cptcmap = (cmap(arange(ncols)))

kprp = 1. / kpprobs

fig = plt.figure(1, figsize=(6,6))
plt.plot([2500, 2500], [0, 12], '--', lw=0.75, c='0.2')
for kpc, c in zip(kpcurves, cptcmap):
    plt.plot(kprp, kpc['kpfact'], '-', lw=1.5, c=c, label=kpc['place'])

# add AS1170.4
asyr = array([20, 25, 50, 100, 200, 250, 500, 800, 1000, 1500, 2000, 2500])
askp = array([0.2, 0.25, 0.35, 0.5, 0.7, 0.75, 1., 1.25, 1.3, 1.5, 1.7, 1.8])
plt.plot(asyr, askp, 'k--', lw=1.5, label='AS1170.4-2007')
    
plt.legend(loc=2, fontsize=10.5)
plt.xlim([0, 10000])
plt.ylim([0, 9])

plt.xlabel('Return Period (Years)')
plt.ylabel('Probability Factor '+r'$k_p$')
plt.grid()
#plt.ylabel(r'\textit{k_p}')

plt.savefig('2018_nsha_kp_factors.png', fmt='png', bbox_inches='tight', dpi=300)

# loop thru all sites and get average kp
nat_kp_2500 = []
# now loop through OQ runs and match lon/lats
for sd2018 in siteDict:
       # interp to as1170.4 probs
       haz2018 = exp(interp(log(kpprobs), log(sd2018[probkey][::-1]), log(imls[::-1]))) #[::-1]
      
       # divide 2500 yr value by 500 year value
       nat_kp_2500.append(haz2018[3] / haz2018[9])
 
nat_kp_2500 = array(nat_kp_2500)

print('\nMean kp factor =', mean(nat_kp_2500))
plt.show()


"""
# loop thru localities and calc differences from 2007
outtxt = 'PLACE,2007PGA_P0.0021,2012PGA_P0.0021,2018PGA_P0.0021,2007/2012 %DIFF,2007/2018 %DIFF, 2012/2018 %DIFF,2007PGA_P0.0004,2012PGA_P0.0004,2018PGA_P0.0004\n'
for hc in hazcmp:
    # calc % difference
    numer = hc['pga2012'][0] - hc['pga2007']
    denom = mean(vstack((hc['pga2007'], hc['pga2012'][0])), axis=0)
    pcdiff_07_12 = 100. * (numer / denom)
    
    numer = hc['pga2018'][0] - hc['pga2007']
    denom = mean(vstack((hc['pga2007'], hc['pga2018'][0])), axis=0)
    pcdiff_07_18 = 100. * (numer / denom)
    
    numer = hc['pga2018'][0] - hc['pga2012'][0]
    denom = mean(vstack((hc['pga2012'][0], hc['pga2018'][0])), axis=0)
    pcdiff_12_18 = 100. * (numer / denom)

    # set % diff text
    outtxt += ','.join((hc['place'], str('%0.2f' % hc['pga2007']), \
                        str('%0.3f' % hc['pga2012'][0]), str('%0.3f' % hc['pga2018'][0]), \
                        str('%0.2f' % pcdiff_07_12), str('%0.2f' % pcdiff_07_18), \
                        str('%0.2f' % pcdiff_12_18), str('%0.3f' % (hc['pga2007']*1.8)), \
                        str('%0.3f' % hc['pga2012'][1]), str('%0.3f' % hc['pga2018'][1]))) + '\n'
            
# write to file
csvfile =  path.join('as1170.4_hazard_comparison_test.csv')
f = open(csvfile, 'wb')
f.write(outtxt)
f.close()
"""           
