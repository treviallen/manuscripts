from __future__ import unicode_literals
from numpy import array, exp, log, interp, ceil, around, arange
from tools.oq_tools import return_annualised_haz_curves
from os import path, mkdir, getcwd
import warnings, sys
from sys import argv
from hazard_tools import get_percent_chance_from_return_period, get_probability_from_percent_chance

def interp_hazard_curves(investigation_time, interp_poe, poe_imls):
    from os import path, mkdir
    from numpy import array, exp, log, interp
    
    # set grid return periods for hazard curve
    return_periods = ['100', '250', '475', '500', '800', '1000', '1500', '2000', '2475', \
                      '2500', '3000', '5000']
    
    probabilityArray = []
    percent_chanceArray = []
    interphazArray = []
    return_periodsArray = []
    
    for return_period in return_periods:
        percent_chanceArray.append(get_percent_chance_from_return_period(float(return_period), investigation_time))
        return_period_num, probability = get_probability_from_percent_chance(percent_chanceArray[-1], investigation_time)
        
        probabilityArray.append(probability)
        
        interphaz = exp(interp(log(probability), log(interp_poe[::-1]), log(poe_imls[::-1])))

        interphazArray.append(interphaz)
        
    return array(interphazArray), return_periods, array(probabilityArray)

###############################################################################
# start main code here
###############################################################################


# string format, either 'PGA' or 'SA(0.1)' for example
period = argv[1]

reload(sys) # for unicode chars
sys.setdefaultencoding("latin-1")
warnings.filterwarnings("ignore")

###############################################################################
# set haz curve file and variables
###############################################################################

# set file
if getcwd().startswith('/nas'):
    hazcurvefile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/source_models/complete_model/final/results_fractilesUHS/hazard_curve-mean-'+period+'_1.csv'
    sitelistfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/shared/nsha_cities.csv'
else:
    hazcurvefile = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/source_models/complete_model/final/results_fractiles/hazard_curve-mean-'+period+'_1.csv'
    sitelistfile = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/shared/nsha_cities.csv'


###############################################################################
# parse site file
###############################################################################

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
# parse first job file to define plotting order
###############################################################################

# get data from first job
siteDict1, imls, investigation_time = return_annualised_haz_curves(hazcurvefile)

# loop thru sites in first job file and plot
interphazArrays = []
writePlaces = []
writeLon = []
writeLat = []
    
for sd1 in siteDict1:
    
    ###############################################################################
    # loops thru places to get title - check if want to plot
    ###############################################################################
    for place, plon, plat in zip(places, place_lon, place_lat):
        if around(plon, decimals=2) == around(sd1['lon'], decimals=2) \
           and around(plat, decimals=2) == around(sd1['lat'], decimals=2):
            if not place.endswith('max'):
                # interpolate curve to standard return periods
                interphazArray, return_periods, probabilityArray = \
                    interp_hazard_curves(investigation_time, sd1['poe_probs_annual'], imls)
                    
                interphazArrays.append(interphazArray)
                
                writePlaces.append(place)
                writeLon.append(sd1['lon'])
                writeLat.append(sd1['lat'])
            
# set header lines
#outtxt = 'RETURN PERIODS,LON,LAT,'+','.join(return_periods) + '\n'
outtxt = 'ANNUAL_PROBABILITY,LON,LAT,'+','.join(['P'+str('%0.6f' % x) for x in probabilityArray]) + '\n'

for place, interpHaz, lon, lat in zip(writePlaces, interphazArrays, writeLon, writeLat):
    outtxt += ','.join((place, str('%0.2f' % lon), str('%0.2f' % lat))) \
              +','+','.join([str('%0.4e' % x) for x in interpHaz]) + '\n'

# set output file
outcsv = path.join('annual_haz_curves', '_'.join(('nsha18_haz_curves', period +'.csv')))

# write file
f = open(outcsv, 'wb')
f.write(outtxt)
f.close()

    
