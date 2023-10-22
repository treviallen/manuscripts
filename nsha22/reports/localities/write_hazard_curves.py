from tools.oq_tools import return_annualised_haz_curves
from numpy import arange, around, array, interp, unique
from os import path, getcwd, sep
from sys import argv

'''
example run:
    run plt_uhs_curves.py results_fractilesUHS/hazard_uhs-mean_1.csv False

'''    
###############################################################################
# set params
###############################################################################

hazcurvefile = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/complete_model/2023_final/results_uhs_fractiles_ta_pref/hazard_curve-mean-PGA_1.csv'   
hazcurvefile = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/complete_model/2023_final/results_uhs_fractiles/hazard_curve-mean-PGA_1.csv'  
hazcurvefile = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/complete_model/2023_final/results_uhs_fractiles_site_class_c/hazard_curve-mean-PGA_1.csv'  

siteDict, imls, investigation_time = return_annualised_haz_curves(hazcurvefile)

###############################################################################
# get IM
###############################################################################

im = path.split(hazcurvefile)[-1].split('-')[-1].split('_')[0]

siteclass = argv[1] # B-E

###################################################################################
# match city name to uhsDict
###################################################################################

# first parse city file
cwd = getcwd()
if cwd.startswith('/Users'): #mac
    citycsv = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/shared/nsha_cities.csv'
elif cwd.startswith('/nas'):
    citycsv = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/shared/nsha_cities.csv'
    
lines = open(citycsv).readlines()
    
# make city dict
cityDict = []
for line in lines:
    dat = line.strip().split(',')
    tmpdict = {'city':dat[2], 'lon':float(dat[0]), 'lat':float(dat[1])} 
    cityDict.append(tmpdict)


# now match cities
for j in range(0, len(siteDict)):
    siteDict[j]['place'] = 'null'
    
    for city in cityDict:
        if around(city['lon'], decimals=2) == around(siteDict[j]['lon'], decimals=2) \
           and around(city['lat'], decimals=2) == around(siteDict[j]['lat'], decimals=2):
           
           # add place
           siteDict[j]['place'] = city['city']
           print(city['city'])
           
###################################################################################
# write hazard curves
###################################################################################
outtxt = 'LOCALITY,LON,LAT,' + ','.join(['poe-'+str(x) for x in imls]) + '\n'
for j, sd in enumerate(siteDict):
    if not siteDict[j]['place'] == 'null':
       outtxt += ','.join((sd['place'], str('%0.2f' % sd['lon']), str('%0.2f' % sd['lat']))) + ',' \
                  + ','.join([str('%0.4e' % x) for x in sd['poe_probs_annual']]) + '\n'
    
outfile = 'nsha23_hazard_curves_SC_'+siteclass+'_'+im+'.csv'
f = open(outfile, 'w')
f.write(outtxt)
f.close()            
