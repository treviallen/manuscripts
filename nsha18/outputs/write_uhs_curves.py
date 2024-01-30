import matplotlib.pyplot as plt
from numpy import arange, array, interp, unique
from os import path, getcwd, sep
from sys import argv
from gmt_tools import cpt2colormap 
from misc_tools import remove_last_cmap_colour, get_log_xy_locs

import matplotlib as mpl
mpl.style.use('classic')

site_class = argv[1].upper()

###############################################################################
# set params
###############################################################################

if getcwd().startswith('/nas'):
    uhsfile = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/complete_model/2023_final/results_uhs_fractiles_SC_'+site_class+'/hazard_uhs-mean_1.csv'

'''
if plt1170 == 'True':
    plt1170 = True
else:
    plt1170 = False
    
    
# get colours
if getcwd().startswith('/nas'):
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/gay-flag-1978.cpt'
else:
    cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//gay-flag-1978.cpt'
ncolours = 9
cmap, zvals = cpt2colormap(cptfile, ncolours)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncolours-1)))


'''

plt1170 = False

altPlaces = True

pltLog = False

###############################################################################
# parse uhs file
###############################################################################

lines = open(uhsfile).readlines()
headers = [x for x in lines[1].strip().split(',')]

# get keys from uhs file
keys = lines[1].strip().split(',')[2:]

# get peridos in keys
periods = []
tmpProb = []
for key in keys:
    tmpProb.append(key.split('~')[0])
    
    if key.startswith('0.1'):
        if key.endswith('PGA'):
            periods.append(0.0)
        else:
            periods.append(float(key.split('(')[-1][:-1]))

# get unique probabilities
probabilities = unique(tmpProb)[::-1] # reorder


# site site data
uhsDict = []
for line in lines[2:]:
    dat = [float(x) for x in line.strip().split(',')]
    tmpdict = {'lon':float(dat[0]), 'lat':float(dat[1])}
    
    for i, prob in enumerate(probabilities):
        startidx = i * len(periods) + 2
        stopidx = startidx + len(periods)
        siteUHS = [float(x) for x in dat[startidx:stopidx]]
        
        
        tmpdict[prob] = array(siteUHS)

    uhsDict.append(tmpdict)

###################################################################################
# match city name to uhsDict
###################################################################################

# first parse city file
cwd = getcwd()
if cwd.startswith('/Users'): #mac
    citycsv = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/shared/nsha_cities.csv'
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
for j in range(0, len(uhsDict)):
    for city in cityDict:
        if city['lon'] == uhsDict[j]['lon'] \
           and city['lat'] == uhsDict[j]['lat']:
           
           # add place
           uhsDict[j]['place'] = city['city']

###################################################################################
# set AS1170.4 shape for site class Be
###################################################################################

per1170 = arange(0, 4.1, 0.1)
shp1170 = []
for t in per1170:
    if t <= 0.1:
        shp1170.append(1.0 + 19.4*t)
    elif t > 0.1 and t <= 1.5:
        shp1170.append(min(0.88/t, 2.94))
    else:
        shp1170.append(1.32 / t**2)
        
shp1170 = array(shp1170)
           
###################################################################################
# write to csvs
###################################################################################
for i, prob in enumerate(probabilities):
    uhstxt = 'GEOSCIENCE AUSTRALIA NSHA18 UHS ACCELERATION VALUES IN UNITS OF G\nPLACE,LON,LAT,' + ','.join(['SA'+str(t) for t in periods]) + '\n'
    
    # loop through places
    for ud in uhsDict:
        uhstxt += ','.join((ud['place'], str('%0.2f' % ud['lon']), str('%0.2f' % ud['lat']), \
                            ','.join([str(h) for h in ud[prob]]))) + '\n'
                            
    # write to file
    uhsOutFile = 'nsha18_uhs_'+prob+'.csv'
    f = open(uhsOutFile, 'wb')
    f.write(uhstxt)
    f.close()
    
    
    
    
    
    
    
    
    
    
        
        
        
        
        
        
        
        
        
        
        
    
    