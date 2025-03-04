from os import path, getcwd, sep
from numpy import array, unique, nan
from misc_tools import checkfloat

###############################################################################
# parse hazard file
###############################################################################

hazardfile = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/complete_model/2023_final/results_uhs_fractiles_SC_B/quantile_map-0.85_1.csv'

lines = open(hazardfile).readlines()
headers = [x for x in lines[1].strip().split(',')]

# get keys from uhs file
keys = lines[1].strip().split(',')[2:]

# get peridos in keys
periods = []
tmpProb = []
for key in keys:
    tmpProb.append(key.split('-')[1])
    
    if key.startswith('PGA'):
        periods.append(0.0)
    else:
        periods.append(float(key.split('(')[1].split(')')[0]))

# get unique probabilities
probabilities = unique(tmpProb)[::-1] # reorder
periods = unique(periods)

# site site data
hazDict = []
for line in lines[2:]:
    dat = [float(x) for x in line.strip().split(',')]
    tmpdict = {'lon':float(dat[0]), 'lat':float(dat[1])}
    
    for i, key in enumerate(keys):
        tmpdict[key] = dat[i+2]
    
    '''
    for i, prob in enumerate(probabilities):
        startidx = i * len(periods) + 2
        stopidx = startidx + len(periods)
        siteUHS = [float(x) for x in dat[startidx:stopidx]]
        
        tmpdict[prob] = array(siteUHS)
    '''
    hazDict.append(tmpdict)
    
"""
###############################################################################
# parse TA preferred hazard file
###############################################################################

#hazardfile = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/complete_model/2023_final/results_uhs_fractiles_ta_pref_site_class_b/hazard_map-mean_1.csv'

lines = open(hazardfile).readlines()
headers = [x for x in lines[1].strip().split(',')]

# site site data
hazDictTA = []
for line in lines[2:]:
    dat = [float(x) for x in line.strip().split(',')]
    tmpdict = {'lon':float(dat[0]), 'lat':float(dat[1])}
    
    for i, key in enumerate(keys):
        tmpdict[key] = dat[i+2]
    
    hazDictTA.append(tmpdict)

# now match locs
for j in range(0, len(hazDict)):
    for hta in hazDictTA:
        if hta['lon'] == hazDict[j]['lon'] \
           and hta['lat'] == hazDict[j]['lat']:
           
           # add place
           hazDict[j]['PGA-0.020-TA_PREF'] = hta['PGA-0.02']
           hazDict[j]['PGA-0.033-TA_PREF'] = hta['PGA-0.033']
"""
###################################################################################
# match city name to hazDict
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
for j in range(0, len(hazDict)):
    hazDict[j]['place'] = 'null'
    
    for city in cityDict:
        if city['lon'] == hazDict[j]['lon'] \
           and city['lat'] == hazDict[j]['lat']:
           
           # add place
           hazDict[j]['place'] = city['city']

###################################################################################
# load AS1170.4 table
###################################################################################

as1170table = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/shared/nsha18_AS1170.4_cmp.csv'

as1170dict = []
lines = open(as1170table).readlines()
for line in lines[1:]:
    dat = line.strip().split(',')
    tmpdict = {'city':dat[0], 'lon':float(dat[1]), 'lat':checkfloat(dat[2]), 'as1170':checkfloat(dat[3]), \
    	         'nshm12_500':checkfloat(dat[4]), 'nsha18_475':float(dat[5]), 'nsha18_500':float(dat[6])} 
    as1170dict.append(tmpdict)

# now match cities
for j in range(0, len(hazDict)):
    hazDict[j]['as1170'] = nan
    hazDict[j]['nsha18_475'] = nan
    for as1170 in as1170dict:
        if as1170['lon'] == hazDict[j]['lon'] \
           and as1170['lat'] == hazDict[j]['lat']:
           
           # add place
           hazDict[j]['as1170'] = as1170['as1170']
           hazDict[j]['nsha18_475'] = as1170['nsha18_475']
           
           
###################################################################################
# write table
###################################################################################

tabtxt = 'LON,LAT,PLACE,AS1170.4-2007,AS1170.4-2018,NSHA18 10% PGA,NSHA23 85th 10% PGA,NSHA23 85th 3.3% NPGA,NSHA23 85th 2% NPGA\n'

for hd in hazDict:
    if not hd['place'].endswith(' max'):
        if not hd['place'] == 'null':
            tabtxt += ','.join((str(hd['lon']), str(hd['lat']), hd['place'], str('%0.2f' % hd['as1170']), \
                               str('%0.2f' % max([hd['as1170'],0.08])), str('%0.3f' % hd['nsha18_475']), \
                               str('%0.3f' % hd['PGA-0.1']), str('%0.3f' % (0.6667 * hd['PGA-0.033'])), \
                               str('%0.3f' % (0.6667 * hd['PGA-0.02'])))) + '\n'
                       
# write file
f = open('draft_nsha23_1170.4_pga_85th.csv', 'w')
f.write(tabtxt)
f.close()