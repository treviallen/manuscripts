from os import path, getcwd, sep
from numpy import array, unique, nan, nanmean, nanmedian, nanstd
from misc_tools import checkfloat

###############################################################################
# parse hazard file
###############################################################################

def parse_haz_dict(hazardfile):

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
        
        hazDict.append(tmpdict)
        
    return hazDict

###################################################################################
# set files and parse
###################################################################################

meanhazardfile = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/complete_model/2023_final/results_uhs_fractiles_SC_B/hazard_map-mean_1.csv'
medianhazardfile = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/complete_model/2023_final/results_uhs_fractiles_SC_B/quantile_map-0.5_1.csv'
fractilefile16 = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/complete_model/2023_final/results_uhs_fractiles_SC_B/quantile_map-0.16_1.csv'
fractilefile84 = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/complete_model/2023_final/results_uhs_fractiles_SC_B/quantile_map-0.84_1.csv'

meanhazDict = parse_haz_dict(meanhazardfile)
medianhazDict = parse_haz_dict(medianhazardfile)
fractileDict16 = parse_haz_dict(fractilefile16)
fractileDict84 = parse_haz_dict(fractilefile84)
    
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
for j in range(0, len(meanhazDict)):
    meanhazDict[j]['place'] = 'null'
    
    for city in cityDict:
        if city['lon'] == meanhazDict[j]['lon'] \
           and city['lat'] == meanhazDict[j]['lat']:
           
           # add place
           meanhazDict[j]['place'] = city['city']

           
###################################################################################
# write 10% table
###################################################################################

tabtxt = 'PLACE,LON,LAT,MEAN,16th PERCENTILE,50th PERCENTILE,84th PERCENTILE \n'

for i, hd in enumerate(meanhazDict):
    if not hd['place'].endswith(' max'):
      if not hd['place'] == 'null':
        tabtxt += ','.join((hd['place'], str(hd['lon']), str(hd['lat']), \
                           str('%0.3f' % hd['PGA-0.1']), \
                           str('%0.3f' % fractileDict16[i]['PGA-0.1']), \
                           str('%0.3f' % medianhazDict[i]['PGA-0.1']), \
                           str('%0.3f' % fractileDict84[i]['PGA-0.1']))) + '\n'
                       
# write file
f = open('nsha23_10pc_fractile.csv', 'w')
f.write(tabtxt)
f.close()

###################################################################################
# write 2% table
###################################################################################

tabtxt = 'PLACE,LON,LAT,MEAN,16th PERCENTILE,50th PERCENTILE,84th PERCENTILE \n'

for i, hd in enumerate(meanhazDict):
    if not hd['place'].endswith(' max'):
      if not hd['place'] == 'null':
        tabtxt += ','.join((hd['place'], str(hd['lon']), str(hd['lat']), \
                           str('%0.3f' % hd['PGA-0.02']), \
                           str('%0.3f' % fractileDict16[i]['PGA-0.02']), \
                           str('%0.3f' % medianhazDict[i]['PGA-0.02']), \
                           str('%0.3f' % fractileDict84[i]['PGA-0.02']))) + '\n'
                       
# write file
f = open('nsha23_2pc_fractile.csv', 'w')
f.write(tabtxt)
f.close()