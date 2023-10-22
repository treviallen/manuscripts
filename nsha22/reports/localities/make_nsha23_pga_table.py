from os import path, getcwd, sep
from numpy import array, unique, nan, nanmean, nanmedian, nanstd
from misc_tools import checkfloat

###############################################################################
# parse hazard file
###############################################################################

hazardfile = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/complete_model/2023_final/results_uhs_fractiles_SC_B/hazard_map-mean_1.csv'

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
    hazDict[j]['nshm12_500'] = nan
    for as1170 in as1170dict:
        if as1170['lon'] == hazDict[j]['lon'] \
           and as1170['lat'] == hazDict[j]['lat']:
           
           # add place
           hazDict[j]['as1170'] = as1170['as1170']
           hazDict[j]['nshm12_500'] = as1170['nshm12_500']
           
###################################################################################
# load NSHA18 table
###################################################################################

nsha18table = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/shared/nsha18_haz_curves_PGA.csv'

nsha18dict = []
lines = open(nsha18table).readlines()
for line in lines[2:]:
    dat = line.strip().split(',')
    tmpdict = {'city':dat[0], 'lon':float(dat[1]), 'lat':checkfloat(dat[2]), \
    	         'nsha18_475':float(dat[5]), 'nsha18_2475':float(dat[11])} 
    nsha18dict.append(tmpdict)

# now match cities
meandiff10 = []
meanpcchange10 = []
meandiff2 = []
meanpcchange2 = []
for j in range(0, len(hazDict)):
    hazDict[j]['nsha18_475'] = nan
    hazDict[j]['nsha18_2475'] = nan
    hazDict[j]['pga_10_percent_change'] = nan
    hazDict[j]['pga_2_percent_change'] = nan
    hazDict[j]['pga_10_difference'] = nan
    hazDict[j]['pga_2_difference'] = nan
    
    for nsha18 in nsha18dict:
        if nsha18['lon'] == hazDict[j]['lon'] \
           and nsha18['lat'] == hazDict[j]['lat']:
           
           # add place
           hazDict[j]['nsha18_475'] = nsha18['nsha18_475']
           hazDict[j]['nsha18_2475'] = nsha18['nsha18_2475']
           
           # calc % change
           hazDict[j]['pga_10_percent_change'] = 100 * (hazDict[j]['PGA-0.1'] - nsha18['nsha18_475']) / nsha18['nsha18_475']
           hazDict[j]['pga_2_percent_change']  = 100 * (hazDict[j]['PGA-0.02'] - nsha18['nsha18_2475']) / nsha18['nsha18_2475']
           
           hazDict[j]['pga_10_difference'] = hazDict[j]['PGA-0.1'] - nsha18['nsha18_475']
           hazDict[j]['pga_2_difference']  = hazDict[j]['PGA-0.02'] - nsha18['nsha18_2475']
           
           meanpcchange10.append(hazDict[j]['pga_10_percent_change'])
           meanpcchange2.append(hazDict[j]['pga_2_percent_change'])
           
           meandiff10.append(hazDict[j]['pga_10_difference'])
           meandiff2.append(hazDict[j]['pga_2_difference'])

print(nanmean(array(meanpcchange10)))
print(nanstd(array(meanpcchange10)))
print(nanmedian(array(meanpcchange2)))
print('\n')
print(nanmean(array(meandiff10)))
print(nanstd(array(meandiff10)))
print(nanmedian(array(meandiff2)))

           
###################################################################################
# write 10% table
###################################################################################

tabtxt = 'PLACE,LON,LAT,AS1170.4-2007,AS1170.4-2018,NSHM13_0.10,NSHA18_PGA_0.10,PGA_0.10,2018 TO 2023 PERCENT CHANGE,ABS DIFFERNCE \n'

for hd in hazDict:
    if not hd['place'].endswith(' max'):
      if not hd['place'] == 'null':
        tabtxt += ','.join((hd['place'], str(hd['lon']), str(hd['lat']), \
                           str('%0.2f' % hd['as1170']), str('%0.2f' % max([hd['as1170'],0.08])), \
                           str('%0.2f' % hd['nshm12_500']), str('%0.3f' % hd['nsha18_475']), str('%0.3f' % hd['PGA-0.1']), \
                           str('%0.1f' % hd['pga_10_percent_change']), str('%0.3f' % hd['pga_10_difference']))) + '\n'
                       
# write file
f = open('nsha23_10pc_pga.csv', 'w')
f.write(tabtxt)
f.close()

###################################################################################
# write 2% table
###################################################################################

tabtxt = 'PLACE,LON,LAT,AS1170.4-2007,NSHA18_PGA_0.02,NSHA23_PGA_0.02,2018 TO 2023 PERCENT CHANGE,ABS DIFFERNCE\n'

for hd in hazDict:
    if not hd['place'].endswith(' max'):
      if not hd['place'] == 'null':
        tabtxt += ','.join((hd['place'], str(hd['lon']), str(hd['lat']), \
                           str('%0.2f' % (1.8 * hd['as1170'])), \
                           str('%0.3f' % hd['nsha18_2475']), str('%0.3f' % hd['PGA-0.02']), \
                           str('%0.1f' % hd['pga_2_percent_change']), str('%0.3f' % hd['pga_2_difference']))) + '\n'
                       
# write file
f = open('nsha23_2pc_pga.csv', 'w')
f.write(tabtxt)
f.close()