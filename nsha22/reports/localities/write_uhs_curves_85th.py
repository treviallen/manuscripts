from numpy import arange, array, interp, unique
from os import path, getcwd, sep
from sys import argv

'''
example run:
    run plt_uhs_curves.py results_fractilesUHS/hazard_uhs-mean_1.csv False

    use: results_uhs_fractiles_SC_B/quantile_uhs-0.85_1.csv

'''    

site_class = argv[1].upper()


###############################################################################
# set params
###############################################################################

uhsfile = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/complete_model/2023_final/results_uhs_fractiles_SC_'+site_class+'/quantile_uhs-0.85_1.csv'   

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
for j in range(0, len(uhsDict)):
    uhsDict[j]['place'] = 'null'
    
    for city in cityDict:
        if city['lon'] == uhsDict[j]['lon'] \
           and city['lat'] == uhsDict[j]['lat']:
           
           # add place
           uhsDict[j]['place'] = city['city']
           
###################################################################################
# plt 10 & 2% hazard curves
###################################################################################

for j, pi in enumerate(probabilities):
    
    outtxt = 'LOCALITY,' + ','.join(['SA('+str(x)+')' for x in periods]) + '\n'
    
    # match city for plotting
    for uhs in uhsDict:
        
        '''
        normSA = uhs[pi] / uhs[pi][0]
        
        outtxt += uhs['place'] + ',' + ','.join([str('%0.4f' % x) for x in normSA]) + '\n'
        '''
        outtxt += uhs['place'] + ',' + ','.join([str('%0.4f' % x) for x in uhs[pi]]) + '\n'
    
    outfile = path.join('uhs_curves', 'nsha23_uhs_'+pi[0:5]+'_SC_'+site_class+'_85th.csv')            
    f = open(outfile, 'w')
    f.write(outtxt)
    f.close()            
	