import matplotlib.pyplot as plt
from numpy import arange, array, interp, unique
from os import path, getcwd, sep
from sys import argv
from gmt_tools import cpt2colormap 
from misc_tools import remove_last_cmap_colour, get_log_xy_locs

import matplotlib as mpl
mpl.style.use('classic')


###############################################################################
# set params
###############################################################################

uhsfile = argv[1]

plt1170 = argv[2]

if plt1170 == 'True':
    plt1170 = True
else:
    plt1170 = False

'''
if pltProb == '10':
    probidx = 0
elif pltProb == '9':
    probidx = 1
elif pltProb == '2':
    probidx = 2
'''

altPlaces = False

pltLog = True

# get colours
#cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//Paired_08.cpt'
cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//gay-flag-1978.cpt'
ncolours = 9
cmap, zvals = cpt2colormap(cptfile, ncolours)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncolours-1)))


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
# plt 10 & 2% hazard curves
###################################################################################

fig = plt.figure(1, figsize=(11, 11))

if altPlaces == False:
    places = ['Perth', 'Darwin', 'Adelaide', 'Melbourne', 'Hobart', 'Canberra', 'Sydney', 'Brisbane']
else:
    places = ['Wongan Hills', 'Kalgoorlie', 'Port Pirie', 'Cooma', 'Yulara', 'Hawker', 'Leongatha', 'Morwell']

if plt1170 == False:
    probidx = [0, 2]
else:
    probidx = [0]

for j, pi in enumerate(probidx):
    ax = plt.subplot(2,1,j+1)
    # match city for plotting
    for i, place in enumerate(places):
        for uhs in uhsDict:
            if place == uhs['place']:
                if pltLog == True:
                    if plt1170 == True:
                        normSA = uhs[probabilities[pi]] / uhs[probabilities[pi]][0]
                        plt.loglog(periods[1:], normSA[1:], lw=2.0, c=cs[i], label=place)
                    else:
                        plt.loglog(periods[1:], uhs[probabilities[pi]][1:], lw=2.0, c=cs[i], label=place)
                else:
                    if plt1170 == True:
                        normSA = uhs[probabilities[pi]] / uhs[probabilities[pi]][0]
                        plt.plot(periods[0:], normSA, lw=2.0, c=cs[i], label=place)
                    else:
                        plt.plot(periods[0:], uhs[probabilities[pi]][0:], lw=2.0, c=cs[i], label=place)
                    
    # plt AS1170.4 spectra
    if plt1170 == True:
        if pltLog == True:
            plt.plot(per1170[1:], shp1170[1:], 'k-', lw=2.0, label='AS1170.4 Be')
        else:
            plt.plot(per1170, shp1170, 'k-', lw=2.0, label='AS1170.4 Be')
        
    plt.xlabel('Period (s)', fontsize=15)
    if plt1170 == True:
        plt.ylabel('Normalised '+str('%0.0f' % (100*float(probabilities[pi])))+'% in 50-year AEP Sa (g)', fontsize=15)
    else:
        plt.ylabel(str('%0.0f' % (100*float(probabilities[pi])))+'% in 50-year AEP Sa (g)', fontsize=15)
    
    if j == 0:
        if pltLog == True:
            plt.legend(loc=3, fontsize=13)
            plt.xlim([0.1, 5])
        else:
            plt.legend(loc=1, fontsize=13)
        


ylims = array(ax.get_ylim())
xlims = array(ax.get_xlim())

#plt.text(2.6, 85, 'a)', fontsize=15, va='top', ha='left')

###############################################################################
# save figs
###############################################################################
if altPlaces == False:
    caps = 'capitals'
else:
    caps = 'regions'
    
if pltLog == True:
    axtype = 'log'
else:
    axtype = 'lin'

if plt1170 == True:
    pngname = '_'.join(('uhs',caps,axtype,'norm.png'))

    plt.savefig(pngname, fmt='png', bbox_inches='tight')
    
else:
    pngname = '_'.join(('uhs',caps,axtype+'.png'))

    plt.savefig(pngname, fmt='png', bbox_inches='tight')

plt.show()