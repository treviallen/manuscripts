import matplotlib.pyplot as plt
from numpy import arange, array, interp, unique
from os import path, getcwd, sep
from sys import argv
from gmt_tools import cpt2colormap 
from misc_tools import remove_last_cmap_colour, get_log_xy_locs
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes

import matplotlib as mpl
mpl.style.use('classic')

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)

'''
example run:
    run plt_uhs_curves.py results_fractilesUHS/hazard_uhs-mean_1.csv False
    /nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/source_models/complete_model/final/results_fractilesUHS/hazard_uhs-mean_1.csv

'''    
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

pltLog = False

# get colours
if getcwd().startswith('/nas'):
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/gay-flag-1978.cpt'
else:
    cptfile = '//Users//trev//Documents//DATA//GMT//cpt//gay-flag-1978.cpt'

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
    citycsv = '/Users/trev/Documents/Geoscience_Australia/NSHA2018/shared/nsha_cities.csv'
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
# plt 10 & 2% hazard curves
###################################################################################

fig = plt.figure(1, figsize=(10, 12))

if altPlaces == False:
    places = ['Perth', 'Darwin', 'Adelaide', 'Melbourne', 'Hobart', 'Canberra', 'Sydney', 'Brisbane']
else:
    places = ['Wongan Hills', 'Karratha', 'Yulara', 'Port Pirie', 'Hawker', 'Cooma', 'Leongatha', 'Morwell']
    #places = ['Wongan Hills', 'Darwin', 'Adelaide', 'Kimba', 'Hawker', 'Canberra', 'Sydney', 'Morwell']

if plt1170 == False:
    probidx = [0, 2]
else:
    probidx = [0]

for j, pi in enumerate(probidx):
    pltlon = []
    pltlat = []
    
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
                        
                pltlon.append(uhs['lon'])
                pltlat.append(uhs['lat'])
                    
    # plt AS1170.4 spectra
    if plt1170 == True:
        if pltLog == True:
            plt.plot(per1170[1:], shp1170[1:], 'k-', lw=2.0, label='AS1170.4 Be')
        else:
            plt.plot(per1170, shp1170, 'k-', lw=2.0, label='AS1170.4 Be')
        
    plt.xlabel('Period (s)', fontsize=18)
    if plt1170 == True:
        plt.ylabel('Normalised '+str('%0.0f' % (100*float(probabilities[pi])))+'% in 50-year AEP SA', fontsize=18)
    else:
        plt.ylabel(str('%0.0f' % (100*float(probabilities[pi])))+'% in 50-year AEP SA (g)', fontsize=18)
    
    if j == 0:
        if pltLog == True:
            plt.legend(loc=3, fontsize=14)
            plt.xlim([0.1, 4])
        else:
            plt.legend(loc=1, fontsize=14)
            plt.xlim([0.0, 4])
    else:
        plt.xlim([0.0, 4])
        
    
    ###############################################################################
    # make map inset
    ###############################################################################
    
    if j == 0:
        
        axins = inset_axes(ax,
                       width="35%",  # width = 30% of parent_bbox
                       height=2.5,  # height : 1 inch
                       loc=9)
        
        m = Basemap(projection='merc',\
                    llcrnrlon=111,llcrnrlat=-45, \
                    urcrnrlon=156,urcrnrlat=-9,\
                    rsphere=6371200.,resolution='l',area_thresh=10000)
                    
        m.drawmapboundary(fill_color='0.8', zorder=0)
        m.fillcontinents(color='w', lake_color='0.8', zorder=1)
        m.drawcoastlines()
        m.drawcountries()
        m.drawstates()
        
        # plt locs
        for i, place in enumerate(places):
            x, y = m(pltlon[i], pltlat[i])
            m.plot(x, y, 's', ms=8, markeredgecolor='k', markerfacecolor=tuple(cs[i]))
            
    

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