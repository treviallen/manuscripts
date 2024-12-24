#from matplotlib.mlab import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource
from numpy import arange, mean, percentile, array, unique, where, argsort, floor, ceil
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from os import path, walk, system
#from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser, CsvCatalogueWriter
from misc_tools import remove_last_cmap_colour, dictlist2array
from mapping_tools import drawshapepoly, labelpolygon, annotate_cities, get_map_polygons, get_field_data, make_street_map
from io_catalogues import parse_ga_event_query
import shapefile

def return_csg_data(jsonfile):
    import json

    with open(jsonfile) as f:
        data = json.load(f)

    csg_data = []
    for feature in data['features']:
        tmp = {'lon':feature['geometry']['coordinates'][0],
               'lat':feature['geometry']['coordinates'][1],
               'status':feature['properties']['CSGSTATUS']}
        
        # append to list
        csg_data.append(tmp)
        
    return csg_data

def parse_iris_stationlist(stationlist):
    lines = open(stationlist).readlines()[3:]
    
    staDict = []
    
    for line in lines:
        dat = line.strip().split('|')
        tmp = {'sta': dat[1], 'lat': float(dat[2]), 'lon': float(dat[3]), \
               'elev': float(dat[4]), 'place': dat[5]}
        staDict.append(tmp)
        
    return staDict  

plt.rcParams['pdf.fonttype'] = 42
mpl.style.use('classic')

##########################################################################################
doLocal = False
# local
clat = -32.34
clon = 150.87
if doLocal == False:
    urcrnrlat = clat + 0.13
    llcrnrlat = clat - 0.13
    urcrnrlon = clon + 0.16
    llcrnrlon = clon - 0.16

##########################################################################################
# set up street map
##########################################################################################

# set map centroid
clon = mean([llcrnrlon, urcrnrlon])
clat = mean([llcrnrlat, urcrnrlat])
            
degrng = urcrnrlon-llcrnrlon
ll_buffer = 0.2
'''
plt, m, ax = make_street_map(clat, clon, service='World_Imagery', ll_buffer = 0.15, \
             xpixels = 1500, plt_inset = False, plt_marker = False)
'''
fig = plt.figure(figsize=(18,10))
#plt.tick_params(labelsize=16)
ax = fig.add_subplot(111)

m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat, epsg=3857) #3112)
#m.drawparallels(arange(-90.,90.,0.5), labels=[1,0,0,0],fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.0)
#m.drawmeridians(arange(0.,360.,0.5), labels=[0,0,0,1], fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.0)
#http://server.arcgisonline.com/arcgis/rest/services

m.arcgisimage(service='World_Imagery', xpixels = 1500, verbose= True)

parSpace = 0.2
merSpace = 0.2
scaleLength = 10
res = 'f'

slon = m.llcrnrlon + 0.18*(m.urcrnrlon - m.llcrnrlon)
slat = m.llcrnrlat + 0.93*(m.urcrnrlat - m.llcrnrlat)
slon0 = mean([m.urcrnrlon, m.llcrnrlon])
slat0 = mean([m.urcrnrlat, m.llcrnrlat])

m.drawmapscale(slon, slat, slon0, slat0, scaleLength, fontsize=14, barstyle='fancy', zorder=10000, fontcolor='w')


'''
    Map Services:
        ESRI_Imagery_World_2D - pref!
        ESRI_StreetMap_World_2D
        NatGeo_World_Map
        NGS_Topo_US_2D
        Ocean_Basemap
        USA_Topo_Maps
        World_Imagery
        World_Physical_Map
        World_Shaded_Relief
        World_Street_Map
        World_Terrain_Base
        World_Topo_Map
        
ESRI_Imagery_World_2D - good sat image
NatGeo_World_Map - nice road map - hard to read
'''

##########################################################################################
# add cities
##########################################################################################
numCities = 9
annotate_cities(numCities, plt, m, markersize=10, markerfacecolor='maroon', markeredgecolor='w', \
                markeredgewidth=1.5, fs=16, weight='medium')
"""
# add extra locs
import matplotlib.patheffects as PathEffects
path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")]
txtoff = 0.1

x, y = m(133.07, -14.92)
plt.plot(x, y, 's', markerfacecolor='k', markeredgecolor='k', markeredgewidth=0.5, markersize=6, zorder=11000)
x, y = m(133.07+txtoff, -14.92+txtoff)
plt.text(x, y, 'Mataranka', size=14, ha='left', weight='normal', path_effects=path_effects)

x, y = m(133.37, -16.25)
plt.plot(x, y, 's', markerfacecolor='k', markeredgecolor='k', markeredgewidth=0.5, markersize=6, zorder=11000)
x, y = m(133.37-txtoff, -16.25+txtoff)
plt.text(x, y, 'Daly Waters', size=14, ha='right', weight='normal', path_effects=path_effects)

x, y = m(133.54, -17.55)
plt.plot(x, y, 's', markerfacecolor='k', markeredgecolor='k', markeredgewidth=0.5, markersize=6, zorder=11000)
x, y = m(133.54-txtoff, -17.55+txtoff)
plt.text(x, y, 'Elliott', size=14, ha='right', weight='normal', path_effects=path_effects)
"""
##########################################################################################
# add DYFI?
##########################################################################################
from gmt_tools import cpt2colormap, remove_last_cmap_colour

cptfile = '//Users//trev//Documents//DATA//GMT//cpt//mi_pop.cpt'
ncols = 10
cmap, zvals = cpt2colormap(cptfile, ncols+1, rev=False)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncols)))

import json
jsonFilePath = 'felt_reports/felt_reports_1km.geojson'

with open(jsonFilePath) as f:
    data = json.load(f)

dyfi_dict = []
for feature in data['features']:
    tmp = {'geomerty':feature['geometry']['coordinates'][0],
           'centroid':feature['properties']['center']['coordinates'],
           'intensity':feature['properties']['intensityFine'],
           'nresp':feature['properties']['nresp']}
    
    # append to list
    dyfi_dict.append(tmp)
    

nresp = 0
min_resp = 0
for dyfi in dyfi_dict:        
    # add to list greater than minObs
    if dyfi['nresp'] > min_resp:
        
        # now plot
        pltx = array(dyfi['geomerty'])[:,0]
        plty = array(dyfi['geomerty'])[:,1]
        
        x, y = m(pltx, plty)
        colidx = int(round(dyfi['intensity']))-1
        c= tuple(cs[colidx][:-1])
        plt.fill(x, y, fc=c, ec='none', lw=0.25, alpha=0.5)
        
    nresp += dyfi['nresp']
        
##########################################################################################
# plt stations
##########################################################################################

def return_all_au_station_data(au_station_file):
    from datetime import datetime
    from os import getcwd
    
    lines = open(au_station_file).readlines()
    
    sta_dict = []
    
    for line in lines:
        dat = line.strip().split('\t')
        #print(line)
        
        if int(dat[5]) < 1:
            dat[5] = 1
        if int(dat[7]) < 1:
            dat[7] = 1
            
        tmp = {'sta':dat[0], 'stlo':float(dat[1]), 'stla':float(dat[2]), 
               'startdate':datetime(int(dat[4]),int(dat[5]),1), 
               'enddate':datetime(int(dat[6]),int(dat[7]),1)}
        
        # append to sta_dict
        sta_dict.append(tmp)
        
    return sta_dict

rdk1_list = return_all_au_station_data('rdk1_stations.txt')
rdk_lon = dictlist2array(rdk1_list, 'stlo')
rdk_lat = dictlist2array(rdk1_list, 'stla')
rdk_sta = dictlist2array(rdk1_list, 'sta')

import matplotlib.patheffects as PathEffects
path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")]


# label stas
'''
urcrnrlat
llcrnrlat
urcrnrlon
llcrnrlon
'''
x,y = m(rdk_lon, rdk_lat)
plt.plot(x, y, '^', c='limegreen', ms=12, zorder=1000, label='Rapid Deployment Kits Phase 1')

for sta, slon, slat in zip(rdk_sta, rdk_lon, rdk_lat):
    if slon >= llcrnrlon-ll_buffer and slon <= urcrnrlon+ll_buffer \
        and slat >= llcrnrlat-ll_buffer and slat <= urcrnrlat+ll_buffer:
            print(sta)
            x,y = m(slon+0.004, slat+0.0015)
            plt.text(x, y, sta, size=12, c='k', va='bottom', ha='left', weight='normal', style='italic', \
                     path_effects=path_effects, zorder=11000)


rdk1_list = return_all_au_station_data('rdk2_stations.txt')
rdk_lon = dictlist2array(rdk1_list, 'stlo')
rdk_lat = dictlist2array(rdk1_list, 'stla')
rdk_sta = dictlist2array(rdk1_list, 'sta')

x,y = m(rdk_lon, rdk_lat)
plt.plot(x, y, '^', c='orange', ms=12, zorder=1000, label='Rapid Deployment Kits Phase 2')

for sta, slon, slat in zip(rdk_sta, rdk_lon, rdk_lat):
    if slon >= llcrnrlon-ll_buffer and slon <= urcrnrlon+ll_buffer \
        and slat >= llcrnrlat-ll_buffer and slat <= urcrnrlat+ll_buffer:
            print(sta)
            x,y = m(slon+0.004, slat+0.0015)
            plt.text(x, y, sta, size=12, c='k', va='bottom', ha='left', weight='normal', style='italic', \
                     path_effects=path_effects, zorder=11000)


rs_list = return_all_au_station_data('rs_station.txt')

x,y = m(rs_list[0]['stlo'], rs_list[0]['stla'])
plt.plot(x, y, '^', c='r', ms=12, zorder=1000, label='Raspberry Shake Network')

x,y = m(rs_list[0]['stlo']+0.004, rs_list[0]['stla']+0.0015)
plt.text(x, y, 'R7AF5', size=12, c='k', va='bottom', ha='left', weight='normal', style='italic', \
         path_effects=path_effects, zorder=11000)

##########################################################################################
# get land & lake polygons for masking
##########################################################################################
from io_catalogues import parse_ga_event_query
gadat = parse_ga_event_query('2024_mswl_earthquakes_export.csv')

la = dictlist2array(gadat, 'lat')
lo = dictlist2array(gadat, 'lon')
mag = dictlist2array(gadat, 'mag_ml')
x,y = m(lo, la)
#scatter = plt.scatter(x,y,s=-75+mag*80, c='yellow', linewidths=0.25, zorder=10000, alpha=0.7, label='Earthquake Epicentres')
scatter = plt.scatter(x,y,s=-75+mag*80, c='crimson', linewidths=0.25, zorder=10000, alpha=0.7, label='Earthquake Epicentres')

plt.legend(loc=4, numpoints=1, fontsize=12)


##########################################################################################
# make map inset
##########################################################################################

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

if doLocal == False:
    axins = zoomed_inset_axes(ax, 0.0018, loc=3)

m2 = Basemap(projection='merc',\
            llcrnrlon=111,llcrnrlat=-45, \
            urcrnrlon=156,urcrnrlat=-9,\
            rsphere=6371200.,resolution='l',area_thresh=10000)
            
m2.drawmapboundary(fill_color='0.8')
m2.fillcontinents(color='w', lake_color='0.8') #, zorder=0)
m2.drawcoastlines()
m2.drawcountries()
m2.drawstates()

# fill main area
xv = mean([llcrnrlon, urcrnrlon])
yv = mean([llcrnrlat, urcrnrlat])
x, y = m2(xv, yv)
plt.plot(x, y, 'rs',ms=6)

##########################################################################################
# label states
##########################################################################################
'''
state = ['WA', 'NT', 'SA', 'QLD', 'NSW', 'VIC', 'TAS']
slat = [-26, -21.0, -29.5, -23.0, -32.5, -37.1, -42.]
slon = [122, 133.5, 135.0, 144.5, 146.5, 143.6, 147.0]
for i, st in enumerate(state):
    x, y = m(slon[i], slat[i])
    plt.text(x, y, st, size=11, horizontalalignment='center', verticalalignment='center', weight='normal')
'''
##########################################################################################
# add colourbar
##########################################################################################

# set colourbar
plt.gcf().subplots_adjust(bottom=0.1)
cax = fig.add_axes([0.34,0.06,0.33,0.03]) # setup colorbar axes.

norm = mpl.colors.Normalize(vmin=0.5, vmax=10.5)#myb
cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal', alpha=0.5)

# set cb labels
ticks = range(1,11)
rom_num = ['I', 'II', 'III', 'IV', 'V', 'VI','VII','VIII','IX','X']
cb.set_ticks(ticks)
cb.set_ticklabels(rom_num)

titlestr = 'Macroseismic Intensity'
cb.set_label(titlestr, fontsize=16)

plt.savefig('muswellbrook_rdks.png', format='png', bbox_inches='tight', dpi=200)
#plt.savefig('camden_csg_boreholes.svg', format='svg', bbox_inches='tight', dpi=150)
plt.show()
