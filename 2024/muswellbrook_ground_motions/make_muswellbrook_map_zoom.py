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
    urcrnrlat = clat + 0.12
    llcrnrlat = clat - 0.12
    urcrnrlon = clon + 0.12
    llcrnrlon = clon - 0.12

##########################################################################################
# set up street map
##########################################################################################

# set map centroid
clon = mean([llcrnrlon, urcrnrlon])
clat = mean([llcrnrlat, urcrnrlat])
            
degrng = urcrnrlon-llcrnrlon
ll_buffer = 0.2
plt, m, ax = make_street_map(clat, clon, service='World_Imagery', ll_buffer = 0.15, \
             xpixels = 1500, plt_inset = False, plt_marker = False)

parSpace = 0.2
merSpace = 0.2
scaleLength = 10
res = 'f'

slon = m.llcrnrlon + 0.08*(m.urcrnrlon - m.llcrnrlon)
slat = m.llcrnrlat + 0.85*(m.urcrnrlat - m.llcrnrlat)
slon0 = mean([m.urcrnrlon, m.llcrnrlon])
slat0 = mean([m.urcrnrlat, m.llcrnrlat])

m.drawmapscale(slon, slat, slon0, slat0, scaleLength, fontsize=12, barstyle='simple', zorder=10000, fontcolor='w')


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
# add shapefiles
##########################################################################################
'''
shpfile = '../NT/shapefiles/PetroleumTitles28August2019.shp'
shpfile = 'PetroleumTitles28August2019.shp'
sf = shapefile.Reader(shpfile)
drawshapepoly(m, plt, sf, col='r',lw=0.75, alpha=0.5, fillshape = True)

jsonfile='NSW_CSG_Boreholes.geojson'
csg_data = return_csg_data(jsonfile)
status = dictlist2array(csg_data, 'status')
csg_lon = dictlist2array(csg_data, 'lon')
csg_lat = dictlist2array(csg_data, 'lat')

idx1 = where(status == 'Permanently Sealed')
x,y = m(csg_lon[idx1], csg_lat[idx1])
plt.plot(x,y,'ro',ms=7,label='Permanently Sealed')

idx1 = where(status == 'Not Producing Gas')
x,y = m(csg_lon[idx1], csg_lat[idx1])
plt.plot(x,y,'o',c='orange',ms=7,label='Not Producing Gas')

idx1 = where(status == 'Producing Gas')
x,y = m(csg_lon[idx1], csg_lat[idx1])
plt.plot(x,y,'o',c='limegreen',ms=7,label='Producing Gas', zorder=10000)
'''
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
scatter = plt.scatter(x,y,s=-75+mag*80, c='yellow', linewidths=0.25, zorder=10000, alpha=0.7, label='Earthquake Epicentres')

plt.legend(loc=4, numpoints=1, fontsize=11)


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

plt.savefig('muswellbrook_rdks.png', format='png', bbox_inches='tight', dpi=150)
#plt.savefig('camden_csg_boreholes.svg', format='svg', bbox_inches='tight', dpi=150)
plt.show()
