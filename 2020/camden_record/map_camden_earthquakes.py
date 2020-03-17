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
if doLocal == False:
    urcrnrlat = -31.6
    llcrnrlat = -36.2
    urcrnrlon = 154
    llcrnrlon = 148.5

##########################################################################################
# set up street map
##########################################################################################

# set map centroid
clon = mean([llcrnrlon, urcrnrlon])
clat = mean([llcrnrlat, urcrnrlat])
            
degrng = urcrnrlon-llcrnrlon
ll_buffer = (urcrnrlon - llcrnrlon) / 2.
plt, m, ax = make_street_map(clat, clon, service='ESRI_Imagery_World_2D', ll_buffer=ll_buffer, \
             xpixels = 1500, plt_inset = False, plt_marker = False, inset_loc=4, inset_multiplier=0.1)

'''
    Map Services:
        ESRI_Imagery_World_2D
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
numCities = 10
#annotate_cities(numCities, plt, m, marker='s')

##########################################################################################
# add earthquakes
##########################################################################################
'''
shpfile = '../NT/shapefiles/PetroleumTitles28August2019.shp'
shpfile = 'PetroleumTitles28August2019.shp'
sf = shapefile.Reader(shpfile)
drawshapepoly(m, plt, sf, col='r',lw=0.75, alpha=0.5, fillshape = True)

gadat = parse_ga_event_query('earthquakes_export_2012-16_250.edit.csv')

la = dictlist2array(gadat, 'lat')
lo = dictlist2array(gadat, 'lon')
mag = dictlist2array(gadat, 'mag_ml')
x,y = m(lo, la)
plt.scatter(x,y,s=mag*5, c='dodgerblue', label='2012-2015', zorder=10000)

legend1 = ax.legend(loc="upper left", title="Year Range")
ax.add_artist(legend1)
handles, labels = scatter.legend_elements(prop="sizes", alpha=0.6)
legend2 = ax.legend(handles, labels, loc="upper right", title="Sizes")
'''
##########################################################################################
# plt stations
##########################################################################################

stationlist = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Networks/AU/gmap-stations-edit.txt'
stationlist = '/Users/trev/Documents/Networks/AU/gmap-stations-edit.txt'
staDict = parse_iris_stationlist(stationlist)

# plot stations
sta_lon = dictlist2array(staDict, 'lon')
sta_lat = dictlist2array(staDict, 'lat')
sta_code = dictlist2array(staDict, 'sta')

cmdnet = ('MABG', 'YARR', 'CATI', 'OKDL', 'WTPK')
cmdlat = []
cmdlon = []
for stla, stlo, sta in zip(sta_lat, sta_lon, sta_code):
    for cn in cmdnet:
        if sta == cn:
            print(cn)
            cmdlat.append(stla)
            cmdlon.append(stlo)

x,y = m(sta_lon, sta_lat)
plt.plot(x, y, '^', c='w', ms=12, zorder=1000, label='ANSN')

x,y = m(cmdlon, cmdlat)
plt.plot(x, y, '^', c='yellow', ms=12, zorder=1000, label='Camden Seismic Network')

# label stas
urcrnrlat
llcrnrlat
urcrnrlon
llcrnrlon
'''
import matplotlib.patheffects as PathEffects
path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")]

for sta, slon, slat in zip(sta_code, sta_lon, sta_lat):
    if slon >= llcrnrlon-ll_buffer and slon <= urcrnrlon+ll_buffer \
        and slat >= llcrnrlat-ll_buffer and slat <= urcrnrlat+ll_buffer:
            print(sta)
            x,y = m(slon-0.005, slat+0.004)
            plt.text(x, y, sta, size=15, c='royalblue', va='bottom', ha='right', weight='normal', \
                     path_effects=path_effects, zorder=11000)
'''
plt.legend(loc=2, numpoints=1, fontsize=11)

##########################################################################################
# get land & lake polygons for masking
##########################################################################################
'''
polys = get_map_polygons(m)

#mask_outside_polygon(polys[1][::-1], ax=None)
#mask_outside_polygons(polys, 'lightskyblue', plt)

# get lake ploygons
polygons = []
for polygon in m.lakepolygons:
    poly = polygon.get_coords()
    plt.fill(poly[:,0], poly[:,1], 'lightblue')
    polygons.append(poly)
'''
##########################################################################################
# make map inset
##########################################################################################

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

if doLocal == False:
    axins = zoomed_inset_axes(ax, 0.03, loc=4)

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

plt.savefig('camden_earthquakes.png', format='png', bbox_inches='tight', dpi=150)
plt.show()
