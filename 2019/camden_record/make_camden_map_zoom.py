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
# parse epicentres
##########################################################################################

'''
# parse HMTK csv
evdict = parse_ga_event_query('nt_earthquakes_export.csv')
year = dictlist2array(evdict, 'year')
month = dictlist2array(evdict, 'month')
mag = dictlist2array(evdict, 'mag')
lat = dictlist2array(evdict, 'lat')
lon = dictlist2array(evdict, 'lon')
'''
##########################################################################################
doLocal = False
# local
if doLocal == False:
    urcrnrlat = -33.85+0.14
    llcrnrlat = -34.36+0.14
    urcrnrlon = 151.
    llcrnrlon = 150.4

"""
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(18,10))
plt.tick_params(labelsize=16)
ax = fig.add_subplot(111)

m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='i',area_thresh=10.)
m.arcgisimage(service='ESRI_Imagery_World_2D')


# draw coastlines, state and country boundaries, edge of map.
#m.shadedrelief()
m.drawcoastlines(color='0.4')
m.drawstates()
m.drawcountries()
if doLocal == False:
    m.drawparallels(arange(-90.,90.,.2), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
    m.drawmeridians(arange(0.,360.,.2), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
else:
    m.drawparallels(arange(-90.,90.,1.), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
    m.drawmeridians(arange(0.,360.,2.), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
#m.drawmapscale(144, -34.8, 146., -38.5, 400, fontsize = 16, barstyle='fancy', zorder=100)


##########################################################################################
# plot gebco
##########################################################################################

print 'Reading netCDF file...'
try:
    nc = NetCDFFile('//Users//tallen//Documents//DATA//GMT//GEBCO//Australia_30c.nc')
except:
    nc = NetCDFFile('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/GEBCO/au_gebco.nc')

zscale =20. #gray
zscale =20. #colour
data = nc.variables['elevation'][:] / zscale
lons = nc.variables['lon'][:]
lats = nc.variables['lat'][:]

# transform to metres
nx = int((m.xmax-m.xmin)/500.)+1
ny = int((m.ymax-m.ymin)/500.)+1

topodat = m.transform_scalar(data,lons,lats,nx,ny)

print 'Getting colormap...'
# get colormap
try:
    cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//wiki-2.0.cpt'
    cmap, zvals = cpt2colormap(cptfile, 30)
except:
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/wiki-2.0.cpt'
    cmap, zvals = cpt2colormap(cptfile, 30)
    
cmap = remove_last_cmap_colour(cmap)
        
# make shading
print 'Making map...'
ls = LightSource(azdeg = 180, altdeg = 5)
norm = mpl.colors.Normalize(vmin=-8000/zscale, vmax=5000/zscale)#myb
norm = mpl.colors.Normalize(vmin=-1000/zscale, vmax=1900/zscale)#wiki

rgb = ls.shade(topodat, cmap=cmap, norm=norm)
im = m.imshow(rgb)
"""
##########################################################################################
# set up street map
##########################################################################################

# set map centroid
clon = mean([llcrnrlon, urcrnrlon])
clat = mean([llcrnrlat, urcrnrlat])
            
degrng = urcrnrlon-llcrnrlon
plt, m, ax = make_street_map(clat, clon, service='ESRI_Imagery_World_2D', ll_buffer = 0.4, \
             xpixels = 1500, plt_inset = False, plt_marker = False)

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
# add epicentres
##########################################################################################
"""
ncols = 12

# get year range
#minyear = 10*floor(year[-1]/10.)
#maxyear = 10*ceil(year[0]/10.)
minyear = 1960
maxyear = 2020
yearrng = float(round(maxyear - minyear))
#cmap = plt.cm.get_cmap('Spectral', ncols)

try:
    cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/temperature.cpt'
    cmap, zvals = cpt2colormap(cptfile, ncols+1, rev=False)
except:
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/temperature.cpt'
    cmap, zvals = cpt2colormap(cptfile, ncols+1, rev=False)
    
cmap = remove_last_cmap_colour(cmap)

cs = (cmap(arange(ncols)))
#cs = vstack((cs[0:1], cs[2:], [1., 1., 1., 1.]))

# get zorder for plotting
'''
sortidx = argsort(argsort(cat.data['magnitude']))
for i in range(0, len(cat.data['magnitude'])): #[0:100])):
    if cat.data['magnitude'][i] >= 2.5:
        #get colour idx
        year = cat.data['year'][i]
        colidx = int(round((ncols-1) * (year-minyear) / yearrng))
        x, y = m(cat.data['longitude'][i], cat.data['latitude'][i])
        zo = sortidx[i] + 20
        plt.plot(x, y, 'o', markerfacecolor=cs[colidx], markeredgecolor='k', markeredgewidth=0.5, \
                 markersize=-4 + cat.data['magnitude'][i]*3.2, zorder=zo, alpha=0.8)
'''
sortidx = argsort(argsort(mag))
for i in range(0, len(mag)): #[0:100])):
    if mag[i] >= 1.0:
        #get colour idx
        evyear = year[i] + month[i]/12.
        colidx = int(floor((ncols) * (evyear-minyear) / yearrng))
        x, y = m(lon[i], lat[i])
        zo = sortidx[i] + 20
        plt.plot(x, y, 'o', mfc=list(cs[colidx][0:3]), markeredgecolor='k', markeredgewidth=0.5, \
                 markersize=-4 + mag[i]*3.2, zorder=zo, alpha=0.8)
    
# make legend
legmag = [3., 5., 7.]
legh = []
for lm in legmag:
    x, y = m(0, 0)
    h = m.plot(x, y, 'ko', mfc='k', markersize=(-4 + lm*3.2), alpha=1., zorder=len(mag)+1, lw=2)
    legh.append(h[0])

mlstr = '$\mathregular{M_L}$'
l = plt.legend(legh, (mlstr+' 3.0', mlstr+' 5.0', mlstr+' 7.0'), loc=2, numpoints=1)
l.set_zorder(len(mag)+5)
"""
##########################################################################################
# add cities
##########################################################################################
numCities = 8
annotate_cities(numCities, plt, m, marker='s')

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

x,y = m(sta_lon, sta_lat)
plt.plot(x, y, '^', c='yellow', ms=12, zorder=12000, label='Camden Seismic Network')

# label stas
urcrnrlat
llcrnrlat
urcrnrlon
llcrnrlon

import matplotlib.patheffects as PathEffects
path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")]

for sta, slon, slat in zip(sta_code, sta_lon, sta_lat):
    if slon >= llcrnrlon and slon <= urcrnrlon \
        and slat >= llcrnrlat and slat <= urcrnrlat:
            print(sta)
            x,y = m(slon-0.005, slat+0.004)
            plt.text(x, y, sta, size=15, c='royalblue', va='bottom', ha='right', weight='normal', \
                     path_effects=path_effects, zorder=11000)

plt.legend(loc=2, numpoints=1, fontsize=12)

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
    axins = zoomed_inset_axes(ax, 0.0045, loc=3)

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
"""
ticks = arange(0, ncols/2.+1)
years = float(minyear) + ticks*10
labels = [str('%0.0f' % x) for x in years]

# normalise
norm = mpl.colors.Normalize(vmin=minyear, vmax=2020)

if doLocal == False:
    cax = fig.add_axes([0.7,0.3,0.02,0.4]) # setup colorbar axes.
else:
    cax = fig.add_axes([0.72,0.3,0.02,0.4]) # setup colorbar axes.
    
cb = colorbar.ColorbarBase(cax, cmap=cmap, orientation='vertical', alpha=0.8, norm=norm) # 
cb.set_ticks(years)
cb.set_ticklabels(labels)
cb.set_label('Year of Earthquake', rotation=270, labelpad=20, fontsize=15)
"""
plt.savefig('camden_csg_boreholes.png', format='png', bbox_inches='tight', dpi=150)
plt.show()
