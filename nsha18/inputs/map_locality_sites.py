from matplotlib.mlab import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource
from numpy import arange, mean, percentile, array, unique, where, argsort, floor, ceil, ones_like, sin, cos, radians, loadtxt
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from os import path, walk, system
from misc_tools import remove_last_cmap_colour
from mapping_tools import drawshapepoly, labelpolygon, get_field_data
from tools.source_shapefile_builder import get_aus_shmax_vectors
from shapely.geometry import Polygon
import shapefile

plt.rcParams['pdf.fonttype'] = 42
mpl.style.use('classic')

##########################################################################################
# parse epicentres
##########################################################################################


##########################################################################################
#108/152/-44/-8
urcrnrlat = -8.
llcrnrlat = -46.
urcrnrlon = 157.
llcrnrlon = 109.

llcrnrlat = -45
urcrnrlat = -5.5
llcrnrlon = 106.
urcrnrlon = 153.
                
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(20,11))
ax = fig.add_subplot(111)
plt.tick_params(labelsize=8)

m = Basemap(projection='lcc',\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            lat_1=lat_1,lat_2=lat_2,lon_0=lon_0, \
            rsphere=6371200.,resolution='l',area_thresh=2000.)

#lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
# draw coastlines, state and country boundaries, edge of map.
m.shadedrelief()
m.drawstates()
#m.drawcountries()
m.drawparallels(arange(-90.,90.,4.), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,6.), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)


##########################################################################################
# parse sites
##########################################################################################
csvfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/shared/nsha_localities.csv'
#shpfile = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/source_models/zones/shapefiles/ARUP_Background/ARUP_background_source_model.shp'
locs = loadtxt(csvfile, delimiter=',')

##########################################################################################
# annotate points
##########################################################################################
# test to make sure arrows in right spot
x, y = m(locs[:,0], locs[:,1])
m.plot(x, y, 'ro',markersize=6, label='AS1170.4 Localities')
plt.legend(loc=3, numpoints=1)


plt.savefig('map_locality_sites.png', format='png', bbox_inches='tight', dpi=200)
plt.show()
