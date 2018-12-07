from matplotlib.mlab import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource
from numpy import arange, mean, percentile, array, unique, where, argsort, floor, ceil
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from os import path, walk, system
#from obspy.imaging.beachball import Beach
from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser, CsvCatalogueWriter
from misc_tools import remove_last_cmap_colour
from mapping_tools import drawshapepoly, labelpolygon
import shapefile

plt.rcParams['pdf.fonttype'] = 42
mpl.style.use('classic')

##########################################################################################
# parse epicentres
##########################################################################################

# parse HMTK csv
'''
hmtk_csv = '/nas/gemd/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/catalogue/data/NSHA18CAT_V0.1_hmtk_declustered.csv'
parser = CsvCatalogueParser(hmtk_csv)    
declcat = parser.read_file()
'''
hmtk_csv = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/catalogue/data/NSHA18CAT_V0.1_hmtk_declustered.csv'
parser = CsvCatalogueParser(hmtk_csv)    
cat = parser.read_file()

##########################################################################################
#108/152/-44/-8
urcrnrlat = -8.
llcrnrlat = -46.
urcrnrlon = 157.
llcrnrlon = 109.
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(12,12))
ax = fig.add_subplot(111)
plt.tick_params(labelsize=12)

m = Basemap(projection='merc',\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='l',area_thresh=2000.)

#lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
# draw coastlines, state and country boundaries, edge of map.
m.shadedrelief()
m.drawstates()
m.drawparallels(arange(-90.,90.,4.), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,6.), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)

##########################################################################################
# plt earthquakes
##########################################################################################

lats = cat.data['latitude']
lons = cat.data['longitude']
mags = cat.data['magnitude']
year = cat.data['year']

idx = where((mags >= 3.0) & (year >= 1977))[0]
x, y = m(lons[idx], lats[idx])
m.plot(x, y, 'o', ms=2.5, mfc='0.4', mec='0.4', zorder=1, label='Mw '+r'$\geq$'+' 3.0 (since 1977)')

idx = where((mags >= 5.6) & (year >= 1880))[0]
x, y = m(lons[idx], lats[idx])
m.plot(x, y, 'o', ms=9, mfc='maroon', mec='w', mew=1, zorder=2, label='Mw '+r'$\geq$'+' 5.6 (since 1880)')
plt.legend(loc=3, fontsize=14, numpoints=1)

##########################################################################################
# add shapefiles
##########################################################################################
shpfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/NSHM_18/Catalogue/Completeness_NSHA18/Leonard_17_mcomp_d.shp'
sf = shapefile.Reader(shpfile)
drawshapepoly(m, plt, sf,col='r',lw=2)

# label shapes
labelpolygon(m, plt, sf, 'NAME', fweight='normal', fsize=20, addOutline=True)

plt.savefig('map_mcomp_model.png', format='png', bbox_inches='tight') #, dpi=200)
plt.show()
