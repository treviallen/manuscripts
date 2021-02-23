#from matplotlib.mlab import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from numpy import arange, mean, percentile, array, unique, where, argsort, floor , loadtxt
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from mapping_tools import drawshapepoly
from os import path, walk, system
from mapping_tools import labelCentroid
import shapefile

plt.rcParams['pdf.fonttype'] = 42
mpl.style.use('classic')

##########################################################################################
#108/152/-44/-8
urcrnrlat = -7.
llcrnrlat = -45.
urcrnrlon = 152.
llcrnrlon = 107.
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(13,8))
plt.tick_params(labelsize=14)
ax = fig.add_subplot(111)
'''
epsg = 3112
service='Ocean_Basemap'
#$service='World_Imagery'
#service='World_Shaded_Relief'
m = Basemap(projection='mill',llcrnrlon=llcrnrlon ,llcrnrlat=llcrnrlat,
            urcrnrlon=urcrnrlon ,urcrnrlat=urcrnrlat, resolution = 'l', epsg = epsg)
    
# xpixels controls the pixels in x direction, and if you leave ypixels
# None, it will choose ypixels based on the aspect ratio
xpixels = 1500
m.arcgisimage(service=service, xpixels = xpixels, verbose= True, dpi=300)
'''

m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='l',area_thresh=2000.)

# draw coastlines, state and country boundaries, edge of map.
m.shadedrelief()
#m.drawcoastlines()
m.drawstates(linewidth=0.5, color='0.3')

m.drawparallels(arange(-90.,90.,6.), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,10.), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
#m.drawmapscale(144, -34.8, 146., -38.5, 400, fontsize = 16, barstyle='fancy', zorder=100)

##########################################################################################
# add ficticious earthquakes
##########################################################################################

import matplotlib.patheffects as path_effects
data = loadtxt('scenario_lat_lons.csv', delimiter=',')
x, y = m(data[:,1], data[:,0])
plt.plot(x, y, '*', mfc='r', mec='k', mew=0.25, markersize=8, label='Scenario Localities')

'''
path_effect=[path_effects.withStroke(linewidth=3, foreground='w')]
x, y = m(data[:,0]+0.75, data[:,1]+0.5)
for i in range(0, len(x)):
    plt.text(x[i], y[i], 'E'+str(i+1), fontsize=14, c='b', \
             va='bottom', ha='left', style='italic', path_effects=path_effect)
'''    
plt.legend(loc=3, fontsize=14, numpoints=1)

# add neotectonic domains
shpfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/source_models/zones/shapefiles/Domains/Domains_Sep2011.shp'
sf = shapefile.Reader(shpfile)
drawshapepoly(m, plt, sf, col='k', lw=.5)

plt.savefig('scenario_locations.png', format='png', dpi=300, bbox_inches='tight')
#plt.savefig('neac_ml_zones.svg', format='svg', bbox_inches='tight')
plt.show()
