#from matplotlib.mlab import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from numpy import arange, mean, percentile, array, unique, where, argsort, floor , loadtxt
from netCDF4 import Dataset as NetCDFFile
#from gmt_tools import cpt2colormap
from os import path, walk, system
from mapping_tools import labelCentroid

plt.rcParams['pdf.fonttype'] = 42
mpl.style.use('classic')

##########################################################################################
bbox = '107.0/153.0/-43.0/-8.0'
urcrnrlat = -7.
llcrnrlat = -44.
urcrnrlon = 154.
llcrnrlon = 105.
lon_0 = mean([llcrnrlon, urcrnrlon])
#lon_0 = 134.
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(13,8))
plt.tick_params(labelsize=12)
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
#m.drawstates()

m.drawparallels(arange(-90.,90.,6.), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,10.), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
#m.drawmapscale(144, -34.8, 146., -38.5, 400, fontsize = 16, barstyle='fancy', zorder=100)

##########################################################################################
# add ficticious earthquakes
##########################################################################################

import matplotlib.patheffects as path_effects
csvfile = 'shakemap_event_list.csv'
#data = loadtxt(csvfile, delimiter=',', skiprows=1)

lines = open(csvfile).readlines()[1:]
lats = []
lons = []
mags = []
for line in lines:
    dat = line.strip().split(',')
    lons.append(float(dat[1]))
    lats.append(float(dat[2]))
    mags.append(float(dat[4]))
    
for lo, la, ma in zip(lons, lats, mags):
    x, y = m(lo, la)
    m.plot(x, y, marker='*', mfc='r', mec='k', mew=0.25, markersize=(5. * ma - 12), alpha=1., zorder=len(mags)+1)
    
# make legend
legmag = [4.5, 5.5, 6.5]
legh = []
for lm in legmag:
    x, y = m(0, 0)
    h = m.plot(x, y, '*', mfc='r', mec='k', mew=0.25, markersize=(5 * lm - 12), alpha=1., zorder=len(mags)+1)
    legh.append(h[0])

l = plt.legend(legh, ('MW 4.5', 'MW 5.5', 'MW 6.5'), loc=3, numpoints=1, fontsize=10, title="Magnitude")
l.set_zorder(len(mags)+5)

plt.savefig('atlas_event_map.png', format='png', bbox_inches='tight')
#plt.savefig('neac_ml_zones.eps', format='eps', bbox_inches='tight')
plt.show()
