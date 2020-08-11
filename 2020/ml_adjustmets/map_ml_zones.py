#from matplotlib.mlab import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from numpy import arange, mean, percentile, array, unique, where, argsort, floor , loadtxt
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from os import path, walk, system
from mapping_tools import labelCentroid

plt.rcParams['pdf.fonttype'] = 42
mpl.style.use('classic')

##########################################################################################
#108/152/-44/-8
urcrnrlat = -5.
llcrnrlat = -49.
urcrnrlon = 157.
llcrnrlon = 98.
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
#m.drawstates()

m.drawparallels(arange(-90.,90.,6.), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,10.), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
#m.drawmapscale(144, -34.8, 146., -38.5, 400, fontsize = 16, barstyle='fancy', zorder=100)
'''
##########################################################################################
# plot gebco
##########################################################################################

print 'Reading netCDF file...'
nc = NetCDFFile('//Users//tallen//Documents//DATA//GMT//GEBCO//Australia_30c.nc')

zscale =20. #gray
zscale =50. #colour
data = nc.variables['elevation'][:] / zscale
lons = nc.variables['lon'][:]
lats = nc.variables['lat'][:]

# transform to metres
nx = int((m.xmax-m.xmin)/500.)+1
ny = int((m.ymax-m.ymin)/500.)+1

topodat = m.transform_scalar(data,lons,lats,nx,ny)

print 'Getting colormap...'
# get colormap
#cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//mby_topo-bath_mod.cpt'
cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//mby_topo-bath.cpt'
#cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//gray.cpt'
cmap, zvals = cpt2colormap(cptfile, 256)
cmap = remove_last_cmap_colour(cmap)

# make shading
print 'Making map...'
ls = LightSource(azdeg = 180, altdeg = 45)
norm = mpl.colors.Normalize(vmin=-8000/zscale, vmax=5000/zscale)#myb
rgb = ls.shade(topodat, cmap=cmap, norm=norm)
im = m.imshow(rgb)
'''

##########################################################################################
# parse mag zones
##########################################################################################

import shapefile
from shapely.geometry import Point, Polygon
from mapping_tools import get_field_data, drawshapepoly, labelpolygon

shpfile = 'shapefile/australia_ml_regions.shp'
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()

drawshapepoly(m, plt, sf, fillcolor='none', edgecolor='r', lw=2.5, zorder=0)
labelpolygon(m, plt, sf, 'ML_REGION', fsize=24, addOutline=True) 

##########################################################################################
# add ficticious earthquakes
##########################################################################################

import matplotlib.patheffects as path_effects
data = loadtxt('ficticious_eq_locs.csv', delimiter=',')
x, y = m(data[:,0], data[:,1])
plt.plot(x, y, '*', mfc='k', mec='k', mew=0, markersize=16, label='Fictitious Epicentres')

path_effect=[path_effects.withStroke(linewidth=3, foreground='w')]
x, y = m(data[:,0]+0.75, data[:,1]+0.5)
for i in range(0, len(x)):
    plt.text(x[i], y[i], 'E'+str(i+1), fontsize=14, c='b', \
             va='bottom', ha='left', style='italic', path_effects=path_effect)
    
plt.legend(loc=3, fontsize=14, numpoints=1)

plt.savefig('neac_ml_zones.png', format='png', bbox_inches='tight')
plt.savefig('neac_ml_zones.svg', format='svg', bbox_inches='tight')
plt.show()
