from matplotlib.mlab import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource
from numpy import arange, mean, percentile, array, unique, where, argsort, floor, ceil, ones_like, sin, cos, radians
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

#GA projection
lon_0 = 134
lat_1 = -36
lat_2 = -18

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
#m.shadedrelief()
m.drawstates()
m.drawcoastlines(linewidth=0.5, linestyle='solid', color='k')
m.drawparallels(arange(-90.,90.,4.), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,6.), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)


##########################################################################################
# add shapefiles
##########################################################################################
#shpfile = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/source_models/zones/shapefiles/Leonard2008/source_model_leonard_2008.shp'
shpfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/source_models/zones/shapefiles/Domains/Domains_Sep2011.shp'
#shpfile = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/source_models/zones/shapefiles/ARUP_Background/ARUP_background_source_model.shp'
#shpfile = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/source_models/zones/shapefiles/NSHA13/NSHA13_regional_source_model_simplified.shp'
sfz = shapefile.Reader(shpfile)
drawshapepoly(m, plt, sfz, col='blue', lw=1.5)

# label shapes
#labelpolygon(m, plt, sf, 'NAME', fweight='normal', fsize=16, addOutline=True)

##########################################################################################
# add simple faults
##########################################################################################
nfsmshp = '/nas/active/ops/community_safety/ehp//georisk_earthquake//neotectonic//Seismicity_Scenario_models//Hazard Map working 2018//ARCGIS//FSM lines//FSD_simple_faults.shp'
    
sf = shapefile.Reader(nfsmshp)
shapes = sf.shapes()
for shape in shapes:
    lons = []
    lats = []
    for xy in shape.points:
        lons.append(xy[0])
        lats.append(xy[1])
    
    x, y = m(lons, lats)
    
    # set colour by 
    plt.plot(x, y, '-', c='g', lw=1.5)


##########################################################################################
# plt stress vectors
##########################################################################################

shpfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/source_models/zones/shapefiles/Other/SHMax_Rajabi_2016.shp'
sfv = shapefile.Reader(shpfile)

lat = get_field_data(sfv, 'LAT', 'float')
lon = get_field_data(sfv, 'LON', 'float')
shmax = get_field_data(sfv, 'SHMAX', 'float')

##########################################################################################
# annotate arrows
##########################################################################################
# test to make sure arrows in right spot
#x, y = m(lon, lat)
#m.plot(x, y, 'r.')

xm, ym = m(lon, lat)
#dx = ones_like(x) * 1000.
#dy = ones_like(y) * 1000.
for x, y, sh in zip(xm, ym, shmax):
    alen = 60000
    head_length = 30000
    dy = alen * cos(radians(sh)) + head_length * cos(radians(sh))
    dx = alen * sin(radians(sh)) + head_length * sin(radians(sh))
    ax.arrow(x-dx, y-dy, dx, dy, head_width=30000, head_length=head_length, lw = 0.3, fc='0.4', fill=True, length_includes_head = True)
    
    # plt opposite
    sh2 = 180. + sh
    dy = alen * cos(radians(sh2)) + head_length * cos(radians(sh2))
    dx = alen * sin(radians(sh2)) + head_length * sin(radians(sh2))
    ax.arrow(x-dx, y-dy, dx, dy, head_width=30000, head_length=head_length, lw=0.3, fc='0.4', fill=True, length_includes_head = True)

##########################################################################################
# get average per zone and plot
##########################################################################################

src_codes = get_field_data(sfz, 'CODE', 'str')
src_shapes = sfz.shapes()

av_shmax, sig_shmax = get_aus_shmax_vectors(src_codes, src_shapes)

# get centroid and plot
alen = 180000
head_length = 60000
   
for i, shape in enumerate(src_shapes):
    centroid = Polygon(shape.points).centroid.wkt
    centroid = centroid.strip('PIONT').replace('(',' ').replace(')',' ').split()
    cx, cy = m(float(centroid[0]),float(centroid[1]))
    
    dy = alen * cos(radians(av_shmax[i])) + head_length * cos(radians(av_shmax[i]))
    dx = alen * sin(radians(av_shmax[i])) + head_length * sin(radians(av_shmax[i]))
    ax.arrow(cx-dx, cy-dy, dx, dy, head_width=60000, head_length=head_length, lw = 2., ec='maroon', fc='maroon', fill=True, length_includes_head = True)
    
    dy = alen * cos(radians(av_shmax[i]+180)) + head_length * cos(radians(av_shmax[i]+180))
    dx = alen * sin(radians(av_shmax[i]+180)) + head_length * sin(radians(av_shmax[i]+180))
    ax.arrow(cx-dx, cy-dy, dx, dy, head_width=60000, head_length=head_length, lw = 2., ec='maroon', fc='maroon', fill=True, length_includes_head = True)



plt.savefig('map_stress_vectors_domains.png', format='png', bbox_inches='tight', dpi=200, transparent=True)
plt.show()
