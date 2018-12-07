from numpy import array, arange, argsort, where, delete, hstack, sqrt, \
                  unique, mean, percentile, log10, ceil, floor, \
                  nan, isnan, around, diff, interp, exp, ones_like
from os import path, sep, mkdir, getcwd, walk
from shapely.geometry import Point, Polygon
#from osgeo import ogr
from datetime import datetime
from sys import argv
import shapefile
import matplotlib.pyplot as plt
from matplotlib import colors, colorbar, style
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
import matplotlib.gridspec as gridspec
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes

style.use('classic')

# import non-standard functions
try:
    from catalogue_tools import weichert_algorithm, aki_maximum_likelihood, bval2beta
    from oq_tools import get_oq_incrementalMFD, beta2bval #, bval2beta
    from mapping_tools import get_field_data, get_field_index, drawoneshapepoly, \
                              drawshapepoly, labelpolygon, get_WGS84_area
    #from catalogue.parsers import parse_ggcat
    from catalogue.writers import ggcat2ascii
    from mag_tools import nsha18_bilin_mw2ml
    from gmt_tools import cpt2colormap 
    from misc_tools import remove_last_cmap_colour, get_log_xy_locs
    
    #from misc_tools import listdir_extension
    #from make_nsha_oq_inputs import write_oq_sourcefile
except:
    cwd = getcwd().split(sep)
    pythonpath = sep.join(pt[0:-3])+sep+'tools'
    print '\nSet environmental variables, e.g.:\n\nexport PYTHONPATH='+pythonpath+':$PYTHONPATH\n'

        
###############################################################################
# set defaults
###############################################################################
fig = plt.figure(1, figsize=(10, 9))
#gs1 = gridspec.GridSpec(1, 2)
#gs1.update(wspace=0.025, hspace=0.05)

if getcwd().startswith('/nas'):
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/Paired_10.cpt'
    shpfile1  = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/GIS/AUS_adm/AUS_adm1.shp'

else:
    cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//Paired_10.cpt'
    
ncolours = 11
cmap, zvals = cpt2colormap(cptfile, ncolours)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncolours-1)))

###############################################################################
# parse shapefile and make shapely objects
###############################################################################

print 'Reading source shapefile...'
sf = shapefile.Reader(shpfile1)
shapes = sf.shapes()
polygons = []
polygonsCopy = []
for poly in shapes:
    polygons.append(Polygon(poly.points))
    polygonsCopy.append(Polygon(poly.points))
    
###############################################################################
# make map1
###############################################################################
# parse data
ax = fig.add_subplot(111)

# set national-scale basemap
llcrnrlat = -45
urcrnrlat = -8
llcrnrlon = 105
urcrnrlon = 155
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

m2 = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,
            projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
            resolution='l',area_thresh=1000.)

# annotate
m2.drawcoastlines(linewidth=0.5,color='0.25')
m2.drawcountries()
m2.drawstates()

# draw parallels and meridians.
ll_space = 6
m2.drawparallels(arange(-90.,90.,ll_space/2.0), labels=[0,0,0,0],fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.5)
m2.drawmeridians(arange(0.,360.,ll_space), labels=[0,0,0,0], fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.5)

# plt source zone boundary
# drawoneshapepoly(m, plt, sf, field, value, **kwargs)
drawoneshapepoly(m2, plt, sf, 'NAME_1', 'New South Wales', col='r', fillshape=True)

'''
xlim = ax.get_xlim()
xtxt = xlim[1] * 0.02
ylim = ax.get_ylim()
ytxt = ylim[1] * 0.02
plt.text(xtxt, ytxt, 'a)', fontsize=17, va='bottom', ha='left')
'''
# label polygons
#labelpolygon(m2, plt, sf, 'CODE')
  
# save figure
plt.savefig('nsw_fill.png', fmt='png', bbox_inches='tight')
    
plt.show()                
