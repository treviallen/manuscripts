#matplotlib quiver example

from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource
from numpy import arange, mean, percentile, array, unique, where, argsort, floor, sqrt
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from os import path, walk, system
from shapely.geometry import Polygon, Point
import pickle
#from obspy.imaging.beachball import Beach
#from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser, CsvCatalogueWriter
from misc_tools import remove_last_cmap_colour, listdir_extension
from mapping_tools import drawshapepoly, get_field_data, drawoneshapepoly, distance, reckon
import shapefile

mpl.style.use('classic')

##########################################################################################
#108/152/-44/-8
urcrnrlat = 1.0
llcrnrlat = -28.
urcrnrlon = 155.
llcrnrlon = 112
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(18,10))
plt.tick_params(labelsize=16)
ax = fig.add_subplot(111)

m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='l',area_thresh=1000.)

# draw coastlines, state and country boundaries, edge of map.
#m.shadedrelief()
#m.bluemarble()
#m.etopo()
m.drawcoastlines()
m.drawstates()
m.drawcountries()
m.fillcontinents(color='w',lake_color='w')
m.drawmapboundary(fill_color='w')
m.drawparallels(arange(-90.,90.,4.), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,6.), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
#m.drawmapscale(144, -34.8, 146., -38.5, 400, fontsize = 16, barstyle='fancy', zorder=100)

##########################################################################################
# plot gebco
##########################################################################################

print('Reading netCDF file...')
nc = NetCDFFile('//Users//trev//Documents//DATA//GMT//GEBCO//Australia_30c.nc')

zscale =20. #gray
zscale =30. #colour
data = nc.variables['z'][:] / zscale
lons = nc.variables['lon'][:]
lats = nc.variables['lat'][:]

# transform to metres
nx = int((m.xmax-m.xmin)/500.)+1
ny = int((m.ymax-m.ymin)/500.)+1

topodat = m.transform_scalar(data,lons,lats,nx,ny)

print('Getting colormap...')
# get colormap
cptfile = '//Users//trev//Documents//DATA//GMT//cpt//wiki-2.0.cpt'
cmap, zvals = cpt2colormap(cptfile, 30)
cmap = remove_last_cmap_colour(cmap)

# make shading
print('Making map...')
ls = LightSource(azdeg = 180, altdeg = 5)
#norm = mpl.colors.Normalize(vmin=-8000/zscale, vmax=5000/zscale)
norm = mpl.colors.Normalize(vmin=-1000/zscale, vmax=1900/zscale)#wiki
rgb = ls.shade(topodat.data, cmap=cmap, norm=norm)
im = m.imshow(rgb, alpha=1.)

##########################################################################################
# add domains
##########################################################################################
import matplotlib.patheffects as PathEffects
path_effects=[PathEffects.withStroke(linewidth=2.5, foreground="w")]

cptfile = '//Users//trev//Documents//DATA//GMT//cpt//keshet.cpt' # qual-dark-06.cpt keshet.cpt #aurora.cpt #cosam.cpt #Set1_06.cpt
ncols = 5
cmap, zvals = cpt2colormap(cptfile, ncols+1)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncols)))

shpfile = 'shapefiles/adj_neotectonic_domains.shp'
sf = shapefile.Reader(shpfile)

shapes = sf.shapes()
trts = get_field_data(sf, 'trt', 'str')
utrts = unique(trts)

labels = ['Precambrian Craton', 'Reactivated Proterozoic Crust', 'Extended Continental Crust', \
          'Phanerozoic Accretionary Crust', 'Passive Margin']

for c, utrt, label in zip(cs, utrts, labels):
    labeLegend = True
    for shape, trt in zip(shapes, trts):
        if trt == utrt:
            x = []
            y = []
            
            p = 0
            parts = shape.parts
            parts.append(len(shape.points))
            
            for prt in parts[1:]:
                #print('part',prt
                while p < prt:
                    x.append(shape.points[p][0])
                    y.append(shape.points[p][1])
                    p += 1
                
                # close polygon if not polyline
                if x[0] != x[-1] or y[0] != y[-1]:
                    x.append(x[0])
                    y.append(y[0])
            
                # plot each polygon                    
                xx, yy = m(x,y)
                #print(edgecolor)
                if labeLegend == True:
                    plt.fill(xx,yy, facecolor=c, edgecolor='none', linewidth=0.75, alpha=0.3, label=label)
                    labeLegend = False
                else:
                    plt.fill(xx,yy, facecolor=c, edgecolor='none', linewidth=0.75, alpha=0.3) 

##########################################################################################
# plot faults
##########################################################################################

def plt_faults(m, plt, ffile, c='k', lw=2):
    lines = open(ffile).readlines()[1:]
    lats = []
    lons = []
    for line in lines:
        dat = line.strip().split()
        lons.append(float(dat[0]))
        lats.append(float(dat[1]))
        
    xx, yy = m(array(lons),array(lats))
    plt.plot(xx, yy, '-', c=c, lw=lw)
'''
ffile = 'fault_geometries/BANDA_model/block_boundaries/TIMO'
plt_faults(m, plt, ffile, c='k', lw=2)
ffile = 'fault_geometries/BANDA_model/block_boundaries/SERA'
plt_faults(m, plt, ffile, c='r', lw=2)
'''
#parse shapefile
shpfile = '/Users/trev/Documents/DATA/GIS/Faults/PB2002_boundaries.shp'
sf = shapefile.Reader(shpfile)

faults2plot = ['SU/AU', 'TI-AU', 'BH\BS', 'BS-AU', 'BH-BS']

discLen = 150. # km
triLen = 40.
halfTriLen = triLen / 2.
for fault in faults2plot:
    flolas = drawoneshapepoly(m, plt, sf, 'Name', fault, lw=1.5, polyline=True)
    
    discLon = []
    discLat = []
    discAzm = []
    remainder = discLen #+ discLen/2.
    for poly in flolas:
        for i in range(1, len(poly[0])):
            # get dist from point 1 to point 2
            rng, az, baz = distance(poly[1][i-1], poly[0][i-1], poly[1][i], poly[0][i])
            
            lens = arange(discLen-remainder, rng, discLen)
            
            # get xy locs for lens
            if len(lens) > 0:
                for l in lens:
                    discPos = reckon(poly[1][i-1], poly[0][i-1], l, az)
                    discLon.append(discPos[0])
                    discLat.append(discPos[1])
                    discAzm.append(az)
                    
                remainder = rng - lens[-1]
            else:
                remainder += rng

    # now, at each point, make dip triangle
    for dlo, dla, daz in zip(discLon, discLat, discAzm):
        # pt1
        triLons = [reckon(dla, dlo, halfTriLen, daz)[0]]
        triLats = [reckon(dla, dlo, halfTriLen, daz)[1]]
        
        # pt2
        triLons.append(reckon(dla, dlo, halfTriLen, daz+180)[0])
        triLats.append(reckon(dla, dlo, halfTriLen, daz+180)[1])
        
        # pt3 - assum equilateral triange
        triHeight = sqrt(triLen**2 - halfTriLen**2)
        triLons.append(reckon(dla, dlo, triHeight, daz-90)[0])
        triLats.append(reckon(dla, dlo, triHeight, daz-90)[1])
    
        xx, yy = m(triLons, triLats)
        #plt.plot(xx, yy, 'o', c='r', lw=2, )
        plt.fill(xx,yy, facecolor='k', edgecolor='none', linewidth=0.75, alpha=1)
    
    
'''
# label polygons
x, y = m(146, -40.1)
plt.text(x, y, 'Bass\nBasin', size=12, c='r', ha='center', va='top', weight='normal', style='italic', path_effects=path_effects)

x, y = m(148.5, -39.)
plt.text(x, y, 'Gippsland\nBasin', size=12, c='r', ha='center', weight='normal', style='italic', path_effects=path_effects)
'''
plt.legend(loc=1, fontsize=13)

plt.savefig('figures/tectonic_setting.png',fmt='png',dpi=150,bbox_inches='tight')
plt.show()