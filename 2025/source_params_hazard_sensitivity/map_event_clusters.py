#matplotlib quiver example

from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource, ListedColormap
from numpy import arange, mean, percentile, array, unique, where, argsort, floor, sqrt, delete, argsort, log10, unique
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap, remove_last_cmap_colour, remove_first_cmap_colour
from os import path, walk, system
from shapely.geometry import Polygon, Point
import pickle
#from obspy.imaging.beachball import Beach
#from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser, CsvCatalogueWriter
from misc_tools import remove_last_cmap_colour, listdir_extension, get_mpl2_colourlist
from mapping_tools import drawshapepoly, get_field_data, drawoneshapepoly, distance, reckon, map_fault_dip_dirn
import shapefile

mpl.style.use('classic')

##########################################################################################
#108/152/-44/-8
urcrnrlat = -1.0
urcrnrlat = -6.0
llcrnrlat = -45.
urcrnrlon = 153.
llcrnrlon = 106
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(18,10))
plt.tick_params(labelsize=13)
ax = fig.add_subplot(111)

m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='l',area_thresh=500.) #i

# draw coastlines, state and country boundaries, edge of map.
#m.shadedrelief()
#m.bluemarble()
#m.etopo()
m.drawmapboundary(fill_color='0.8')
m.fillcontinents(color='w', lake_color='0.8') #, zorder=0)
m.drawcoastlines()
m.drawstates()
m.drawcountries()
m.drawparallels(arange(-90.,90.,4.), labels=[1,0,0,0],fontsize=13, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,6.), labels=[0,0,0,1], fontsize=13, dashes=[2, 2], color='0.5', linewidth=0.75)
#m.drawmapscale(144, -34.8, 146., -38.5, 400, fontsize = 16, barstyle='fancy', zorder=100)

##########################################################################################
# plot gebco
##########################################################################################
'''
print( 'Reading netCDF file...')
try:
    nc = NetCDFFile('//Users//trev//Documents//DATA//GMT//GEBCO//au_indo_gebco_2020.nc')
except:
    nc = NetCDFFile('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/GEBCO/au_gebco.nc')

zscale =30. #gray
#zscale =30. #colour
data = nc.variables['elevation'][:] / zscale
lons = nc.variables['lon'][:]
lats = nc.variables['lat'][:]

# transform to metres
nx = int((m.xmax-m.xmin)/500.)+1
ny = int((m.ymax-m.ymin)/500.)+1

topodat = m.transform_scalar(data,lons,lats,nx,ny)

print( 'Getting colormap...')
# get colormap
try:
    cptfile = '//Users//trev//Documents//DATA//GMT//cpt//wiki-2.0.cpt'
    cptfile = '//Users//trev//Documents//DATA//GMT//cpt//lightgrey.cpt'
    cptfile = '//Users//trev//Documents//DATA//GMT//cpt//grey_fade_2.cpt'
except:
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/wiki-2.0.cpt'

#cmap, zvals = cpt2colormap(cptfile, 30) # wiki
cmap, zvals = cpt2colormap(cptfile, 31) # grey
#cmap = remove_last_cmap_colour(cmap)

#cmap = ListedColormap
        
# make shading
print( 'Making map...')
ls = LightSource(azdeg = 135, altdeg = 10)
#norm = mpl.colors.Normalize(vmin=-8000/zscale, vmax=5000/zscale)#myb
#norm = mpl.colors.Normalize(vmin=-1000/zscale, vmax=1900/zscale)#wiki
norm = mpl.colors.Normalize(vmin=-2000/zscale, vmax=3500/zscale)#wiki

rgb = ls.shade(topodat.data, cmap=cmap, norm=norm)
im = m.imshow(rgb, alpha=1.0)
'''
##########################################################################################
# add domains
##########################################################################################
'''
import matplotlib.patheffects as PathEffects
path_effects=[PathEffects.withStroke(linewidth=2.5, foreground="w")]

cptfile = '//Users//trev//Documents//DATA//GMT//cpt//keshet.cpt' # qual-dark-06.cpt keshet.cpt #aurora.cpt #cosam.cpt #Set1_06.cpt
ncols = 7
cmap, zvals = cpt2colormap(cptfile, ncols+1)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncols)))
cs = delete(cs, 2, 0) # get rid of yellow

shpfile = '../../2019/nac_attenuation/shapefiles/adj_neotectonic_domains.shp'
sf = shapefile.Reader(shpfile)

shapes = sf.shapes()
trts = get_field_data(sf, 'trt', 'str')
utrts = unique(trts)

labels = ['Precambrian Craton', 'Reactivated Proterozoic Crust', 'Extended Continental Crust', \
          'Phanerozoic Accretionary Crust', 'Passive Margin', 'Sprigg Orogen']

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
                    plt.fill(xx,yy, facecolor=c, edgecolor='none', linewidth=0.75, alpha=0.2, label=label)
                    labeLegend = False
                else:
                    plt.fill(xx,yy, facecolor=c, edgecolor='none', linewidth=0.75, alpha=0.2) 

l1 = plt.legend(fontsize=11, loc=3, bbox_to_anchor=(0.125, 0))
'''
##########################################################################################
# plot faults
##########################################################################################


##########################################################################################
# annotate text
##########################################################################################
import matplotlib.patheffects as PathEffects
path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")]
"""
# label polygons
x, y = m(135, -18)
plt.text(x, y, 'North Australian\nCraton', size=11, c='maroon', ha='center', va='top', weight='normal', style='italic', path_effects=path_effects)

x, y = m(126, -14.5)
plt.text(x, y, 'Kimberly\nCraton', size=11, c='maroon', ha='center', va='top', weight='normal', style='italic', path_effects=path_effects)

x, y = m(118.75, -21.5)
plt.text(x, y, 'Pilbara\nCraton', size=11, c='maroon', ha='center', va='center', weight='normal', style='italic', path_effects=path_effects)

x, y = m(134.75, -32)
plt.text(x, y, 'Gawler\nCraton', size=11, c='maroon', ha='center', va='center', weight='normal', style='italic', path_effects=path_effects)

x, y = m(119.5, -29.5)
plt.text(x, y, 'Yilgarn\nCraton', size=11, c='maroon', ha='center', va='center', weight='normal', style='italic', path_effects=path_effects)

x, y = m(126, -12.5)
plt.text(x, y, 'Western\nExtended Crust', size=11, c='g', ha='center', weight='normal', style='italic', path_effects=path_effects)

x, y = m(148, -39.5)
plt.text(x, y, 'Eastern\nExtended Crust', size=11, c='g', ha='center', va='center', weight='normal', style='italic', path_effects=path_effects)

x, y = m(130, -24.5)
plt.text(x, y, 'Central Australia Reactivated\nProterozoic Crust', size=11, c='orangered', \
         va='center', ha='center', weight='normal', style='italic', path_effects=path_effects)

x, y = m(146.5, -28)
plt.text(x, y, 'Tasmanides', size=11, c='b', ha='center', weight='normal', style='italic', path_effects=path_effects)

x, y = m(129, -5)
plt.text(x, y, 'Banda Sea', size=11, c='k', ha='center', va='center', weight='normal', style='italic', path_effects=path_effects)

# add faults

x, y = m(129, -9.8)
plt.text(x, y, 'Timor Trough', size=10, c='k', va='center', ha='center', weight='normal', \
         rotation=10 , path_effects=path_effects)
         
x, y = m(117, -11.9)
plt.text(x, y, 'Sunda Arc', size=10, c='k', va='center', ha='center', weight='normal', \
         rotation=-0.5 , path_effects=path_effects)
         
'''x, y = m(122, -7)
plt.text(x, y, 'Flores-Wetar\nThrust', size=10, c='k', va='center', ha='center', weight='normal', \
         rotation=-0 , path_effects=path_effects)
'''         
x, y = m(141, -4.)
plt.text(x, y, 'New Guinea\nHighlands', size=10, c='k', va='center', ha='center', weight='normal', \
         rotation=-0 , path_effects=path_effects)
"""         


#plt.legend(loc=3, fontsize=12)
##########################################################################################
# set colours
##########################################################################################

# get stress range
logstressbins = arange(-.7, 1.7, 0.2)
logstressrng = logstressbins[-1] - logstressbins[0]

ncols = len(logstressbins) - 1

#logstressrng = float(round(maxstress - minstress))

cptfile = '/Users/trev/Documents/DATA/GMT/cpt/temperature.cpt'
cptfile = '//Users//trev//Documents//DATA//GMT//cpt//keshet.cpt'
#cptfile = '//Users//trev//Documents//DATA//GMT//cpt//plasma.cpt'

cmap, zvals = cpt2colormap(cptfile, ncols+1, rev=True)
cmap = remove_first_cmap_colour(cmap)

cs = (cmap(arange(ncols)))
'''
cmap = plt.get_cmap('viridis_r', len(logstressbins))
cs = (cmap(arange(len(logstressbins))))
'''
##########################################################################################
# add earthquakes
##########################################################################################

#csvfile = '../../2023/au_stress_drop/brune_stats.csv'
csvfile = 'brune_stats_cluster.csv'
#data = loadtxt(csvfile, delimiter=',', skiprows=1)

print(csvfile)
lines = open(csvfile).readlines()[1:]
lats = []
lons = []
mags = []
qual = []
stressdrops = []
cluster = []

for line in lines:
    dat = line.strip().split(',')
    lons.append(float(dat[2]))
    lats.append(float(dat[3]))
    mags.append(float(dat[8]))
    qual.append(float(dat[-2]))
    stressdrops.append(float(dat[10]))
    cluster.append(int(float(dat[-1])))

logstress = log10(array(stressdrops))
cluster = array(cluster)
lons = array(lons)
lats = array(lats)

unique_clusters = unique(cluster)

# map clusters
cols = get_mpl2_colourlist()
for i, uc in enumerate(unique_clusters):
    idx = where(cluster == uc)[0]
    x, y = m(lons[idx], lats[idx])
    m.plot(x, y, 'o', ms=8, c=cols[i], label='Cluster '+str(i+1))

plt.legend(loc=3, numpoints=1, fontsize=12	)

#########################################################################################
# finish
##########################################################################################

plt.savefig('mapped_event_clusters.png',fmt='png',dpi=300,bbox_inches='tight')
#plt.savefig('figures/fig_1.eps',fmt='pdf',dpi=300,bbox_inches='tight')
plt.show()