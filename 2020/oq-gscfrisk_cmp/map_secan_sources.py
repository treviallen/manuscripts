from scipy.interpolate import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patheffects as PathEffects
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource
from numpy import arange, mean, percentile, array, unique, where, argsort, floor
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from os import path, walk, system
#from obspy.imaging.beachball import Beach
from openquake.hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser, CsvCatalogueWriter
from misc_tools import remove_last_cmap_colour
from mapping_tools import drawshapepoly, get_map_polygons, mask_outside_polygons
import matplotlib as mpl
mpl.style.use('classic')
import importlib
import codecs, sys
importlib.reload(sys) # for unicode chars


plt.rcParams['pdf.fonttype'] = 42
'''
##########################################################################################
# parse epicentres
##########################################################################################
hmtk_csv = '../catalogues/2010SHEEF/SHEEF2010_slab_hmtk.csv'
# parse HMTK csv
parser = CsvCatalogueParser(hmtk_csv)
cat = parser.read_file()
'''

##########################################################################################
# SW Can
urcrnrlat = 49
llcrnrlat = 41.0
urcrnrlon = -60.5
llcrnrlon = -81.5

lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(18,10))
plt.tick_params(labelsize=16)
ax = fig.add_subplot(111)

m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='h',area_thresh=500.)

# draw coastlines, state and country boundaries, edge of map.
#m.shadedrelief()
m.drawcoastlines(color='0.5', linewidth=0.5)
#m.fillcontinents(color=None, lake_color='lightskyblue')
m.drawstates(color='0.5', linewidth=0.5)
m.drawcountries(color='0.5', linewidth=0.5)
m.drawparallels(arange(-90.,90.,1.), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.5)
m.drawmeridians(arange(0.,360.,2.), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.5)
#m.drawmapscale(144, -34.8, 146., -38.5, 400, fontsize = 16, barstyle='fancy', zorder=100)

##########################################################################################
# plot gebco
##########################################################################################

print('Reading netCDF file...')
nc = NetCDFFile('//Users//trev//Documents//DATA//GMT//GEBCO//southern_canada_30c.nc')

zscale =20. #gray
zscale =50. #colour
try:
    data = nc.variables['elevation'][:] / zscale
    lons = nc.variables['lon'][:]
    lats = nc.variables['lat'][:]
except:
    data = nc.variables['z'][:] / zscale
    lons = nc.variables['x'][:]
    lats = nc.variables['y'][:]

# transform to metres
nx = int((m.xmax-m.xmin)/500.)+1
ny = int((m.ymax-m.ymin)/500.)+1

topodat = m.transform_scalar(data,lons,lats,nx,ny)

print('Getting colormap...')
# get colormap
cptfile = '//Users//trev//Documents//DATA//GMT//cpt//mby_topo-bath.cpt'
#cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//gray.cpt'
cmap, zvals = cpt2colormap(cptfile, 256)
cmap = remove_last_cmap_colour(cmap)
        
# make shading
print('Making map...')
ls = LightSource(azdeg = 180, altdeg = 45)
norm = mpl.colors.Normalize(vmin=-8000/zscale, vmax=5000/zscale)#myb
rgb = ls.shade(topodat, cmap=cmap, norm=norm)
im = m.imshow(rgb, alpha=0.6)

'''
##########################################################################################
# add epicentres
##########################################################################################

ncols = 18

# get year range
minyear = 10*floor(cat.data['year'][0]/10.)
maxyear = cat.data['year'][-1]
yearrng = float(round(maxyear - minyear))
#cmap = plt.cm.get_cmap('Spectral', ncols)

cptfile = '/Users/trev/Documents/DATA/GMT/cpt/temperature.cpt'
#cptfile = '/Users/trev/Documents/DATA/GMT/cpt/grayscale08.cpt'
cmap, zvals = cpt2colormap(cptfile, ncols+1, rev=False)
cmap = remove_last_cmap_colour(cmap)

cs = (cmap(arange(ncols)))
#cs = vstack((cs[0:1], cs[2:], [1., 1., 1., 1.]))

# get zorder for plotting
sortidx = argsort(argsort(cat.data['magnitude']))
for i in range(0, len(cat.data['magnitude'])): #[0:100])):
    if cat.data['magnitude'][i] >= 2.5:
        #get colour idx
        year = cat.data['year'][i]
        colidx = int(round((ncols-1) * (year-minyear) / yearrng))
        x, y = m(cat.data['longitude'][i], cat.data['latitude'][i])
        zo = sortidx[i] + 20
        plt.plot(x, y, 'o', markerfacecolor=cs[colidx], markeredgecolor='k', markeredgewidth=0.5, \
                 markersize=-4 + cat.data['magnitude'][i]*3.2, zorder=zo, alpha=0.95)
    
# make legend
legmag = [3., 5., 7.]
legh = []
for lm in legmag:
    x, y = m(0, 0)
    h = m.plot(x, y, 'ko', mfc='k', markersize=(-4 + lm*3.2), alpha=1., zorder=len(cat.data['magnitude'])+1, lw=2)
    legh.append(h[0])

l = plt.legend(legh, ('MW 3.0', 'MW 5.0', 'MW 7.0'), loc=1, numpoints=1)
l.set_zorder(len(cat.data['magnitude'])+5)
'''
##########################################################################################
# get land & lake polygons for masking
##########################################################################################
polys = get_map_polygons(m)

#mask_outside_polygon(polys[1][::-1], ax=None)
#mask_outside_polygons(polys, 'lightskyblue', plt)

# get lake ploygons
polygons = []
for polygon in m.lakepolygons:
    poly = polygon.get_coords()
    plt.fill(poly[:,0], poly[:,1], 'lightskyblue')
    polygons.append(poly)

##########################################################################################
# add source zones
##########################################################################################
import shapefile
shpfile = 'SECan_H2.shp'
sf = shapefile.Reader(shpfile)

drawshapepoly(m, plt, sf, lw=1.5, col='r')
'''
shpfile = 'swcan_crustal.shp'
sf = shapefile.Reader(shpfile)
drawshapepoly(m, plt, sf, lw=2., col='b', ls='--', polyline=True)

shpfile = 'swcan_slab.shp'
sf = shapefile.Reader(shpfile)
drawshapepoly(m, plt, sf, lw=2., col='k', ls='-', polyline=True)
'''
# set legend
m.plot([0,0],[10, 10], 'r-', lw=1.5, label='Shallow Area Sources') 
plt.legend(loc=1, fontsize=16)
'''
# label thrusts
x, y = m(-132.5, 52.5)
txt = plt.text(x, y, 'HGT', color='b', size=16, ha='right', va='top', weight='normal')
txt.set_path_effects([PathEffects.withStroke(linewidth=3, foreground='w')])

x, y = m(-129.3, 50.5)
txt = plt.text(x, y, 'WIN', color='b', size=16, ha='right', va='top', weight='normal')
txt.set_path_effects([PathEffects.withStroke(linewidth=3, foreground='w')])

x, y = m(-128.3, 49.6)
txt = plt.text(x, y, 'EXP', color='b', size=16, ha='right', va='top', weight='normal')
txt.set_path_effects([PathEffects.withStroke(linewidth=3, foreground='w')])

x, y = m(-127, 47.8)
txt = plt.text(x, y, 'JdF', color='b', size=16, ha='right', va='top', weight='normal')
txt.set_path_effects([PathEffects.withStroke(linewidth=3, foreground='w')])
'''
##########################################################################################
# add cities
##########################################################################################
	
clat = [46.500, 43.650, 45.310, 45.509, 46.800, 47.830, 45.950, 44.650, 46.233]#, 47.570]
clon = [-81.000, -79.380, -75.910, -73.554, -71.230, -69.530, -66.650, -63.600, -63.133]#, -52.720]
cloc = ['Sudbury', 'Toronto (city hall)', 'Ottawa (Kanata)', 'Montral (city hall)','Qubec City', \
        'Rivire-du-Loup','Fredericton','Halifax','Charlottetown']#, "St. John's"] # latin

cloc = ['Sudbury', 'Toronto (city hall)', 'Ottawa (Kanata)', 'Montreal (city hall)','Quebec City', \
        'Riviere-du-Loup','Fredericton','Halifax','Charlottetown']#, "St. John's"] # latin

x, y = m(clon, clat)
plt.plot(x, y, 'o', markerfacecolor='maroon', markeredgecolor='w', markeredgewidth=2, markersize=14)

# label cities
off = 0.25
for i, loc in enumerate(cloc):
    if i == 0 or i == 1 or i == 7 or i == 3:
        x, y = m(clon[i]+0.12, clat[i]+0.11)
        txt = plt.text(x, y, loc, size=18, horizontalalignment='left', weight='normal')
    elif i == 12:
        x, y = m(clon[i]-0.11, clat[i]-0.11)
        txt = plt.text(x, y, loc, size=18, ha='right', va='top', weight='normal')
    else:
        x, y = m(clon[i]-0.15, clat[i]+0.11)
        txt = plt.text(x, y, loc, size=18, horizontalalignment='right', weight='normal')
        
    txt.set_path_effects([PathEffects.withStroke(linewidth=3, foreground='w')])

'''
##########################################################################################
# add stations for 4.4
##########################################################################################

folder = '/Users/trev/Documents/Earthquake_Data/20120619.Moe/Moe_4.4/psa'
stns = []
for root, dirnames, filenames in walk(folder):
    for filename in filenames:
        if filename.endswith('.psa'):
            stns.append(filename.split('.')[1])
            
stns = unique(array(stns))

stlat = []
stlon = []
stnet = []

lines = open('/Users/trev/Documents/Code/process_waves/stationlist.dat').readlines()
for stn in stns:
    for line in lines:
        if line.startswith(stn):
            dat = line.strip().split('\t')
            lat = float(dat[5])
            lon = float(dat[4])
            net = dat[6]
            
    stlat.append(lat)
    stlon.append(lon)
    stnet.append(net)
    print(stn, net, lat)
stlat = array(stlat)
stlon = array(stlon)
#unet = unique(array(stnet))
#print(stns, stnet


# loop thru networks and plot
unet = ['AU', 'MEL', 'S', 'UOM']
sym = ['^', 'H', 'd', 's'] 
ms = [12, 13, 12, 12]
hnd = []
for i, u in enumerate(unet):
    idx = where(array(stnet) == u)[0]
    x, y = m(stlon[idx], stlat[idx])
    h = plt.plot(x, y, sym[i], markerfacecolor = 'w', markeredgecolor='k', markeredgewidth=1, markersize=ms[i])
    hnd.append(h[0])

##########################################################################################
# add stations for 5.4
##########################################################################################

folder = 'Moe_5.4/psa'
stns = []
for root, dirnames, filenames in walk(folder):
    for filename in filenames:
        if filename.endswith('.psa'):
            stns.append(filename.split('.')[1])
            
stns = unique(array(stns))

stlat = []
stlon = []
stnet = []

lines = open('/Users/trev/Documents/Code/process_waves/stationlist.dat').readlines()
for stn in stns:
    for line in lines:
        if line.startswith(stn+'\t'):
            dat = line.strip().split('\t')
            lat = float(dat[5])
            lon = float(dat[4])
            net = dat[6]
    stlat.append(lat)
    stlon.append(lon)
    stnet.append(net)
stlat = array(stlat)
stlon = array(stlon)

# loop thru networks and plot
hnd = []
for i, u in enumerate(unet):
    idx = where(array(stnet) == u)[0]
    x, y = m(stlon[idx], stlat[idx])
    h = plt.plot(x, y, sym[i], markerfacecolor = '0.1', markeredgecolor='w', markeredgewidth=1.5, markersize=ms[i])
    hnd.append(h[0])
    
plt.legend(hnd, list(unet), fontsize=14, loc=4, numpoints=1)

##########################################################################################
# annotate earthquake
##########################################################################################

eqlat = -38.252
eqlon = 146.234
blats = [-39.75]
blons = [146.5]

# line connecting epi to beachball
x, y = m([blons[0], eqlon], [blats[0], eqlat])
plt.plot(x, y, 'k-', lw=1.5)

x, y = m(eqlon, eqlat)
plt.plot(x, y, '*', markerfacecolor='r', markeredgecolor='w', markeredgewidth=.5, markersize=25)

# Add beachballs for two events
x, y = m(blons, blats)

# Two focal mechanisms for beachball routine, specified as [strike, dip, rake]
focmecs = [[122, 22, 167]]
for i in range(len(focmecs)):
    b = Beach(focmecs[i], xy=(x[i], y[i]), width=100000, linewidth=1, facecolor='k')
    b.set_zorder(10)
    ax.add_collection(b)

##########################################################################################
# draw box for fig 2
##########################################################################################

# draw box
f2lat = [-38.55, -37.75, -37.75, -38.55, -38.55]
f2lon = [145.7, 145.7, 146.75, 146.75, 145.7]
x, y = m(f2lon, f2lat)
plt.plot(x, y, 'k--', lw=1.5)

# write text
ft2lat = min(f2lat)+0.
ft2lon = max(f2lon)+0.1
x, y = m(ft2lon, ft2lat)
plt.text(x, y, 'Figure 2', size=14, horizontalalignment='left', weight='bold', style='italic')
'''
##########################################################################################
# make map inset
##########################################################################################

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
axins = zoomed_inset_axes(ax, 0.055, loc=4)

m2 = Basemap(width=6000000,height=5500000,
        rsphere=(6378137.00,6356752.3142), \
        resolution='l',area_thresh=20000.,projection='lcc', \
        lat_1=49.,lat_2=77,lat_0=63,lon_0=-90.)
            
m2.drawmapboundary(fill_color='0.8')
m2.fillcontinents(color='w', lake_color='0.8')
m2.drawcoastlines()
m2.drawcountries()
m2.drawstates()

# fill main area
xv = array([llcrnrlon, llcrnrlon, urcrnrlon, urcrnrlon, llcrnrlon])
yv = array([llcrnrlat, urcrnrlat, urcrnrlat, llcrnrlat, llcrnrlat])
x, y = m2(xv, yv)
plt.plot(x, y, 'r-', lw=1.5)

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
'''
#cb = plt.colorbar()# ticks=ticks,
ticks = arange(0, 19, 2)
years = float(minyear) + ticks*10
labels = [str('%0.0f' % x) for x in years]

# normalise
norm = mpl.colors.Normalize(vmin=minyear, vmax=maxyear)

cax = fig.add_axes([0.78,0.3,0.02,0.4]) # setup colorbar axes.
cb = colorbar.ColorbarBase(cax, cmap=cmap, orientation='vertical', alpha=0.95, norm=norm) # 
cb.set_ticks(years)
cb.set_ticklabels(labels)
cb.set_label('Year of Earthquake', rotation=270, labelpad=20, fontsize=15)
'''

plt.savefig('secan_source_map.png', format='png', bbox_inches='tight', dpi=800)
plt.show()
