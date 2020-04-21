from scipy.interpolate import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource
from numpy import arange, mean, percentile, array, unique, where
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from misc_tools import remove_last_cmap_colour
from mapping_tools import drawmanualshapepoly, drawshapepoly
from os import path, walk, system
from obspy.imaging.beachball import beach, beachball
import shapefile

plt.rcParams['pdf.fonttype'] = 42
import matplotlib as mpl
mpl.style.use('classic')


##########################################################################################
# set up map
##########################################################################################

urcrnrlat = -34.1
llcrnrlat = -43.0
urcrnrlon = 151.
llcrnrlon = 141.0
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(18,10))
plt.tick_params(labelsize=16)
ax = fig.add_subplot(111)

m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='i',area_thresh=150)

# draw coastlines, state and country boundaries, edge of map.
m.drawcoastlines()
m.drawstates()
m.drawcountries()
m.drawparallels(arange(-90.,90.,2.), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,2.), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmapscale(144, -34.8, 146., -38.5, 400, fontsize = 16, barstyle='fancy', zorder=100)

##########################################################################################
# plot gebco
##########################################################################################

print('Reading netCDF file...')
nc = NetCDFFile('//Users//trev//Documents//DATA//GMT//GEBCO//se_aus_gebco2014.nc')

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
rgb = ls.shade(topodat, cmap=cmap, norm=norm)
im = m.imshow(rgb, alpha=1.)

##########################################################################################
# add basins
##########################################################################################
import matplotlib.patheffects as PathEffects
path_effects=[PathEffects.withStroke(linewidth=2.5, foreground="w")]

shpfile = 'gis/gippslans_bass_basin_wgs84.shp'
sf = shapefile.Reader(shpfile)
drawshapepoly(m, plt, sf, col='r', fillcolor='r', edgecolor='r', lw=0.75, alpha=0.2, fillshape = True)

# label polygons
x, y = m(146, -40.1)
plt.text(x, y, 'Bass\nBasin', size=12, c='r', ha='center', va='top', weight='normal', style='italic', path_effects=path_effects)

x, y = m(148.5, -39.)
plt.text(x, y, 'Gippsland\nBasin', size=12, c='r', ha='center', weight='normal', style='italic', path_effects=path_effects)

##########################################################################################
# add cities
##########################################################################################
path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")]

clat = [-37.814, -42.882, -35.282]
clon = [144.964, 147.324, 149.129]
cloc = ['Melbourne', 'Hobart', 'Canberra']
x, y = m(clon, clat)
plt.plot(x, y, 'o', markerfacecolor='maroon', markeredgecolor='w', markeredgewidth=2, \
         markersize=14)

# label cities
off = 0.25
for i, loc in enumerate(cloc):
    if i == 1:
        x, y = m(clon[i]+0.12, clat[i]+0.1)
        plt.text(x, y, loc, size=18, horizontalalignment='left', weight='normal', path_effects=path_effects)
    else:
        x, y = m(clon[i]-0.12, clat[i]+0.1)
        plt.text(x, y, loc, size=18, horizontalalignment='right', weight='normal', path_effects=path_effects)

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

lines = open('/Users/trev/Documents/Code/my_codes/stationlist.dat').readlines()
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
    print(stn, net, lat)
stlat = array(stlat)
stlon = array(stlon)
#unet = unique(array(stnet))
#print(stns, stnet


# loop thru networks and plot
unet = ['AU', 'MEL', 'S1', 'UM']
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

folder = '/Users/trev/Documents/Earthquake_Data/20120619.Moe/Moe_5.4/psa'
stns = []
for root, dirnames, filenames in walk(folder):
    for filename in filenames:
        if filename.endswith('.psa'):
            stns.append(filename.split('.')[1])
            
stns = unique(array(stns))

stlat = []
stlon = []
stnet = []

lines = open('/Users/trev/Documents/Code/my_codes/stationlist.dat').readlines()
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
blats = [-39.6]
blons = [147.]

# line connecting epi to beachball
x, y = m([blons[0], eqlon], [blats[0], eqlat])
plt.plot(x, y, 'k-', lw=1.5)

x, y = m(eqlon, eqlat)
plt.plot(x, y, '*', markerfacecolor='r', markeredgecolor='w', markeredgewidth=.5, markersize=25)

# Add beachballs for two events
x, y = m(blons, blats)

# Two focal mechanisms for beachball routine, specified as [strike, dip, rake]
#beachball(mt, size=200, linewidth=2, facecolor='dodgerblue',outfile=outfile)
focmecs = [[122, 22, 167]]
for i in range(len(focmecs)):
    b = beach(focmecs[i], xy=(x[i], y[i]), width=100000, linewidth=1, facecolor='k')
    b.set_zorder(1000)
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

##########################################################################################
# make map inset
##########################################################################################

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
axins = zoomed_inset_axes(ax, 0.06, loc=3)

m2 = Basemap(projection='merc',\
            llcrnrlon=111,llcrnrlat=-45, \
            urcrnrlon=156,urcrnrlat=-9,\
            rsphere=6371200.,resolution='l',area_thresh=10000)
            
m2.drawmapboundary(fill_color='0.8')
m2.fillcontinents(color='w', lake_color='0.8') #, zorder=0)
m2.drawcoastlines()
m2.drawcountries()
m2.drawstates()

# fill main area
xv = array([llcrnrlon, llcrnrlon, urcrnrlon, urcrnrlon, llcrnrlon])
yv = array([llcrnrlat, urcrnrlat, urcrnrlat, llcrnrlat, llcrnrlat])
x, y = m2(xv, yv)
plt.plot(x, y, 'k-', lw=1.5)

##########################################################################################
# label states
##########################################################################################

state = ['WA', 'NT', 'SA', 'QLD', 'NSW', 'VIC', 'TAS']
slat = [-26, -21.0, -29.5, -23.0, -32.5, -37.1, -44.2]
slon = [122, 133.5, 135.0, 144.5, 146.5, 143.6, 150.0]
for i, st in enumerate(state):
    x, y = m2(slon[i], slat[i])
    plt.text(x, y, st, size=11, horizontalalignment='center', verticalalignment='center', weight='normal')

plt.savefig('moe_map.png', format='png', bbox_inches='tight', dpi=300)
plt.show()
