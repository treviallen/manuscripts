from matplotlib.mlab import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource
from numpy import arange, mean, percentile, array, unique, where
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from os import path, walk, system
#from obspy.imaging.beachball import Beach

plt.rcParams['pdf.fonttype'] = 42


##########################################################################################
# set up map
##########################################################################################

f2lat = [-38.55, -37.75, -37.75, -38.55, -38.55]
f2lon = [145.7, 145.7, 146.75, 146.75, 145.7]

urcrnrlat = max(f2lat)
llcrnrlat = min(f2lat)
urcrnrlon = max(f2lon)
llcrnrlon = min(f2lon)
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(18,11))
plt.tick_params(labelsize=16)
ax = fig.add_subplot(111)
ax.tick_params(axis='x', labelsize=16)
ax.tick_params(axis='y', labelsize=16)

#lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
m = Basemap(projection='merc', \
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='f',area_thresh=100)

# draw coastlines, state and country boundaries, edge of map.
m.drawcoastlines()
m.drawstates()
m.drawcountries()
m.drawmapboundary(fill_color='lightblue')
#m.drawmapboundary(fill_color='0.8')
m.drawparallels(arange(-90.,90.,.2), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.45', linewidth=0.75)
m.drawmeridians(arange(0.,360.,.2), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.45', linewidth=0.75)
m.drawmapscale(146.5, -38.5, 146.5, -38.5, 40, fontsize = 17, barstyle='fancy', zorder=100)

##########################################################################################
# plot gebco
##########################################################################################

print 'Reading netCDF file...'
netcdffile = '//Users//tallen//Documents//DATA//SRTM03//srtm_66_20//srtm_66_20.grd'
nc = NetCDFFile(netcdffile)

zscale =2. #gray
zscale =10. #colour
data = nc.variables['z'][:] / zscale
lons = nc.variables['x'][:]
lats = nc.variables['y'][:]

# transform to metres
nx = int((m.xmax-m.xmin)/30.)+1
ny = int((m.ymax-m.ymin)/30.)+1

topodat = m.transform_scalar(data,lons,lats,nx,ny)

print 'Getting colormap...'
# get colormap
cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//mby_topo-bath_mod.cpt'
#cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//gray.cpt'
cmap, zvals = cpt2colormap(cptfile, 256)
#cmap = cm.get_cmap('terrain', 256)

# make shading
print 'Making map...'
ls = LightSource(azdeg = 180, altdeg = 45)
norm = mpl.colors.Normalize(vmin=-8000/zscale, vmax=5000/zscale)#myb
rgb = ls.shade(topodat, cmap=cmap, norm=norm)
im = m.imshow(rgb)

##########################################################################################
# add cities
##########################################################################################

clat = [-38.197, -38.125, -38.235, -38.161] #-37.814, 
clon = [146.541, 146.263, 146.395, 145.932] #144.964, 
cloc = ['Traralgon', 'Moe', 'Morewell', 'Warragul'] # 'Melbourne', 
x, y = m(clon, clat)
plt.plot(x, y, 'o', markerfacecolor='maroon', markeredgecolor='w', markeredgewidth=2, markersize=14)

# label cities
for i, loc in enumerate(cloc):
    x, y = m(clon[i]+0.012, clat[i]+0.008)
    plt.text(x, y, loc, size=18, horizontalalignment='left', weight='normal')

##########################################################################################
# add powerstations
##########################################################################################

plat = [-38.255, -38.272, -38.177]
plon = [146.578, 146.391, 146.342]
x, y = m(plon, plat)
ph = plt.plot(x, y, 'p', markerfacecolor='b', markeredgecolor='w', markeredgewidth=1, markersize=17)

# label power plants
ploc = ['Loy Yang', 'Hazelwood', 'Yallourn']
for i, loc in enumerate(ploc):
    if i == 2:
        x, y = m(plon[i]+0.012, plat[i]+0.01)
    else:
        x, y = m(plon[i]+0.012, plat[i]-0.02)
    plt.text(x, y, loc, size=16, horizontalalignment='left', weight='normal', style='italic')

##########################################################################################
# add stations for 4.4
##########################################################################################

folder = 'Moe_4.4/psa'
stns = []
for root, dirnames, filenames in walk(folder):
    for filename in filenames:
        if filename.endswith('.psa'):
            stns.append(filename.split('.')[1])
            
stns = unique(array(stns))

stlat = []
stlon = []
stnet = []

lines = open('/Users/tallen/Documents/Code/process_waves/stationlist.dat').readlines()
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
#unet = unique(array(stnet))
#print stns, stnet


# loop thru networks and plot
unet = ['AU', 'MEL', 'S', 'UOM']
sym = ['^', 'H', 'd', 's'] 
ms = [14, 15, 14, 13]
for i, u in enumerate(unet):
    idx = where(array(stnet) == u)[0]
    x, y = m(stlon[idx], stlat[idx])
    h = plt.plot(x, y, sym[i], markerfacecolor = 'w', markeredgecolor='k', markeredgewidth=1, markersize=ms[i])

# label stations
for i, stn in enumerate(stns):
    if stlat[i] <= urcrnrlat and stlat[i] >= llcrnrlat \
       and stlon[i] <= urcrnrlon and stlon[i] >= llcrnrlon:
        if stn == 'NARR':
            x, y = m(stlon[i]-0.015, stlat[i]-0.017)
            plt.text(x, y, stn, size=14, ha='right', weight='normal', style='normal')
        elif stn == 'MOE3':
            x, y = m(stlon[i]-0.01, stlat[i]+0.01)
            plt.text(x, y, stn, size=14, ha='right', weight='normal', style='normal')
        elif stn == 'MOE4':
            x, y = m(stlon[i]-0.01, stlat[i]+0.01)
            plt.text(x, y, stn, size=14, ha='right', weight='normal', style='normal')
        else:
            x, y = m(stlon[i]+0.015, stlat[i]-0.017)
            plt.text(x, y, stn, size=14, ha='left', weight='normal', style='normal')
        
print stns

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

lines = open('/Users/tallen/Documents/Code/process_waves/stationlist.dat').readlines()
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
stlat = array(stlat)
stlon = array(stlon)

# loop thru networks and plot
hnd = []
for i, u in enumerate(unet):
    idx = where(array(stnet) == u)[0]
    x, y = m(stlon[idx], stlat[idx])
    h = plt.plot(x, y, sym[i], markerfacecolor = '0.1', markeredgecolor='w', markeredgewidth=1.5, markersize=ms[i])
    hnd.append(h[0])

# label stations
for i, stn in enumerate(stns):
    if stlat[i] <= urcrnrlat and stlat[i] >= llcrnrlat \
       and stlon[i] <= urcrnrlon and stlon[i] >= llcrnrlon:
        if stn == 'NARR':
            x, y = m(stlon[i]-0.015, stlat[i]-0.017)
            plt.text(x, y, stn, size=14, ha='right', weight='normal', style='normal')
        else:
            x, y = m(stlon[i]+0.015, stlat[i]-0.017)
            plt.text(x, y, stn, size=14, ha='left', weight='normal', style='normal')

# add power stations?
plt.legend(hnd, list(unet), fontsize=16, loc=1, numpoints=1)

##########################################################################################
# annotate earthquake
##########################################################################################

# plt M5.4
eqlat = -38.229
eqlon = 146.218
blats = [-39.75]
blons = [146.5]

x, y = m(eqlon, eqlat)
plt.plot(x, y, '*', markerfacecolor='darkred', markeredgecolor='w', markeredgewidth=1.5, markersize=10*4.2)

# plt M5.4
eqlat = -38.252
eqlon = 146.234
blats = [-39.75]
blons = [146.5]

x, y = m(eqlon, eqlat)
plt.plot(x, y, '*', markerfacecolor='r', markeredgecolor='w', markeredgewidth=1.5, markersize=10*5.4)
##########################################################################################
# make map inset
##########################################################################################

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
axins = zoomed_inset_axes(ax, 0.0093, loc=2)

m2 = Basemap(projection='merc',\
            llcrnrlon=111,llcrnrlat=-45, \
            urcrnrlon=156,urcrnrlat=-9,\
            rsphere=6371200.,resolution='l',area_thresh=10000)
            
#m2 = Basemap(llcrnrlon=-20,llcrnrlat=3,urcrnrlon=0,urcrnrlat=18, ax=axins)
m2.drawmapboundary(fill_color='0.8')
m2.fillcontinents(color='w', lake_color='0.8', zorder=0)
m2.drawcoastlines()
m2.drawcountries()
m2.drawstates(color='0.4', linewidth=0.75)

# add torrens island
plat = -34.806
plon = 138.525
x, y = m2(plon, plat)
ph = plt.plot(x, y, 'p', markerfacecolor='b', markeredgecolor='w', markeredgewidth=1, markersize=14)
x, y = m2(plon-1.5, plat-4)
plt.text(x, y, 'Torrens\nIsland', size=13, horizontalalignment='right', weight='normal', style='italic')


# fill main area
xv = array([llcrnrlon, llcrnrlon, urcrnrlon, urcrnrlon, llcrnrlon])
yv = array([llcrnrlat, urcrnrlat, urcrnrlat, llcrnrlat, llcrnrlat])
x, y = m2(xv, yv)
plt.fill(x, y, 'r', edgecolor='r')

##########################################################################################
# label states
##########################################################################################

state = ['WA', 'NT', 'SA', 'QLD', 'NSW', 'VIC', 'TAS']
slat = [-26, -21.0, -29.5, -23, -32.5, -37.1, -43.5]
slon = [122, 133.5, 135.0, 145, 146.5, 143.5, 150.0]
for i, st in enumerate(state):
    x, y = m2(slon[i], slat[i])
    plt.text(x, y, st, size=11, horizontalalignment='center', verticalalignment='center', weight='normal')

plt.savefig('moe_map_zoom.png', format='png', bbox_inches='tight', dpi=150)
plt.show()
