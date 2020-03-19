from matplotlib.mlab import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource
from numpy import arange, mean, percentile, array, unique, where, argsort, vstack
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from os import path, walk, system
#from obspy.imaging.beachball import Beach
from datetime import datetime
from gmt_tools import cpt2colormap

plt.rcParams['pdf.fonttype'] = 42

##########################################################################################
# parse epicentres
##########################################################################################
epifile = 'Moe_Events.csv'
lines = open(epifile).readlines()[1:]
evdt = []
evla = []
evlo = []
evdp = []
evml = []
for line in lines:
    dat = line.split(',')
    evdt.append(datetime.strptime(dat[0][0:-2], '%Y-%m-%d %H%M %S'))
    evla.append(float(dat[1]))
    evlo.append(float(dat[2]))
    evdp.append(float(dat[3]))
    evml.append(float(dat[5]))

##########################################################################################
# set up map
##########################################################################################

f2lat = [-38.39, -38.1, -38.1, -38.39, -38.39]
f2lon = [146.1, 146.1, 146.45, 146.45, 146.1]

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
#m.drawcoastlines()
#m.drawstates()
#m.drawcountries()
#m.drawmapboundary(fill_color='lightblue')
#m.drawmapboundary(fill_color='0.8')
m.drawparallels(arange(-90.,90.,.1), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.45', linewidth=0.75)
m.drawmeridians(arange(0.,360.,.1), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.45', linewidth=0.75)
m.drawmapscale(146.16, -38.12, 146.16, -38.12, 10, fontsize = 16, barstyle='fancy', zorder=500)

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
# plot epicentres
##########################################################################################
fig = plt.figure(1, figsize=(10, 10))
ncols = 10
#cmap = plt.cm.get_cmap('Spectral', ncols)

cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/grayscale08.cpt'
cmap, zvals = cpt2colormap(cptfile, ncols, rev=True)

cs = (cmap(arange(ncols)))
cs = cs[::-1]
#cs = vstack((cs[0:3], cs[4:]))

zmin = 8 
zmax = 18
zrng = 10
# get zorder for plotting
sortidx = argsort(argsort(array(evml)))
for i, ed in enumerate(evdp):
    if ed != 10.:
       #get colour idx
       #month = int(ed.strftime('%m'))
       colidx = int(ed - zmin)
       x, y = m(evlo[i], evla[i])
       zo = sortidx[i] + 20
       if colidx > 9:
           colidx = 9
       if colidx < 0:
           colidx == 0
       color = cs[colidx]
       p = m.plot(x, y, 'o', color=color, markersize=(6 + evml[i]*2.2), alpha=1., zorder=zo)

# make legend
legmag = [1.5, 3.5, 5.5]
legh = []
for lm in legmag:
    x, y = m(0, 0)
    h = m.plot(x, y, 'ko', mfc='w', markersize=(6 + lm*2.2), alpha=1., zorder=zo, lw=2)
    legh.append(h[0])

plt.legend(legh, ('ML 1.5', 'ML 4.0', 'ML 5.5'), loc=1, numpoints=1)
    
##########################################################################################
# add cities
##########################################################################################

clat = [-38.125, -38.235] #-37.814, 
clon = [146.263, 146.395] #144.964, 
cloc = ['Moe', 'Morwell'] # 'Melbourne', 
x, y = m(clon, clat)
plt.plot(x, y, 'o', markerfacecolor='k', markeredgecolor='w', markeredgewidth=2, markersize=14)

# label cities
for i, loc in enumerate(cloc):
    x, y = m(clon[i]+0.005, clat[i]+0.005)
    plt.text(x, y, loc, size=18, horizontalalignment='left', weight='normal', zorder=1000)

##########################################################################################
# add powerstations
##########################################################################################
"""
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
"""

##########################################################################################
# annotate earthquake
##########################################################################################
"""
eqlat = -38.252
eqlon = 146.234
blats = [-39.75]
blons = [146.5]
'''
# line connecting epi to beachball
x, y = m([blons[0], eqlon], [blats[0], eqlat])
plt.plot(x, y, 'k-', lw=1.5)
'''
x, y = m(eqlon, eqlat)
plt.plot(x, y, '*', markerfacecolor='r', markeredgecolor='w', markeredgewidth=.5, markersize=50)
'''
# Add beachballs for two events
x, y = m(blons, blats)

# Two focal mechanisms for beachball routine, specified as [strike, dip, rake]
focmecs = [[122, 22, 167]]
for i in range(len(focmecs)):
    b = Beach(focmecs[i], xy=(x[i], y[i]), width=100000, linewidth=1, facecolor='k')
    b.set_zorder(10)
    ax.add_collection(b)
'''
"""
##########################################################################################
# make map inset
##########################################################################################
"""
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
"""
#cb = plt.colorbar()# ticks=ticks,

cmap2 = cmap.from_list('custom', cs, N=ncols) #, ncols
#norm = mpl.colors.BoundaryNorm(ticks, cmap.N)
 
cax = fig.add_axes([0.76,0.3,0.02,0.4]) # setup colorbar axes.
norm = mpl.colors.Normalize(vmin=zmin, vmax=zmax)
cb = colorbar.ColorbarBase(cax, cmap=cmap2, orientation='vertical', norm=norm) # norm=norm
cb.ax.tick_params(labelsize=16)
cb.set_label('Earthquake Depth (km)', fontsize=18, rotation=270, labelpad=25)

plt.savefig('moe_map_depth.png', format='png', bbox_inches='tight', dpi=300)
plt.show()
