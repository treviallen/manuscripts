from scipy.interpolate import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource
from misc_tools import remove_last_cmap_colour
from numpy import arange, mean, percentile, array, unique, where, argsort
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from os import path, walk, system, getcwd
import shapefile
#from obspy.imaging.beachball import Beach

plt.rcParams['pdf.fonttype'] = 42
import matplotlib as mpl
mpl.style.use('classic')

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
m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='i',area_thresh=150)

# draw coastlines, state and country boundaries, edge of map.
#m.drawcoastlines()
m.drawstates()
m.drawcountries()
#m.drawmapboundary(fill_color='lightblue')
#m.drawmapboundary(fill_color='0.8')
m.drawparallels(arange(-90.,90.,.2), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.45', linewidth=0.75)
m.drawmeridians(arange(0.,360.,.2), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.45', linewidth=0.75)
m.drawmapscale(146.5, -38.5, 146.5, -38.5, 40, fontsize = 17, barstyle='fancy', zorder=100)

##########################################################################################
# plot gebco
##########################################################################################

print('Reading netCDF file...')
if getcwd().startswith('/nas'):
    netcdffile = 'topo/srtm_66_20.grd'
else:
    netcdffile = '//Users//trev//Documents//DATA//SRTM03//srtm_66_20//srtm_66_20.grd'
nc = NetCDFFile(netcdffile)

zscale =20. #colour
data = nc.variables['z'][:] / zscale
lons = nc.variables['x'][:]
lats = nc.variables['y'][:]
'''
print('Reading netCDF file...')
nc = NetCDFFile('//Users//trev//Documents//DATA//GMT//GEBCO//se_aus_gebco2014.nc')

zscale =20. #gray
zscale =30. #colour
data = nc.variables['z'][:] / zscale
lons = nc.variables['lon'][:]
lats = nc.variables['lat'][:]
'''
# transform to metres      
nx = int((m.xmax-m.xmin)/30.)+1
ny = int((m.ymax-m.ymin)/30.)+1

topodat = m.transform_scalar(data,lons,lats,nx,ny)

print('Getting colormap...')
# get colormap
if getcwd().startswith('/nas'):
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/wiki-2.0.cpt'
else:
    cptfile = '//Users//trev//Documents//DATA//GMT//cpt//wiki-2.0.cpt'
cmap, zvals = cpt2colormap(cptfile, 117)
print(len(zvals))
cmap = remove_last_cmap_colour(cmap)
        
# make shading
print('Making map...')
ls = LightSource(azdeg = 180, altdeg = 5)
#norm = mpl.colors.Normalize(vmin=-8000/zscale, vmax=5000/zscale)
norm = mpl.colors.Normalize(vmin=-1000/zscale, vmax=1900/zscale) #wiki
rgb = ls.shade(topodat, cmap=cmap, norm=norm)
im = m.imshow(rgb, alpha=1.)

##########################################################################################
# add simple faults
##########################################################################################
if getcwd().startswith('/nas'):
    nfsmshp = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/source_models/faults/FSM/FSD_simple_faults.shp'
else:
    nfsmshp = '/Users/trev/Documents/Geoscience_Australia/NSHA2018/source_models/faults/FSM/FSD_simple_faults.shp'
    
sf = shapefile.Reader(nfsmshp)
shapes = sf.shapes()
records = sf.records()
'''
lt_rates = get_field_data(sf, 'SL_RT_LT', 'float') # long term slip rate
lt_rates = array(lt_rates)
st_rates = get_field_data(sf, 'SL_RT_ST', 'float') # short term slip rate
st_rates = array(st_rates)

minslip = 0.
maxslip = 160.
slip_diff = maxslip - minslip

sortidx = argsort(argsort(lt_rates))
'''
i = 0
for shape in shapes:
    lons = []
    lats = []
    for xy in shape.points:
        lons.append(xy[0])
        lats.append(xy[1])
    
    x, y = m(lons, lats)
    
    '''
    if ltr > maxslip:
        ltr = maxslip
    
    # set colour by 
    colidx = int(round((cmap.N-1) * (ltr-minslip) / slip_diff))
    '''
    if i == 0:
        #plt.plot(x, y, 'r-', lw=2.5, zorder=sortidx[i], label='Neotectonic Features')
        plt.plot(x, y, 'r-', lw=2.5)
    else:
        plt.plot(x, y, 'r-', lw=2.5)
    
    i += 1


##########################################################################################
# add cities
##########################################################################################
import matplotlib.patheffects as PathEffects
pe = [PathEffects.withStroke(linewidth=2.5, foreground="w")]


clat = [-38.197, -38.125, -38.235, -38.161] #-37.814, 
clon = [146.541, 146.263, 146.395, 145.932] #144.964, 
cloc = ['Traralgon', 'Moe', 'Morewell', 'Warragul'] # 'Melbourne', 
x, y = m(clon, clat)
plt.plot(x, y, 'o', markerfacecolor='maroon', markeredgecolor='w', markeredgewidth=2, markersize=14)

# label cities
for i, loc in enumerate(cloc):
    x, y = m(clon[i]+0.012, clat[i]+0.008)
    plt.text(x, y, loc, size=18, horizontalalignment='left', weight='normal', path_effects=pe)

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
    plt.text(x, y, loc, size=16, horizontalalignment='left', weight='normal', style='italic', path_effects=pe)

##########################################################################################
# add stations for 4.4
##########################################################################################

folder = 'psa/20120720/'
stns = []
for root, dirnames, filenames in walk(folder):
    for filename in filenames:
        if filename.endswith('.psa'):
            stns.append(filename.split('.')[1])
            
stns = unique(array(stns))

stlat = []
stlon = []
stnet = []

if getcwd().startswith('/nas'):
    lines = open('/nas/users/u56903/unix/Code/my_codes/stationlist.dat').readlines()
else:
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
#unet = unique(array(stnet))
#print(stns, stnet

# loop thru networks and plot
unet = ['AU', 'MEL', 'S', 'UM']
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
            plt.text(x, y, stn, size=14, ha='right', weight='normal', style='normal', path_effects=pe)
        elif stn == 'MOE3':
            x, y = m(stlon[i]-0.01, stlat[i]+0.01)
            plt.text(x, y, stn, size=14, ha='right', weight='normal', style='normal', path_effects=pe)
        elif stn == 'MOE4':
            x, y = m(stlon[i]-0.01, stlat[i]+0.01)
            plt.text(x, y, stn, size=14, ha='right', weight='normal', style='normal', path_effects=pe)
        else:
            x, y = m(stlon[i]+0.015, stlat[i]-0.017)
            plt.text(x, y, stn, size=14, ha='left', weight='normal', style='normal', path_effects=pe)
        
print(stns)

##########################################################################################
# add stations for 5.4
##########################################################################################

folder = 'psa/20120619/'
stns = []
for root, dirnames, filenames in walk(folder):
    for filename in filenames:
        if filename.endswith('.psa'):
            stns.append(filename.split('.')[1])
            
stns = unique(array(stns))

stlat = []
stlon = []
stnet = []

if getcwd().startswith('/nas'):
    lines = open('/nas/users/u56903/unix/Code/my_codes/stationlist.dat').readlines()
else:
    lines = open('/Users/trev/Documents/Code/my_codes/stationlist.dat').readlines()
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
    h = plt.plot(x, y, sym[i], markerfacecolor = '0.1', mec='w', mew=1.5, markersize=ms[i], label=unet[i])
    hnd.append(h[0])

# label stations
for i, stn in enumerate(stns):
    if stlat[i] <= urcrnrlat and stlat[i] >= llcrnrlat \
       and stlon[i] <= urcrnrlon and stlon[i] >= llcrnrlon:
        if stn == 'NARR':
            x, y = m(stlon[i]-0.015, stlat[i]-0.017)
            plt.text(x, y, stn, size=14, ha='right', weight='normal', style='normal', path_effects=pe)
        else:
            x, y = m(stlon[i]+0.015, stlat[i]-0.017)
            plt.text(x, y, stn, size=14, ha='left', weight='normal', style='normal', path_effects=pe)

# add power stations?
plt.legend(fontsize=16, loc=1, numpoints=1)

##########################################################################################
# annotate earthquake
##########################################################################################

# plt M4.4 - using SRC locs
eqlat = -38.229
eqlon = 146.218
blats = [-39.75]
blons = [146.5]

x, y = m(eqlon, eqlat)
plt.plot(x, y, '*', markerfacecolor='darkred', markeredgecolor='w', markeredgewidth=1.5, markersize=10*4.2)

# plt M5.4
eqlat = -38.279
eqlon = 146.269
blats = [-39.75]
blons = [146.5]

x, y = m(eqlon, eqlat)
plt.plot(x, y, '*', markerfacecolor='r', markeredgecolor='w', markeredgewidth=1.5, markersize=10*5.4)
##########################################################################################
# make map inset
##########################################################################################

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
axins = zoomed_inset_axes(ax, 0.0065, loc=2)

m2 = Basemap(projection='merc',\
            llcrnrlon=111,llcrnrlat=-45, \
            urcrnrlon=156,urcrnrlat=-9,\
            rsphere=6371200.,resolution='l',area_thresh=10000)
            
m2.drawmapboundary(fill_color='0.8')
m2.fillcontinents(color='w', lake_color='0.8') #, zorder=0)
m2.drawcoastlines()
m2.drawcountries()
m2.drawstates()

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
