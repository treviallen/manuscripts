#matplotlib quiver example

from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource, ListedColormap
from numpy import arange, mean, percentile, array, unique, where, argsort, floor, sqrt, delete, argsort
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from os import path, walk, system
from shapely.geometry import Polygon, Point
import pickle
#from obspy.imaging.beachball import Beach
#from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser, CsvCatalogueWriter
from misc_tools import remove_last_cmap_colour, listdir_extension, remove_first_cmap_colour, dictlist2array
from mapping_tools import drawshapepoly, get_field_data, drawoneshapepoly, distance, reckon, map_fault_dip_dirn
import shapefile
from io_catalogues import parse_ga_event_query

mpl.style.use('classic')

##########################################################################################
#108/152/-44/-8
urcrnrlat = -29.
llcrnrlat = -35.5
urcrnrlon = 121
llcrnrlon = 114
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(18,10))
plt.tick_params(labelsize=13)
ax = fig.add_subplot(111)

m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='f',area_thresh=100.) #i

# draw coastlines, state and country boundaries, edge of map.
#m.shadedrelief()
#m.bluemarble()
#m.etopo()
m.drawcoastlines()
m.drawstates()
m.drawcountries()
#m.fillcontinents(color='w',lake_color='w')
#m.drawmapboundary(fill_color='w')
m.drawparallels(arange(-90.,90.,1.), labels=[1,0,0,0],fontsize=13, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,1.), labels=[0,0,0,1], fontsize=13, dashes=[2, 2], color='0.5', linewidth=0.75)
#m.drawmapscale(144, -34.8, 146., -38.5, 400, fontsize = 16, barstyle='fancy', zorder=100)

##########################################################################################
# plot gebco
##########################################################################################

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
ls = LightSource(azdeg = 100, altdeg = 10)
#norm = mpl.colors.Normalize(vmin=-8000/zscale, vmax=5000/zscale)#myb
#norm = mpl.colors.Normalize(vmin=-1000/zscale, vmax=1900/zscale)#wiki
norm = mpl.colors.Normalize(vmin=-2000/zscale, vmax=3500/zscale)#wiki

rgb = ls.shade(topodat.data, cmap=cmap, norm=norm)
im = m.imshow(rgb, alpha=1.0)

##########################################################################################
# add earthquakes
##########################################################################################
import matplotlib.patheffects as PathEffects
path_effects=[PathEffects.withStroke(linewidth=2.5, foreground="w")]

##########################################################################################
# set colours
##########################################################################################

ncols = 12

# get year range
minyear = 1960
maxyear = 2020
yearrng = float(round(maxyear - minyear))

#cptfile = '/Users/trev/Documents/DATA/GMT/cpt/temperature.cpt'
cptfile = '//Users//trev//Documents//DATA//GMT//cpt//keshet.cpt'
cmap, zvals = cpt2colormap(cptfile, ncols+1, rev=True)
cmap = remove_first_cmap_colour(cmap)

cs = (cmap(arange(ncols)))

##########################################################################################
# add earthquakes
##########################################################################################

evdict = parse_ga_event_query('earthquakes_export.csv')

lons = dictlist2array(evdict, 'lon')
lats = dictlist2array(evdict, 'lat')
mags = dictlist2array(evdict, 'mag')
years = dictlist2array(evdict, 'year')

# get zorder for plotting
sortidx = argsort(argsort(mags))
for i in range(0, len(mags)): #[0:100])):
    #get colour idx
    colidx = int(floor((ncols) * (years[i]-minyear+0.1) / yearrng))
    if colidx > ncols-1:
        colidx = ncols-1
        
    x, y = m(lons[i], lats[i])
    zo = sortidx[i] + 20
    plt.plot(x, y, 'o', mfc=list(cs[colidx][0:3]), markeredgecolor='k', markeredgewidth=0.25, \
             markersize=(3. * mags[i] - 2), zorder=zo, alpha=0.8)
 
# make legend
legmag = [2.5, 4.5, 6.5]
legh = []
for lm in legmag:
    x, y = m(0, 0)
    h = m.plot(x, y, 'o', mfc='k', mec='w', mew=0.25, markersize=(3 * lm - 2), alpha=1., zorder=len(mags)+1)
    legh.append(h[0])

legtxt = ('$\mathregular{M}$ 2.5', '$\mathregular{M}$ 4.5', '$\mathregular{M}$ 6.5')
plt.rcParams['legend.title_fontsize'] = 15
l = plt.legend(legh, legtxt, loc=3, numpoints=1, fontsize=13, title="Magnitude", labelspacing=0.75)
l.set_zorder(len(mags)+5)

##########################################################################################
# add stations
##########################################################################################
from data_fmt_tools import parse_iris_stationlist

netfiles = ['/Users/trev/Documents/Networks/AU/au-gmap-stations-no_nwoa.txt', \
            '/Users/trev/Documents/Networks/S1/s1-gmap-stations.txt', \
            '/Users/trev/Documents/Networks/IU/iu-gmap-stations.txt', \
            '/Users/trev/Documents/Networks/SWAN/2p-gmap-stations.txt', \
            '/Users/trev/Documents/Networks/SWAN/2p-gmap-stations_ar.txt']
            
syms = ['^', 's', 'H', 'D', 'v']
symlab = ['AU', 'S1', 'IU', '2P', '2P(AR)']
msize = [9, 8, 11, 8, 9]
zorder = len(mags)+100
for i, netfile in enumerate(netfiles):
# map stas
    stas = parse_iris_stationlist(netfile)

    stlons = dictlist2array(stas, 'lon')
    stlats = dictlist2array(stas, 'lat')
    
    x,y = m(stlons, stlats)
    plt.plot(x, y, syms[i], mfc='w', ms=msize[i], mew=1, label=symlab[i], zorder=zorder+i)

# add legend
plt.legend(loc=4, fontsize=13, title="Networks", numpoints=1)

#re-add first legend
plt.gca().add_artist(l)
##########################################################################################
# make map inset
##########################################################################################

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
axins = zoomed_inset_axes(ax, 0.035, loc=1)

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
# add colourbar
##########################################################################################

#cb = plt.colorbar()# ticks=ticks,

ticks = arange(0, ncols+1)
years = float(minyear) + ticks*5
labels = [str('%0.0f' % x) for x in years]
labels[-1] += '+'

# normalise
norm = mpl.colors.Normalize(vmin=minyear, vmax=maxyear)

cax = fig.add_axes([0.315,0.035,0.4,0.025]) # setup colorbar axes.
cb = colorbar.ColorbarBase(cax, cmap=cmap, orientation='horizontal', alpha=0.8, norm=norm) # 
cb.set_ticks(years)
cb.set_ticklabels(labels)
cb.set_label('Year of Earthquake', rotation=0, fontsize=15, labelpad=5) # 

##########################################################################################
# finish
##########################################################################################

#plt.savefig('pre_swan_eq_map.png',fmt='png',dpi=300,bbox_inches='tight')
plt.savefig('pre_swan_eq_map.eps',fmt='eps',dpi=300,bbox_inches='tight')
plt.show()