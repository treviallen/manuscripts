#matplotlib quiver example

from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource, ListedColormap
from numpy import arange, mean, percentile, array, unique, where, argsort, floor, sqrt, delete, argsort
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from data_fmt_tools import parse_iris_stationlist
from mapping_tools import annotate_cities
from io_catalogues import parse_ga_event_query
from os import path, walk, system
from shapely.geometry import Polygon, Point
import pickle
#from obspy.imaging.beachball import Beach
#from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser, CsvCatalogueWriter
from misc_tools import remove_last_cmap_colour, listdir_extension, remove_first_cmap_colour, dictlist2array
from mapping_tools import drawshapepoly, get_field_data, drawoneshapepoly, distance, reckon, map_fault_dip_dirn, make_street_map
import shapefile
#from io_catalogues import parse_ga_event_query
from tools.mfd_tools import parse_hmtk_cat
from obspy import UTCDateTime

mpl.style.use('classic')

##########################################################################################
#108/152/-44/-8
eqlo = 150.85
eqla = -32.34
urcrnrlat = eqla + 1.25
llcrnrlat = eqla - 1.25
urcrnrlon = eqlo + 1.25
llcrnrlon = eqlo - 1.25
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(18,10))
plt.tick_params(labelsize=13)
ax = fig.add_subplot(111)
'''
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
m.drawparallels(arange(-90.,90.,1.),fontsize=13, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,1.), fontsize=13, dashes=[2, 2], color='0.5', linewidth=0.75)
'''
#m.drawmapscale(144, -34.8, 146., -38.5, 50, fontsize = 12, barstyle='fancy', zorder=100)

##########################################################################################
# set up street map
##########################################################################################

# set map centroid
clon = mean([llcrnrlon, urcrnrlon])
clat = mean([llcrnrlat, urcrnrlat])
            
degrng = urcrnrlon-llcrnrlon
ll_buffer = (urcrnrlon - llcrnrlon) / 2.
'''
plt, m, ax = make_street_map(clat, clon, service='World_Imagery', \
             xpixels = 150, plt_inset = False, plt_marker = False, inset_loc=4, inset_multiplier=0.1)
'''
from mpl_toolkits.basemap import Basemap

m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat, epsg=3112)
m.drawparallels(arange(-90.,90.,0.5), labels=[1,0,0,0],fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.0)
m.drawmeridians(arange(0.,360.,0.5), labels=[0,0,0,1], fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.0)
#http://server.arcgisonline.com/arcgis/rest/services

m.arcgisimage(service='World_Imagery', xpixels = 1500, verbose= True)

'''
    Map Services:
        ESRI_Imagery_World_2D
        ESRI_StreetMap_World_2D
        NatGeo_World_Map
        NGS_Topo_US_2D
        Ocean_Basemap
        USA_Topo_Maps
        World_Imagery
        World_Physical_Map
        World_Shaded_Relief
        World_Street_Map
        World_Terrain_Base
        World_Topo_Map
        
ESRI_Imagery_World_2D - good sat image
NatGeo_World_Map - nice road map - hard to read
'''

##########################################################################################
# add earthquakes
##########################################################################################
import matplotlib.patheffects as PathEffects
#path_effects=[PathEffects.withStroke(linewidth=2.5, foreground="w")]

##########################################################################################
# set colours
##########################################################################################

ncols = 15

# get year range
minyear = 1880
maxyear = 2030
yearrng = float(round(maxyear - minyear))

#cptfile = '/Users/trev/Documents/DATA/GMT/cpt/temperature.cpt'
cptfile = '//Users//trev//Documents//DATA//GMT//cpt//keshet.cpt'
cmap, zvals = cpt2colormap(cptfile, ncols+1, rev=True)
cmap = remove_first_cmap_colour(cmap)

cs = (cmap(arange(ncols)))

##########################################################################################
# add earthquakes
##########################################################################################
'''
hmtk_csv = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/catalogue/data/merged_NSHA23-ISCGEM_hmtk.csv'

evdict = []
lines = open(hmtk_csv).readlines()[1:]
	
for line in lines:
    dat = line.strip().split(',')
    tmp = {'lon':float(dat[7]), 'lat':float(dat[8]), 'mag':float(dat[10]), 'year':float(dat[1])}
    	
    evdict.append(tmp) 
'''
evdict = parse_ga_event_query('mswl_maps_earthquakes_export.csv')

lons = dictlist2array(evdict, 'lon')
lats = dictlist2array(evdict, 'lat')
mags = dictlist2array(evdict, 'mag')
years = dictlist2array(evdict, 'year')

# get zorder for plotting
sortidx = argsort(argsort(mags))
for i in range(0, len(mags)): #[0:100])):
    #get colour idx
    if mags[i] >= 2.:
        colidx = int(floor((ncols) * (years[i]-minyear+0.1) / yearrng))
        x, y = m(lons[i], lats[i])
        zo = sortidx[i] + 20
        plt.plot(x, y, 'o', mfc=list(cs[colidx][0:3]), markeredgecolor='k', markeredgewidth=0.25, \
                 markersize=(2.75 * mags[i] - 1), zorder=zo, alpha=0.8)

# make legend
legmag = [2., 3.5, 5.]
legh = []
for lm in legmag:
    x, y = m(0, 0)
    h = m.plot(x, y, 'o', mfc='w', mec='k', mew=0.25, markersize=(2.75 * lm - 1), alpha=1., zorder=len(mags)+1)
    legh.append(h[0])

legtxt = ('$\mathregular{M}$ 2.0', '$\mathregular{M}$ 3.5', '$\mathregular{M}$ 5.0')
l = plt.legend(legh, legtxt, loc=1, numpoints=1, fontsize=12, title="Magnitude", title_fontsize=16, \
               labelspacing=0.75, frameon=False)
for text in l.get_texts():
    text.set_color("w")
l.set_zorder(len(mags)+5)
'''
import matplotlib.font_manager as font_manager
font = font_manager.FontProperties(weight='bold',
                                   style='normal')
'''
plt.setp(l.get_title(), color='white')

##########################################################################################
# add 50 km buffer
##########################################################################################

from mapping_tools import distance, reckon

azims = arange(0,360,1)
rng = 15

azim_lon = []
azim_lat = []
eqlo = 150.86
eqla = -32.34
for azim in azims:
    loc = reckon(eqla, eqlo, rng, azim)
    azim_lon.append(loc[0])
    azim_lat.append(loc[1])
    
azim_lon = array(azim_lon)
azim_lat = array(azim_lat)

x,y = m(azim_lon, azim_lat)
plt.plot(x, y, '--', c='c', lw=1.5) #, zorder=41000)

##########################################################################################
# plt stations
##########################################################################################
import matplotlib.patheffects as PathEffects
path_effects=[PathEffects.withStroke(linewidth=2, foreground="w")]

stationlist = '/Users/trev/Documents/Networks/AU/gmap-stations-edit.txt'
staDict = parse_iris_stationlist(stationlist)

# plot stations
sta_lon = dictlist2array(staDict, 'lon')
sta_lat = dictlist2array(staDict, 'lat')
sta_code = dictlist2array(staDict, 'sta')

cmdnet = ('NTLS', 'NTLH', 'MGCD')
cmdlat = []
cmdlon = []
for stla, stlo, sta in zip(sta_lat, sta_lon, sta_code):
    for cn in cmdnet:
        if sta == cn:
            print(cn)
            cmdlat.append(stla)
            cmdlon.append(stlo)
            print(cmdlat)
            if sta == 'NTLS':
                x,y = m(stlo-0.015, stla-0.024)
                plt.text(x, y, sta, size=11, c='k', va='top', ha='right', weight='normal', style='italic', \
                     path_effects=path_effects, zorder=40000)
            else:
                x,y = m(stlo-0.015, stla-0.005)
                plt.text(x, y, sta, size=11, c='k', va='bottom', ha='right', weight='normal', style='italic', \
                     path_effects=path_effects, zorder=40000)
            
'''
x,y = m(slon+0.004, slat+0.0015)
            plt.text(x, y, sta, size=12, c='k', va='bottom', ha='left', weight='normal', style='italic', \
                     path_effects=path_effects, zorder=11000)
'''
x,y = m(sta_lon, sta_lat)
plt.plot(x, y, '^', c='lime', ms=10, zorder=41000)
'''
x,y = m(cmdlon, cmdlat)
plt.plot(x, y, '^', c='yellow', ms=12, zorder=1000)
'''
# label stas
urcrnrlat
llcrnrlat
urcrnrlon
llcrnrlon

##########################################################################################
# add cities
##########################################################################################
numCities = 15
blacklist = ['Umina', 'Bateau Bay']
annotate_cities(numCities, plt, m, markersize=8, markerfacecolor='k', markeredgecolor='w', \
                markeredgewidth=1., fs=12, weight='medium', marker='s', blacklist=blacklist, splittext=True)

##########################################################################################
# add colourbar
##########################################################################################

#cb = plt.colorbar()# ticks=ticks,

ticks = arange(0, (ncols+1)/2)
years = float(minyear) + ticks*20
labels = [str('%0.0f' % x) for x in years]

# normalise
norm = mpl.colors.Normalize(vmin=minyear, vmax=maxyear)

cax = fig.add_axes([0.33,0.045,0.34,0.025]) # setup colorbar axes.
cb = colorbar.ColorbarBase(cax, cmap=cmap, orientation='horizontal', alpha=0.8, norm=norm) # 
cb.set_ticks(years)
cb.set_ticklabels(labels)
cb.ax.tick_params(color='k', labelcolor='k')
for spine in cb.ax.spines.values():
    spine.set_edgecolor('k')

cb.set_label('Year of Earthquake', rotation=0, fontsize=15, labelpad=5, color='k')
cb.outline.set_edgecolor('k') # 

##########################################################################################
# make map inset
##########################################################################################

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

axins = zoomed_inset_axes(ax, 0.012, loc=2)

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
xv = mean([llcrnrlon, urcrnrlon])
yv = mean([llcrnrlat, urcrnrlat])
x, y = m2(xv, yv)
plt.plot(x, y, 'rs',ms=6)

##########################################################################################
# finish
##########################################################################################

plt.savefig('hunter_eq_map_satellite.png',fmt='png',dpi=300,bbox_inches='tight')
#plt.savefig('figures/fig_1.eps',fmt='pdf',dpi=300,bbox_inches='tight')
plt.show()