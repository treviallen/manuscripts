from matplotlib.mlab import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource
from numpy import arange, mean, percentile, array, unique, where, argsort, floor, ceil
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from os import path, walk, system
#from obspy.imaging.beachball import Beach
from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser, CsvCatalogueWriter
from misc_tools import remove_last_cmap_colour

plt.rcParams['pdf.fonttype'] = 42
mpl.style.use('classic')

##########################################################################################
# parse epicentres
##########################################################################################

# parse HMTK csv
hmtk_csv = '/nas/gemd/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/catalogue/data/merged_NSHA18-ISCGEM_hmtk.csv'
parser = CsvCatalogueParser(hmtk_csv)    
fullcat = parser.read_file()

hmtk_csv = '/nas/gemd/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/catalogue/data/NSHA18CAT_V0.1_hmtk_declustered.csv'
parser = CsvCatalogueParser(hmtk_csv)    
declcat = parser.read_file()

##########################################################################################
#108/152/-44/-8
urcrnrlat = -29.5
llcrnrlat = -33.3
urcrnrlon = 118.75
llcrnrlon = 115.25

lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(18,10))
plt.tick_params(labelsize=16)
ax = fig.add_subplot(111)

m = Basemap(projection='merc',\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='h',area_thresh=500.)

# lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,

# draw coastlines, state and country boundaries, edge of map.
#m.shadedrelief()
m.drawcoastlines()
m.drawstates()
m.drawcountries()
m.drawparallels(arange(-90.,90.,1.), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,1.), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
#m.drawmapscale(144, -34.8, 146., -38.5, 400, fontsize = 16, barstyle='fancy', zorder=100)

##########################################################################################
# plot gebco
##########################################################################################

print 'Reading netCDF file...'
#nc = NetCDFFile('//Users//tallen//Documents//DATA//GMT//GEBCO//Australia_30c.nc')
nc = NetCDFFile('/nas/gemd/ehp/georisk_earthquake/hazard/DATA/GEBCO/au_gebco.nc')

zscale =20. #gray
zscale =50. #colour
data = nc.variables['elevation'][:] / zscale
lons = nc.variables['lon'][:]
lats = nc.variables['lat'][:]

# transform to metres
nx = int((m.xmax-m.xmin)/500.)+1
ny = int((m.ymax-m.ymin)/500.)+1

topodat = m.transform_scalar(data,lons,lats,nx,ny)

print 'Getting colormap...'
# get colormap
try:
    cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//wiki-2.0.cpt'
    cmap, zvals = cpt2colormap(cptfile, 30)
except:
    cptfile = '/nas/gemd/ehp/georisk_earthquake/hazard/DATA/cpt/wiki-2.0.cpt'
    cmap, zvals = cpt2colormap(cptfile, 30)
    
cmap = remove_last_cmap_colour(cmap)

# make shading
print 'Making map...'
ls = LightSource(azdeg = 180, altdeg = 5)
#norm = mpl.colors.Normalize(vmin=-8000/zscale, vmax=5000/zscale)#myb
norm = mpl.colors.Normalize(vmin=-1000/zscale, vmax=1900/zscale)#wiki
rgb = ls.shade(topodat, cmap=cmap, norm=norm)
im = m.imshow(rgb)


##########################################################################################
# add epicentres
##########################################################################################

ncols = 18

# get year range
minyear = 10*floor(fullcat.data['year'][0]/10.)
maxyear = fullcat.data['year'][-1]
yearrng = float(round(maxyear - minyear))
#cmap = plt.cm.get_cmap('Spectral', ncols)

try:
    cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/temperature.cpt'
    cmap, zvals = cpt2colormap(cptfile, ncols+1, rev=False)
except:
    cptfile = '/nas/gemd/ehp/georisk_earthquake/hazard/DATA/cpt/temperature.cpt'
    cmap, zvals = cpt2colormap(cptfile, ncols+1, rev=False)
    
cmap = remove_last_cmap_colour(cmap)

cs = (cmap(arange(ncols)))

# get zorder for plotting - plt full cat
sortidx = argsort(argsort(fullcat.data['magnitude']))
for i in range(0, len(fullcat.data['magnitude'])): #[0:100])):
    if fullcat.data['magnitude'][i] >= 2.5:
       x, y = m(fullcat.data['longitude'][i], fullcat.data['latitude'][i])
       zo = sortidx[i] + 20
       plt.plot(x, y, 'o', mfc='limegreen', markeredgecolor='k', markeredgewidth=0.5, \
                 markersize=-4 + fullcat.data['magnitude'][i]*3.2, zorder=zo)
        

# get zorder for plotting - plt dec cat
sortidx = argsort(argsort(declcat.data['magnitude']))
for i in range(0, len(declcat.data['magnitude'])): #[0:100])):
    if declcat.data['magnitude'][i] >= 2.5:
       x, y = m(declcat.data['longitude'][i], declcat.data['latitude'][i])
       zo = sortidx[i] + 20 + len(fullcat.data['magnitude'])
       plt.plot(x, y, 'o', mfc='orangered', markeredgecolor='k', markeredgewidth=0.5, \
                 markersize=-4 + declcat.data['magnitude'][i]*3.2, zorder=zo)

    
# make 1st legend
legmag = [3., 5., 7.]
legh1 = []
for lm in legmag:
    x, y = m(0, 0)
    h = m.plot(x, y, 'o', markeredgecolor='k', mfc='k', markersize=(-4 + lm*3.2), alpha=1., zorder=10*len(fullcat.data['magnitude'])+1, lw=2)
    legh1.append(h[0])

l1 = plt.legend(legh1, ('MW 3.0', 'MW 5.0', 'MW 7.0'), loc=2, numpoints=1, fontsize=13)
l1.set_zorder(10*len(fullcat.data['magnitude'])+1)

# make 2nd legend
legh2 = []
x, y = m(0, 0)
lm = 5.5
h = m.plot(x, y, 'o', markeredgecolor='k', mfc='limegreen', markersize=(-4 + lm*3.2), alpha=1., zorder=10*len(fullcat.data['magnitude'])+2, lw=2)
legh2.append(h[0])
h = m.plot(x, y, 'o', markeredgecolor='k', mfc='orangered', markersize=(-4 + lm*3.2), alpha=1., zorder=10*len(fullcat.data['magnitude'])+2, lw=2)
legh2.append(h[0])


l2 = plt.legend(legh2, ('Full Catalogue', 'Declustered'), loc=1, numpoints=1, fontsize=13)
l2.set_zorder(10*len(fullcat.data['magnitude'])+2)

# re-add l1
plt.gca().add_artist(l1)

##########################################################################################
# add cities
##########################################################################################

clat = [-31.95]
clon = [115.86]
cloc = ['Perth']
x, y = m(clon, clat)
plt.plot(x, y, 'o', markerfacecolor='maroon', markeredgecolor='w', markeredgewidth=2, markersize=18)

# label cities
off = 0.09
for i, loc in enumerate(cloc):
    if i == 1:
        x, y = m(clon[i]+0.08, clat[i]+0.08)
        plt.text(x, y, loc, size=20, ha='left', weight='normal')
    else:
        x, y = m(clon[i]-0.08, clat[i]+0.08)
        plt.text(x, y, loc, size=20, ha='right', weight='normal')


##########################################################################################
# make map inset
##########################################################################################

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes

axins = zoomed_inset_axes(ax, 0.022, loc=3)

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
plt.plot(x, y, 'r-', lw=1.5)




plt.savefig('map_swsz_decluster.png', format='png', bbox_inches='tight', dpi=200)
plt.show()
