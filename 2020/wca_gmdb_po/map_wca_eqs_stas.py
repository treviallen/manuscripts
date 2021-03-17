from scipy.interpolate import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource
from numpy import arange, mean, percentile, array, unique, where, argsort, floor
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from io_catalogues import parse_ga_event_query
from misc_tools import remove_last_cmap_colour, dictlist2array
from mapping_tools import drawshapepoly #, labelpolygon, get_field_data
import shapefile

plt.rcParams['pdf.fonttype'] = 42
mpl.style.use('classic')

##########################################################################################
# parse epicentres
##########################################################################################
#hmtk_csv = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/catalogue/data/AUSTCAT_V0.12_hmtk_deblast.csv'
wca_evlist = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/Ground_Motion/Data/western-central_australia/unique_event_list.csv'
# parse  csv

lines = open(wca_evlist).readlines()

mag = []
lon = []
lat = []
nrec = []

for line in lines[1:]:
    dat = line.strip().split(',')
    mag.append(float(dat[6]))
    lon.append(float(dat[1]))
    lat.append(float(dat[2]))
    nrec.append(int(dat[-1]))
    
##########################################################################################
# parse stations
##########################################################################################
#hmtk_csv = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/catalogue/data/AUSTCAT_V0.12_hmtk_deblast.csv'
wca_stalist = 'sta_vs30_list.csv'
# parse  csv

lines = open(wca_stalist).readlines()

sta = []
stalon = []
stalat = []

for line in lines[1:]:
    dat = line.strip().split(',')
    sta.append(dat[0])
    stalon.append(float(dat[1]))
    stalat.append(float(dat[2]))

stalon = array(stalon)
stalat = array(stalat)  
##########################################################################################
#108/152/-44/-8
urcrnrlat = -16.
llcrnrlat = -36.
urcrnrlon = 138.
llcrnrlon = 111.
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(18,10))
plt.tick_params(labelsize=16)
ax = fig.add_subplot(111)

m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='h',area_thresh=2000.)

# draw coastlines, state and country boundaries, edge of map.
#m.shadedrelief()
m.drawcoastlines()
m.drawstates()
m.drawcountries()
m.drawparallels(arange(-90.,90.,4.), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,6.), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
#m.drawmapscale(144, -34.8, 146., -38.5, 400, fontsize = 16, barstyle='fancy', zorder=100)

##########################################################################################
# plot gebco
##########################################################################################

print('Reading netCDF file...')
#nc = NetCDFFile('//Users//trev//Documents//DATA//GMT//GEBCO//au_gebco.grd')
nc = NetCDFFile('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/GEBCO/au_gebco.nc')

zscale =20. #gray
zscale =30. #colour
data = nc.variables['elevation'][:] / zscale
lons = nc.variables['lon'][:]
lats = nc.variables['lat'][:]

# transform to metres
nx = int((m.xmax-m.xmin)/500.)+1
ny = int((m.ymax-m.ymin)/500.)+1

topodat = m.transform_scalar(data,lons,lats,nx,ny)

print('Getting colormap...')
# get colormap
cptfile = '//Users//trev//Documents//DATA//GMT//cpt//wiki-2.0.cpt'
cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/wiki-2.0.cpt'
cmap, zvals = cpt2colormap(cptfile, 30)
cmap = remove_last_cmap_colour(cmap)

# make shading
print('Making map...')
ls = LightSource(azdeg = 180, altdeg = 5)
norm = mpl.colors.Normalize(vmin=-1000/zscale, vmax=1900/zscale)#wiki
rgb = ls.shade(topodat, cmap=cmap, norm=norm)
im = m.imshow(rgb)

##########################################################################################
# add stations
##########################################################################################

x, y = m(stalon, stalat)
m.plot(x, y, 'w^', ms=8)

print 'Get station network!'

##########################################################################################
# add epicentres
##########################################################################################

ncols = 5

# get year range
#minrecs = 10*floor(cat.data['year'][0]/10.)
minrecs = 0
maxrecs = 25 #cat.data['year'][-1]
recsrng = 25
#cmap = plt.cm.get_cmap('Spectral', ncols)

try:
    cptfile = '/Users/tallen/Documents/DATA/GMT/cpt/Sandy_Scars_Aus.cpt'
    cmap, zvals = cpt2colormap(cptfile, ncols+1, rev=False)
except:
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/Sandy_Scars_Aus.cpt'
    cmap, zvals = cpt2colormap(cptfile, ncols+1, rev=False)
    
cmap = remove_last_cmap_colour(cmap)

cs = (cmap(arange(ncols)))
#cs = vstack((cs[0:1], cs[2:], [1., 1., 1., 1.]))

# get zorder for plotting
sortidx = argsort(argsort(mag))
for i in range(0, len(mag)):
    if mag[i] >= 3.5:
        #get colour idx
        n = nrec[i]
        if n > 25:
            n = 25
        colidx = int(round((ncols-1) * n / recsrng))
        x, y = m(lon[i], lat[i])
        zo = sortidx[i] + 20
        plt.plot(x, y, 'o', mfc=list(cs[colidx]), markeredgecolor='k', markeredgewidth=0.5, \
                 markersize=-10 + mag[i]*5, zorder=zo, alpha=0.8)
    
# make legend
legmag = [4., 5., 6.]
legh = []
for lm in legmag:
    x, y = m(0, 0)
    h = m.plot(x, y, 'ko', mfc='k', markersize=(-10 + lm*5), alpha=1., zorder=len(mag)+1, lw=2)
    legh.append(h[0])

l = plt.legend(legh, ('MW 4.0', 'MW 5.0', 'MW 6.0'), loc=2, numpoints=1)
l.set_zorder(len(mag)+5)

##########################################################################################
# add domains
##########################################################################################
shpfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/source_models/zones/shapefiles/Domains/Domains_Sep2011.shp'
sfz = shapefile.Reader(shpfile)
drawshapepoly(m, plt, sfz, edgecolor='r', col='r', lw=1., zorder=1)

##########################################################################################
# add colourbar
##########################################################################################

#cb = plt.colorbar()# ticks=ticks,
ticks = arange(0, 26, 5)
#years = float(minrecs) + ticks*10
labels = [str(x) for x in ticks]
labels[-1] = labels[-1]+'+'

# normalise
norm = mpl.colors.Normalize(vmin=minrecs, vmax=maxrecs)

cax = fig.add_axes([0.79,0.3,0.015,0.4]) # setup colorbar axes.
cb = colorbar.ColorbarBase(cax, cmap=cmap, orientation='vertical', alpha=0.8, norm=norm)
cb.set_ticks(ticks)
cb.set_ticklabels(labels)
cb.set_label('Number of Records', rotation=270, labelpad=20, fontsize=18)

plt.savefig('wca_eqs_sta.png', format='png', bbox_inches='tight', dpi=150)
plt.show()
