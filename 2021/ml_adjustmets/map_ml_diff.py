#from matplotlib.mlab import griddata
from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource
from numpy import arange, mean, percentile, array, unique, where, argsort, floor, loadtxt
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap
from os import path, walk, system
#from obspy.imaging.beachball import Beach
#from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser, CsvCatalogueWriter
from misc_tools import remove_last_cmap_colour
from os import getcwd

plt.rcParams['pdf.fonttype'] = 42
mpl.style.use('classic')

##########################################################################################
# parse mag zones
##########################################################################################
if getcwd().startswith('/nas'):
    mlRegFile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/catalogue/magnitude/ml/australia_ml_regions.txt'
else:
    mlRegFile = '/Users/trev/Documents/Geoscience_Australia/NSHA2018/catalogue/magnitude/ml/australia_ml_regions.txt'
lines = open(mlRegFile).readlines()

# get WA polygons
walon = []
walat = []
for line in lines:
    if line.startswith('SWAustralia'):
        
        walon.append(float(line.split()[3]))
        walat.append(float(line.split()[2]))

walat = array(walat)
walon = array(walon)

# get SEA polygons
ealon = []
ealat = []
for line in lines:
    if line.startswith('SEAustralia'):
        
        ealon.append(float(line.split()[3]))
        ealat.append(float(line.split()[2]))
        
ealat = array(ealat)
ealon = array(ealon)

# get SA polygons
salon = []
salat = []
for line in lines:
    if line.startswith('SAustralia'):
        
        salon.append(float(line.split()[3]))
        salat.append(float(line.split()[2]))
        
salat = array(salat)
salon = array(salon)


##########################################################################################
# parse epicentres
##########################################################################################
if getcwd().startswith('/nas'):
    magdiffcsv = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/catalogue/matlab/rev_mag_diff.csv'
else:
    magdiffcsv = '/Users/trev/Documents/Geoscience_Australia/NSHA2018/catalogue/matlab/rev_mag_diff.csv'

data = loadtxt(magdiffcsv, delimiter=',')

evlo = data[:,0]
evla = data[:,1]
mlrv = data[:,3]
mldf = data[:,4]

##########################################################################################
#108/152/-44/-8
urcrnrlat = -8.
llcrnrlat = -44.5
urcrnrlon = 152.5
llcrnrlon = 106.5
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(8,4.5))
plt.tick_params(labelsize=8)
ax = fig.add_subplot(111)

m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='l',area_thresh=2000.)

# draw coastlines, state and country boundaries, edge of map.
#m.shadedrelief()
m.drawcoastlines(linewidth=0.3)
m.drawstates(linewidth=0.3)
m.drawcountries(linewidth=0.3)
m.drawparallels(arange(-90.,90.,6.), labels=[1,0,0,0],fontsize=8, dashes=[2, 2], color='0.5', linewidth=0.5)
m.drawmeridians(arange(0.,360.,10.), labels=[0,0,0,1], fontsize=8, dashes=[2, 2], color='0.5', linewidth=0.5)
#m.drawmapscale(144, -34.8, 146., -38.5, 400, fontsize = 16, barstyle='fancy', zorder=100)

##########################################################################################
# plot gebco
##########################################################################################

print( 'Reading netCDF file...')
try:
    nc = NetCDFFile('//Users//trev//Documents//DATA//GMT//GEBCO//Australia_30c.nc')
except:
    nc = NetCDFFile('/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/GEBCO/au_gebco.nc')

zscale =20. #gray
zscale =20. #colour
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
    cmap, zvals = cpt2colormap(cptfile, 30)
except:
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/wiki-2.0.cpt'
    cmap, zvals = cpt2colormap(cptfile, 30)
    
cmap = remove_last_cmap_colour(cmap)
        
# make shading
print( 'Making map...')
ls = LightSource(azdeg = 180, altdeg = 5)
#norm = mpl.colors.Normalize(vmin=-8000/zscale, vmax=5000/zscale)#myb
norm = mpl.colors.Normalize(vmin=-1000/zscale, vmax=1900/zscale)#wiki

rgb = ls.shade(topodat.data, cmap=cmap, norm=norm)
im = m.imshow(rgb, alpha=1.0)


##########################################################################################
# plt ML boundaries
##########################################################################################

x, y = m(walon, walat)
plt.plot(x, y, '-', c='k',lw=.75)
#labelCentroid(plt, m, 'WA', 16, walon, walat, 0)

x, y = m(ealon, ealat)
plt.plot(x, y, '-', c='k',lw=.75)
#labelCentroid(plt, m, 'EA', 16, ealon, ealat, 0)

x, y = m(salon, salat)
plt.plot(x, y, '-', c='k', lw=.75)
#labelCentroid(plt, m, 'SA', 16, salon, salat, 1)

##########################################################################################
# add epicentres
##########################################################################################

ncols = 16
if getcwd().startswith('/nas'):
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/BlueWhiteOrangeRed.cpt'
else:
    cptfile = '/Users/trev/Documents/DATA/GMT/cpt/BlueWhiteOrangeRed.cpt'

cmap, zvals = cpt2colormap(cptfile, ncols+1, rev=True)
cmap = remove_last_cmap_colour(cmap)

cs = (cmap(arange(ncols)))
#cs = vstack((cs[0:1], cs[2:], [1., 1., 1., 1.]))

# get zorder for plotting
sortidx = argsort(argsort(mlrv))
for i in range(0, len(mlrv)): 
    if mlrv[i] >= 0.0:
        #get colour idx
        colidx = int(round((ncols-1) * (mldf[i]+1) / 2.0))
        if colidx < 0:
        	colidx = 0
        x, y = m(evlo[i], evla[i])
        zo = sortidx[i] + 20
        plt.plot(x, y, 'o', markerfacecolor=list(cs[colidx][0:3]), markeredgecolor='k', markeredgewidth=0.25, \
                 markersize=-2 + mlrv[i]*1.5, zorder=zo, alpha=0.8)
    
# make legend
legmag = [3., 5., 7.]
legh = []
for lm in legmag:
    x, y = m(0, 0)
    h = m.plot(x, y, 'ko', mfc='k', markersize=(-2 + lm*1.5), alpha=1., markeredgewidth=0.25, zorder=len(mlrv)+1, lw=2)
    legh.append(h[0])

l = plt.legend(legh, ('$\mathregular{M_L}$ 3.0', '$\mathregular{M_L}$ 5.0', '$\mathregular{M_L}$ 7.0'), loc=3, numpoints=1, fontsize=8)
l.set_zorder(len(mlrv)+5)

##########################################################################################
# label states
##########################################################################################
'''
state = ['WA', 'NT', 'SA', 'QLD', 'NSW', 'VIC', 'TAS']
slat = [-26, -21.0, -30.0, -23.0, -32.5, -37.1, -42.]
slon = [122, 133.5, 135.0, 144.5, 146.5, 143.6, 147.0]
for i, st in enumerate(state):
    x, y = m(slon[i], slat[i])
    plt.text(x, y, st, size=6, horizontalalignment='center', verticalalignment='center', weight='normal')
'''    
##########################################################################################
# add colourbar
##########################################################################################

#cb = plt.colorbar()# ticks=ticks,
ticks = arange(0, 10, 1)
mags = -1 + ticks*0.25
labels = [str('%0.2f' % x) for x in mags]

# normalise
norm = mpl.colors.Normalize(vmin=-1, vmax=1)

cax = fig.add_axes([0.3,0.01,0.4,0.03]) # setup colorbar axes.
cb = colorbar.ColorbarBase(cax, cmap=cmap, orientation='horizontal', alpha=0.8, norm=norm) # 
cb.set_ticks(mags)
cb.set_ticklabels(labels)
cb.ax.tick_params(labelsize=8)
cb.set_label('$\mathregular{M_L}$ Difference', rotation=0, labelpad=2, fontsize=10)

plt.savefig('aus_ml_diff.png', format='png', bbox_inches='tight', dpi=200)
plt.savefig('aus_ml_diff.eps', format='eps', bbox_inches='tight', dpi=200)
plt.show()
