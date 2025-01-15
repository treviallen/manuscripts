#matplotlib quiver example

from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from matplotlib.colors import LightSource, ListedColormap
from numpy import arange, mean, percentile, array, unique, where, argsort, floor, sqrt, delete, argsort, log10
from netCDF4 import Dataset as NetCDFFile
from gmt_tools import cpt2colormap, remove_last_cmap_colour, remove_first_cmap_colour
from os import path, walk, system
from shapely.geometry import Polygon, Point
import pickle
#from obspy.imaging.beachball import Beach
#from hmtk.parsers.catalogue.csv_catalogue_parser import CsvCatalogueParser, CsvCatalogueWriter
from misc_tools import remove_last_cmap_colour, listdir_extension
from mapping_tools import drawshapepoly, get_field_data, drawoneshapepoly, distance, reckon, map_fault_dip_dirn
import shapefile

mpl.style.use('classic')

##########################################################################################
#108/152/-44/-8
urcrnrlat = -1.0
urcrnrlat = -6.0
llcrnrlat = -45.
urcrnrlon = 153.
llcrnrlon = 106
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(figsize=(18,10))
plt.tick_params(labelsize=13)
ax = fig.add_subplot(111)

m = Basemap(projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='l',area_thresh=500.) #i

# draw coastlines, state and country boundaries, edge of map.
#m.shadedrelief()
#m.bluemarble()
#m.etopo()
m.drawcoastlines()
m.drawstates()
m.drawcountries()
#m.fillcontinents(color='w',lake_color='w')
#m.drawmapboundary(fill_color='w')
m.drawparallels(arange(-90.,90.,4.), labels=[1,0,0,0],fontsize=13, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,6.), labels=[0,0,0,1], fontsize=13, dashes=[2, 2], color='0.5', linewidth=0.75)
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
ls = LightSource(azdeg = 135, altdeg = 10)
#norm = mpl.colors.Normalize(vmin=-8000/zscale, vmax=5000/zscale)#myb
#norm = mpl.colors.Normalize(vmin=-1000/zscale, vmax=1900/zscale)#wiki
norm = mpl.colors.Normalize(vmin=-2000/zscale, vmax=3500/zscale)#wiki

rgb = ls.shade(topodat.data, cmap=cmap, norm=norm)
im = m.imshow(rgb, alpha=1.0)

##########################################################################################
# add domains
##########################################################################################
import matplotlib.patheffects as PathEffects
path_effects=[PathEffects.withStroke(linewidth=2.5, foreground="w")]

cptfile = '//Users//trev//Documents//DATA//GMT//cpt//keshet.cpt' # qual-dark-06.cpt keshet.cpt #aurora.cpt #cosam.cpt #Set1_06.cpt
ncols = 7
cmap, zvals = cpt2colormap(cptfile, ncols+1)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncols)))
cs = delete(cs, 2, 0) # get rid of yellow

shpfile = '../../2019/nac_attenuation/shapefiles/adj_neotectonic_domains.shp'
sf = shapefile.Reader(shpfile)

shapes = sf.shapes()
trts = get_field_data(sf, 'trt', 'str')
utrts = unique(trts)

labels = ['Precambrian Craton', 'Reactivated Proterozoic Crust', 'Extended Continental Crust', \
          'Phanerozoic Accretionary Crust', 'Passive Margin', 'Sprigg Orogen']

for c, utrt, label in zip(cs, utrts, labels):
    labeLegend = True
    for shape, trt in zip(shapes, trts):
        if trt == utrt:
            x = []
            y = []
            
            p = 0
            parts = shape.parts
            parts.append(len(shape.points))
            
            for prt in parts[1:]:
                #print('part',prt
                while p < prt:
                    x.append(shape.points[p][0])
                    y.append(shape.points[p][1])
                    p += 1
                
                # close polygon if not polyline
                if x[0] != x[-1] or y[0] != y[-1]:
                    x.append(x[0])
                    y.append(y[0])
            
                # plot each polygon                    
                xx, yy = m(x,y)
                #print(edgecolor)
                if labeLegend == True:
                    plt.fill(xx,yy, facecolor=c, edgecolor='none', linewidth=0.75, alpha=0.2, label=label)
                    labeLegend = False
                else:
                    plt.fill(xx,yy, facecolor=c, edgecolor='none', linewidth=0.75, alpha=0.2) 

l1 = plt.legend(fontsize=11, loc=3, bbox_to_anchor=(0.125, 0))
##########################################################################################
# plot faults
##########################################################################################


##########################################################################################
# annotate text
##########################################################################################
import matplotlib.patheffects as PathEffects
path_effects=[PathEffects.withStroke(linewidth=3, foreground="w")]
"""
# label polygons
x, y = m(135, -18)
plt.text(x, y, 'North Australian\nCraton', size=11, c='maroon', ha='center', va='top', weight='normal', style='italic', path_effects=path_effects)

x, y = m(126, -14.5)
plt.text(x, y, 'Kimberly\nCraton', size=11, c='maroon', ha='center', va='top', weight='normal', style='italic', path_effects=path_effects)

x, y = m(118.75, -21.5)
plt.text(x, y, 'Pilbara\nCraton', size=11, c='maroon', ha='center', va='center', weight='normal', style='italic', path_effects=path_effects)

x, y = m(134.75, -32)
plt.text(x, y, 'Gawler\nCraton', size=11, c='maroon', ha='center', va='center', weight='normal', style='italic', path_effects=path_effects)

x, y = m(119.5, -29.5)
plt.text(x, y, 'Yilgarn\nCraton', size=11, c='maroon', ha='center', va='center', weight='normal', style='italic', path_effects=path_effects)

x, y = m(126, -12.5)
plt.text(x, y, 'Western\nExtended Crust', size=11, c='g', ha='center', weight='normal', style='italic', path_effects=path_effects)

x, y = m(148, -39.5)
plt.text(x, y, 'Eastern\nExtended Crust', size=11, c='g', ha='center', va='center', weight='normal', style='italic', path_effects=path_effects)

x, y = m(130, -24.5)
plt.text(x, y, 'Central Australia Reactivated\nProterozoic Crust', size=11, c='orangered', \
         va='center', ha='center', weight='normal', style='italic', path_effects=path_effects)

x, y = m(146.5, -28)
plt.text(x, y, 'Tasmanides', size=11, c='b', ha='center', weight='normal', style='italic', path_effects=path_effects)

x, y = m(129, -5)
plt.text(x, y, 'Banda Sea', size=11, c='k', ha='center', va='center', weight='normal', style='italic', path_effects=path_effects)

# add faults

x, y = m(129, -9.8)
plt.text(x, y, 'Timor Trough', size=10, c='k', va='center', ha='center', weight='normal', \
         rotation=10 , path_effects=path_effects)
         
x, y = m(117, -11.9)
plt.text(x, y, 'Sunda Arc', size=10, c='k', va='center', ha='center', weight='normal', \
         rotation=-0.5 , path_effects=path_effects)
         
'''x, y = m(122, -7)
plt.text(x, y, 'Flores-Wetar\nThrust', size=10, c='k', va='center', ha='center', weight='normal', \
         rotation=-0 , path_effects=path_effects)
'''         
x, y = m(141, -4.)
plt.text(x, y, 'New Guinea\nHighlands', size=10, c='k', va='center', ha='center', weight='normal', \
         rotation=-0 , path_effects=path_effects)
"""         
#plt.legend(loc=3, fontsize=12)
##########################################################################################
# set colours
##########################################################################################

# get stress range
logstressbins = arange(-.7, 1.6, 0.2)
logstressrng = logstressbins[-1] - logstressbins[0]

ncols = len(logstressbins) - 1

#logstressrng = float(round(maxstress - minstress))

cptfile = '/Users/trev/Documents/DATA/GMT/cpt/temperature.cpt'
cptfile = '//Users//trev//Documents//DATA//GMT//cpt//keshet.cpt'
#cptfile = '//Users//trev//Documents//DATA//GMT//cpt//plasma.cpt'

cmap, zvals = cpt2colormap(cptfile, ncols+1, rev=True)
cmap = remove_first_cmap_colour(cmap)

cs = (cmap(arange(ncols)))
'''
cmap = plt.get_cmap('viridis_r', len(logstressbins))
cs = (cmap(arange(len(logstressbins))))
'''
##########################################################################################
# add earthquakes
##########################################################################################

csvfile = 'brune_stats.csv'
#data = loadtxt(csvfile, delimiter=',', skiprows=1)

print(csvfile)
lines = open(csvfile).readlines()[1:]
lats = []
lons = []
mags = []
qual = []
stressdrops = []

for line in lines:
    dat = line.strip().split(',')
    lons.append(float(dat[2]))
    lats.append(float(dat[3]))
    mags.append(float(dat[8]))
    qual.append(float(dat[-1]))
    stressdrops.append(float(dat[9]))

logstress = log10(array(stressdrops))
        
# get zorder for plotting
sortidx = argsort(argsort(mags))
for i in range(0, len(mags)): #[0:100])):
    if qual[i] == 1:
        #get colour idx
        colidx = int(floor((ncols) * (logstress[i]-logstressbins[0]) / logstressrng))
        if colidx < 0:
            colidx = 0
        if colidx > ncols-1:
            colidx = ncols-1

        x, y = m(lons[i], lats[i])
        zo = sortidx[i] + 20
        plt.plot(x, y, 'o', mfc=list(cs[colidx][0:3]), markeredgecolor='k', markeredgewidth=0.25, \
                 markersize=(5. * mags[i] - 12), zorder=zo, alpha=0.8)
  
# make legend
legmag = [3.5, 4.5, 5.5, 6.5]
legh = []
for lm in legmag:
    x, y = m(0, 0)
    h = m.plot(x, y, 'o', mfc='k', mec='w', mew=0.25, markersize=(5 * lm - 12), alpha=1., zorder=len(mags)+1)
    legh.append(h[0])

legtxt = ('$\mathregular{M_W}$ 3.5', '$\mathregular{M_W}$ 4.5', '$\mathregular{M_W}$ 5.5', '$\mathregular{M_W}$ 6.5')
l = plt.legend(legh, legtxt, loc=3, numpoints=1, fontsize=10, title="Magnitude", labelspacing=0.75)
l.set_zorder(len(mags)+5)

# re-add tectonic reg leg
plt.gca().add_artist(l1)


##########################################################################################
# add colourbar
##########################################################################################
from matplotlib import colorbar

plt.gcf().subplots_adjust(bottom=0.1)
cax = fig.add_axes([0.315,0.035,0.4,0.025]) # setup colorbar axes.

# normalise
norm = mpl.colors.Normalize(vmin=logstressbins[0], vmax=logstressbins[-1])

cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal', alpha=0.9)

ticks = arange(0, ncols+1)
#stressdrops = float(minstress) + ticks*5
labels = [str('%0.1f' % 10**x) for x in logstressbins]

#ticks = arange(0,len(logstressbins))
cb.set_ticks(logstressbins)
cb.set_ticklabels(labels)

cb.set_label('Brune Stress Drop (MPa)', rotation=0, fontsize=15, labelpad=5) # 

##########################################################################################
# finish
##########################################################################################

plt.savefig('mapped_brune_stress_drops.png',fmt='png',dpi=300,bbox_inches='tight')
#plt.savefig('figures/fig_1.eps',fmt='pdf',dpi=300,bbox_inches='tight')
plt.show()