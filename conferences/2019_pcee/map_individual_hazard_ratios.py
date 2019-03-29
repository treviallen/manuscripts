from matplotlib import colors, colorbar
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from numpy import arange, mean, percentile, array, unique, where, argsort, floor, ceil, ones_like, sin, cos, radians, log10, isnan, linspace
from gmt_tools import cpt2colormap
from misc_tools import remove_last_cmap_colour
import matplotlib.gridspec as gridspec

plt.rcParams['pdf.fonttype'] = 42
mpl.style.use('classic')

##########################################################################################
# parse hazard data
##########################################################################################

#hazFactFile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/NSHM_18/reporting/NSHA18-grids_and_maps/data_pkg/3_locality_data/3.1_hazard_factors/nsha18_hazard_factor_comparison.csv'
hazFactFile = 'nsha18_hazard_factor_comparison_edit.csv'

locLats = []
locLons = []
hazFactAS1170 = []
hazFactNSHA10 = []
hazFactNSHA2 = []

lines = open(hazFactFile).readlines()[1:]

for line in lines:
    data = line.strip().split(',')
    
    locLons.append(float(data[1]))
    locLats.append(float(data[2]))
    hazFactAS1170.append(float(data[3]))
    hazFactNSHA10.append(float(data[6]))
    hazFactNSHA2.append(float(data[-5]))
    
hazFactAS1170 = array(hazFactAS1170)
hazFactNSHA10 = array(hazFactNSHA10)
hazFactNSHA2 = array(hazFactNSHA2)   

# set min AS1170.4 value to 0.08
idx = hazFactAS1170 < 0.08
hazFactAS1170[idx] = 0.08

##########################################################################################
# get cmap
##########################################################################################
ncols = 17
try:
    cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//temp_18lev.cpt'
    cmap, zvals = cpt2colormap(cptfile, ncols)
except:
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/temp_19lev.cpt'
    cmap, zvals = cpt2colormap(cptfile, ncols)
    
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncols-1)))

##########################################################################################
# make first plot AS1170.4-2018 / NSHA 10%
##########################################################################################


llcrnrlat = -44
urcrnrlat = -8.
llcrnrlon = 107.
urcrnrlon = 152.
                
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)

fig = plt.figure(1, figsize=(15,11))
gs1 = gridspec.GridSpec(1, 2)
hspace = 0.01
gs1.update(wspace=0.11, hspace=hspace) # negative looks bad in "show", but ok in pngs

ax = fig.add_subplot(111)
plt.tick_params(labelsize=8)

m = Basemap(projection='lcc',\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            lat_1=lat_1,lat_2=lat_2,lon_0=lon_0, \
            rsphere=6371200.,resolution='l',area_thresh=2000.)

#lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
# draw coastlines, state and country boundaries, edge of map.
m.shadedrelief()
m.drawstates()
#m.drawcountries()
m.drawparallels(arange(-90.,90.,6.), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,8.), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)

logRatios = log10(hazFactNSHA10 / hazFactAS1170)

minscale = log10(0.25)
maxscale = log10(4.0)

for lat, lon, lograt in zip(locLats, locLons, logRatios):
    if not isnan(lograt):
        colidx = int(round((ncols-2) * (lograt-minscale) / (maxscale-minscale)))
        if colidx < 0:
            colidx = 0
        if colidx > ncols-2:
            colidx = ncols-2
            
        x, y = m(lon, lat)
        plt.plot(x, y, 'o', mfc=list(cs[colidx][0:3]), markeredgecolor='k', markeredgewidth=0.5, \
                 markersize=9, alpha=1.)

plt.title('Ratio of NSHA18 1/500 AEP and\nAS1170.4 Hazard Factors', fontsize=22)
xlim = ax.get_xlim()
xtxt = xlim[1] * 0.02
ylim = ax.get_ylim()
ytxt = ylim[1] * 0.02

'''
###########################################################################################
make colourbar
###########################################################################################
'''    

bounds = linspace(log10(0.25),log10(4),9)
colbounds = linspace(log10(0.25),log10(4),17)
norm = colors.BoundaryNorm(boundaries=colbounds, ncolors=ncols-1)
# set colourbar
plt.gcf().subplots_adjust(bottom=0.1)
cby = 0.025
cbh = 0.035
cax = fig.add_axes([0.3,cby,0.4,cbh]) # setup colorbar axes.
#norm = colors.Normalize(vmin=vmin, vmax=vmax)
cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')

cb.set_ticks(bounds)
labels = [str(round(10**x, 2)) for x in bounds]
labels[0] = '0.0'
cb.set_ticklabels(labels)
cb.ax.tick_params(labelsize=16)

# set title
titlestr = 'Hazard Ratio (NSHA18 / AS1170.4-2017 [R2108])'
cb.set_label(titlestr, fontsize=20)

plt.savefig('map_design_ratios_500.png', format='png', bbox_inches='tight', dpi=200)
#plt.show()

##########################################################################################
# make second plot AS1170.4-2018 / NSHA 2%
##########################################################################################
fig = plt.figure(2, figsize=(15,11))
ax = fig.add_subplot(111)
plt.tick_params(labelsize=8)

m = Basemap(projection='lcc',\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            lat_1=lat_1,lat_2=lat_2,lon_0=lon_0, \
            rsphere=6371200.,resolution='l',area_thresh=2000.)

#lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
# draw coastlines, state and country boundaries, edge of map.
m.shadedrelief()
m.drawstates()
#m.drawcountries()
m.drawparallels(arange(-90.,90.,6.), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,8.), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)

logRatios = log10(hazFactNSHA2 / hazFactAS1170)

minscale = log10(0.25)
maxscale = log10(4.0)

for lat, lon, lograt in zip(locLats, locLons, logRatios):
    if not isnan(lograt):
        colidx = int(round((ncols-2) * (lograt-minscale) / (maxscale-minscale)))
        if colidx < 0:
            colidx = 0
        if colidx > ncols-2:
            colidx = ncols-2
            
        x, y = m(lon, lat)
        plt.plot(x, y, 'o', mfc=list(cs[colidx][0:3]), markeredgecolor='k', markeredgewidth=0.5, \
                 markersize=9, alpha=1.)

plt.title('Ratio of NSHA18 1/2475 AEP and\nAS1170.4 Hazard Factors', fontsize=22)
#plt.text(xtxt, ytxt, 'b)', fontsize=22, va='bottom', ha='left')

'''
###########################################################################################
make colourbar
###########################################################################################
'''    

bounds = linspace(log10(0.25),log10(4),9)
colbounds = linspace(log10(0.25),log10(4),17)
norm = colors.BoundaryNorm(boundaries=colbounds, ncolors=ncols-1)
# set colourbar
plt.gcf().subplots_adjust(bottom=0.1)
cby = 0.025
cbh = 0.035
cax = fig.add_axes([0.3,cby,0.4,cbh]) # setup colorbar axes.
#norm = colors.Normalize(vmin=vmin, vmax=vmax)
cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')

cb.set_ticks(bounds)
labels = [str(round(10**x, 2)) for x in bounds]
labels[0] = '0.0'
cb.set_ticklabels(labels)
cb.ax.tick_params(labelsize=16)

# set title
titlestr = 'Hazard Ratio (NSHA18 / AS1170.4-2017 [R2108])'
cb.set_label(titlestr, fontsize=20)

plt.savefig('map_design_ratios_2475.png', format='png', bbox_inches='tight', dpi=200)
plt.show()
