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
from mapping_tools import drawshapepoly, labelpolygon
import shapefile
import matplotlib.gridspec as gridspec

plt.rcParams['pdf.fonttype'] = 42
mpl.style.use('classic')

fig = plt.figure(figsize=(24,30))

gs1 = gridspec.GridSpec(3, 2)
gs1.update(wspace=-0.06, hspace=0.05) # negative looks bad in "show", but ok in pngs


##########################################################################################
# parse epicentres
##########################################################################################

# parse HMTK csv
hmtk_csv = '/nas/gemd/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/catalogue/data/NSHA18CAT_V0.1_hmtk_declustered.csv'
hmtk_csv = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/catalogue/data/NSHA18CAT_V0.1_hmtk_declustered.csv'
parser = CsvCatalogueParser(hmtk_csv)    
cat = parser.read_file()

##########################################################################################
#108/152/-44/-8
urcrnrlat = -8.
llcrnrlat = -46.
urcrnrlon = 157.
llcrnrlon = 109.
lon_0 = mean([llcrnrlon, urcrnrlon])
lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
lat_2 = percentile([llcrnrlat, urcrnrlat], 75)


ax = fig.add_subplot(gs1[0])
plt.tick_params(labelsize=8)

m = Basemap(projection='merc',\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='l',area_thresh=2000.)

#lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
# draw coastlines, state and country boundaries, edge of map.
m.shadedrelief()
#m.drawcoastlines()
m.drawstates()
#m.drawcountries()
m.drawparallels(arange(-90.,90.,4.), labels=[0,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,6.), labels=[0,0,0,0], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)

##########################################################################################
# plt earthquakes
##########################################################################################

lats = cat.data['latitude']
lons = cat.data['longitude']
mags = cat.data['magnitude']
year = cat.data['year']

idx = where((mags >= 3.0) & (year >= 1977))[0]
x, y = m(lons[idx], lats[idx])
m.plot(x, y, 'o', ms=2, mfc='0.4', mec='0.4', zorder=1)

idx = where((mags >= 5.6) & (year >= 1880))[0]
x, y = m(lons[idx], lats[idx])
m.plot(x, y, 'o', ms=7, mfc='maroon', mec='w', mew=1, zorder=2)

##########################################################################################
# add shapefiles
##########################################################################################
shpfile = 'ARUP/ARUP_source_model.shp'
sf = shapefile.Reader(shpfile)
drawshapepoly(m, plt, sf,col='r',lw=1.5)

# label shapes
#labelpolygon(m, plt, sf, 'NAME', fweight='normal', fsize=16, addOutline=True)

# add subplot letter
x, y = m(llcrnrlon+1, urcrnrlat-1)
plt.text(x, y, '(a)', size=24, ha='left', va='top', weight='normal')

# add model name
x, y = m(llcrnrlon+1, llcrnrlat+1)
plt.text(x, y, 'Arup', size=36, ha='left', va='bottom', weight='normal')


##########################################################################################
# add 2nd axis
##########################################################################################

ax = fig.add_subplot(gs1[1])
plt.tick_params(labelsize=8)

m = Basemap(projection='merc',\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='l',area_thresh=2000.)

#lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
# draw coastlines, state and country boundaries, edge of map.
m.shadedrelief()
#m.drawcoastlines()
m.drawstates()
#m.drawcountries()
m.drawparallels(arange(-90.,90.,4.), labels=[0,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,6.), labels=[0,0,0,0], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)

##########################################################################################
# plt earthquakes
##########################################################################################

idx = where((mags >= 3.0) & (year >= 1977))[0]
x, y = m(lons[idx], lats[idx])
m.plot(x, y, 'o', ms=2, mfc='0.4', mec='0.4', zorder=1)

idx = where((mags >= 5.6) & (year >= 1880))[0]
x, y = m(lons[idx], lats[idx])
m.plot(x, y, 'o', ms=7, mfc='maroon', mec='w', mew=1, zorder=2)

##########################################################################################
# add shapefiles
##########################################################################################
shpfile = 'ARUP_Background/ARUP_background_source_model.shp'
sf = shapefile.Reader(shpfile)
drawshapepoly(m, plt, sf,col='r',lw=1.5)

# label shapes
#labelpolygon(m, plt, sf, 'NAME', fweight='normal', fsize=16, addOutline=True)

# add subplot letter
x, y = m(llcrnrlon+1, urcrnrlat-1)
plt.text(x, y, '(b)', size=24, ha='left', va='top', weight='normal')

# add model name
x, y = m(llcrnrlon+1, llcrnrlat+1)
plt.text(x, y, 'Arup Background', size=36, ha='left', va='bottom', weight='normal')

##########################################################################################
# add 3rd axis
##########################################################################################

ax = fig.add_subplot(gs1[2])
plt.tick_params(labelsize=8)

m = Basemap(projection='merc',\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='l',area_thresh=2000.)

#lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
# draw coastlines, state and country boundaries, edge of map.
m.shadedrelief()
#m.drawcoastlines()
m.drawstates()
#m.drawcountries()
m.drawparallels(arange(-90.,90.,4.), labels=[0,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,6.), labels=[0,0,0,0], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)

##########################################################################################
# plt earthquakes
##########################################################################################

idx = where((mags >= 3.0) & (year >= 1977))[0]
x, y = m(lons[idx], lats[idx])
m.plot(x, y, 'o', ms=2, mfc='0.4', mec='0.4', zorder=1)

idx = where((mags >= 5.6) & (year >= 1880))[0]
x, y = m(lons[idx], lats[idx])
m.plot(x, y, 'o', ms=7, mfc='maroon', mec='w', mew=1, zorder=2)

##########################################################################################
# add shapefiles
##########################################################################################
shpfile = 'Domains/Domains_Sep2011.shp'
sf = shapefile.Reader(shpfile)
drawshapepoly(m, plt, sf,col='r',lw=1.5)

# label shapes
#labelpolygon(m, plt, sf, 'NAME', fweight='normal', fsize=16, addOutline=True)

# add subplot letter
x, y = m(llcrnrlon+1, urcrnrlat-1)
plt.text(x, y, '(c)', size=24, ha='left', va='top', weight='normal')

# add model name
x, y = m(llcrnrlon+1, llcrnrlat+1)
plt.text(x, y, 'Neotectonic Domains', size=36, ha='left', va='bottom', weight='normal')

##########################################################################################
# add 4th axis
##########################################################################################

ax = fig.add_subplot(gs1[3])
plt.tick_params(labelsize=8)

m = Basemap(projection='merc',\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='l',area_thresh=2000.)

#lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
# draw coastlines, state and country boundaries, edge of map.
m.shadedrelief()
#m.drawcoastlines()
m.drawstates()
#m.drawcountries()
m.drawparallels(arange(-90.,90.,4.), labels=[0,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,6.), labels=[0,0,0,0], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)

##########################################################################################
# plt earthquakes
##########################################################################################

idx = where((mags >= 3.0) & (year >= 1977))[0]
x, y = m(lons[idx], lats[idx])
m.plot(x, y, 'o', ms=2, mfc='0.4', mec='0.4', zorder=1)

idx = where((mags >= 5.6) & (year >= 1880))[0]
x, y = m(lons[idx], lats[idx])
m.plot(x, y, 'o', ms=7, mfc='maroon', mec='w', mew=1, zorder=2)

##########################################################################################
# add shapefiles
##########################################################################################
shpfile = 'Leonard2008/source_model_leonard_2008.shp'
sf = shapefile.Reader(shpfile)
drawshapepoly(m, plt, sf,col='r',lw=1.5)

# label shapes
#labelpolygon(m, plt, sf, 'NAME', fweight='normal', fsize=16, addOutline=True)

# add subplot letter
x, y = m(llcrnrlon+1, urcrnrlat-1)
plt.text(x, y, '(d)', size=24, ha='left', va='top', weight='normal')

# add model name
x, y = m(llcrnrlon+1, llcrnrlat+1)
plt.text(x, y, 'Leonard (2008)', size=36, ha='left', va='bottom', weight='normal')

##########################################################################################
# add 5th axis
##########################################################################################

ax = fig.add_subplot(gs1[4])
plt.tick_params(labelsize=8)

m = Basemap(projection='merc',\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='l',area_thresh=2000.)

#lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
# draw coastlines, state and country boundaries, edge of map.
m.shadedrelief()
#m.drawcoastlines()
m.drawstates()
#m.drawcountries()
m.drawparallels(arange(-90.,90.,4.), labels=[0,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,6.), labels=[0,0,0,0], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)

##########################################################################################
# plt earthquakes
##########################################################################################

idx = where((mags >= 3.0) & (year >= 1977))[0]
x, y = m(lons[idx], lats[idx])
m.plot(x, y, 'o', ms=2, mfc='0.4', mec='0.4', zorder=1)

idx = where((mags >= 5.6) & (year >= 1880))[0]
x, y = m(lons[idx], lats[idx])
m.plot(x, y, 'o', ms=7, mfc='maroon', mec='w', mew=1, zorder=2)

##########################################################################################
# add shapefiles
##########################################################################################
shpfile = 'NSHA13_Background/NSHA13_background_source_model.shp'
sf = shapefile.Reader(shpfile)
drawshapepoly(m, plt, sf,col='r',lw=1.5)

# label shapes
#labelpolygon(m, plt, sf, 'NAME', fweight='normal', fsize=16, addOutline=True)

# add subplot letter
x, y = m(llcrnrlon+1, urcrnrlat-1)
plt.text(x, y, '(e)', size=24, ha='left', va='top', weight='normal')

# add model name
x, y = m(llcrnrlon+1, llcrnrlat+1)
plt.text(x, y, 'NSHM12 Background', size=36, ha='left', va='bottom', weight='normal')

##########################################################################################
# add 6th axis
##########################################################################################

ax = fig.add_subplot(gs1[5])
plt.tick_params(labelsize=8)

m = Basemap(projection='merc',\
            llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
            urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
            rsphere=6371200.,resolution='l',area_thresh=2000.)

#lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
# draw coastlines, state and country boundaries, edge of map.
m.shadedrelief()
#m.drawcoastlines()
m.drawstates()
#m.drawcountries()
m.drawparallels(arange(-90.,90.,4.), labels=[0,0,0,0],fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)
m.drawmeridians(arange(0.,360.,6.), labels=[0,0,0,0], fontsize=16, dashes=[2, 2], color='0.5', linewidth=0.75)

##########################################################################################
# plt earthquakes
##########################################################################################

idx = where((mags >= 3.0) & (year >= 1977))[0]
x, y = m(lons[idx], lats[idx])
m.plot(x, y, 'o', ms=2, mfc='0.4', mec='0.4', zorder=1)

idx = where((mags >= 5.6) & (year >= 1880))[0]
x, y = m(lons[idx], lats[idx])
m.plot(x, y, 'o', ms=7, mfc='maroon', mec='w', mew=1, zorder=2)

##########################################################################################
# add shapefiles
##########################################################################################
shpfile = 'SinMcC2016/SIN_MCC_NSHA18_EDIT.shp'
sf = shapefile.Reader(shpfile)
drawshapepoly(m, plt, sf,col='r',lw=1.5)

# label shapes
#labelpolygon(m, plt, sf, 'NAME', fweight='normal', fsize=16, addOutline=True)

# add subplot letter
x, y = m(llcrnrlon+1, urcrnrlat-1)
plt.text(x, y, '(f)', size=24, ha='left', va='top', weight='normal')

# add model name
x, y = m(llcrnrlon+1, llcrnrlat+1)
modname = 'Sinadinovski &'+'\n'+'McCue (2016)'
plt.text(x, y, modname, size=36, ha='left', va='bottom', weight='normal')

##########################################################################################
# finish
##########################################################################################

plt.savefig('mapped_background_sources.png', format='png', bbox_inches='tight', dpi=150)
plt.show()
