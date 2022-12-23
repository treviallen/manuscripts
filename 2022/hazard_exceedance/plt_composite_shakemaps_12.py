# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 11:47:02 2015

convert gridded hazard to map

Usage:
    python map_nsha18.py <path to csv file>
    

@author: tallen
"""
from sys import argv
from matplotlib import colors, colorbar #, cm
from os import path, mkdir, getcwd
#import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset as NetCDFFile
from numpy import arange, array, log10, mean, mgrid, ogrid, percentile, ma, isnan, nan, where, delete, floor, ceil
from mapping_tools import get_map_polygons, mask_outside_polygons, cpt2colormap, labelCentroid # drawshapepoly, labelpolygon, 
import shapefile
from scipy.constants import g
import matplotlib.gridspec as gridspec

#from gmt_tools import cpt2colormap
#from shapely.geometry import Point, Polygon

##############################################################################
# set some default values here
##############################################################################
#mpl.use('Agg')
mpl.style.use('classic')
mpl.rcParams['pdf.fonttype'] = 42

drawshape = False # decides whether to overlay seismic sources   
pltGSHAP = False # don't use this here!

# set map resolution
res = 'f' 
#res = 'l'

# get current working directory (and computer!)
cwd = getcwd()

# set figure size based on number of models
maprows = 1
mapcols = 2
figure = plt.figure(1,figsize=(15,8))

gs1 = gridspec.GridSpec(maprows, mapcols)
#gs1.update(wspace=0.025, hspace=0.05)

hspace = -0.09 * maprows - 0.05
gs1.update(wspace=0.025, hspace=hspace) # negative looks bad in "show", but ok in pngs

#ncfiles = ['max_pga_grid_au.grd', 'max_pga_grid_ceus.grd']
ncfiles = ['max_pga_grid_au_banda.grd', 'max_pga_grid_ceus_banda.grd']
           
pltlett = ['a)', 'b)', 'c)', 'd)', 'e)']
modnames = ['AU GMMs', 'AU & CEUS GMMs']

##############################################################################
# loop thru files and parse hazard grids
##############################################################################

for ii, ncf in enumerate(ncfiles):

    # set axes
    ax = figure.add_subplot(gs1[ii])
                
    ##############################################################################    
    # now make maps
    ##############################################################################
    
    # get map bounds
    bbox = '107.0/153.0/-43.0/-8.0'
        
    bbox = bbox.split('/')
    minlon = float(bbox[0])
    maxlon = float(bbox[1])
    minlat = float(bbox[2])
    maxlat = float(bbox[3])
        
    llcrnrlat = minlat
    urcrnrlat = maxlat
    llcrnrlon = minlon
    urcrnrlon = maxlon
    lon_0 = mean([llcrnrlon, urcrnrlon])
    lon_0 = 134.
    lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
    lat_2 = percentile([llcrnrlat, urcrnrlat], 75)
    
    # set map
    # Projection used for National Mapping
    m = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
                urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,
                projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
                resolution=res,area_thresh=2000.)
                
    #m.drawmapboundary(fill_color='lightgray')
    #m.fillcontinents(color='white',lake_color='lightgray',zorder=0)
    m.drawcoastlines(linewidth=0.5,color='k')
    m.drawcountries(color='0.2')
    m.drawstates(color='0.2')
                    
    m.drawparallels(arange(-90.,90.,6), labels=[0,0,0,0],fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.5)
    m.drawmeridians(arange(0.,360.,6), labels=[0,0,0,0], fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.5)
    
    # get colormap from cpt file
    cptfile = '/Users/trev/Documents/Geoscience_Australia/NSHA2018/postprocessing/maps/cw1-013_mod.cpt'
    ncolours = 12
    cmap, zvals = cpt2colormap(cptfile, ncolours, rev=True)
    
    print('Making map...')
    cmap.set_under('w', 1.0)
    
    bounds = array([1E-3, 0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.12, 0.16, 0.24])
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=ncolours)
    
    ##############################################################################    
    # read netcdf
    ##############################################################################
    
    print( 'Reading netCDF file...')
    nc = NetCDFFile(ncf)
    
    data = nc.variables['z'][:]
    lons = nc.variables['x'][:]
    lats = nc.variables['y'][:]
    
    # transform to metres
    nx = int((m.xmax-m.xmin)/500.)+1
    ny = int((m.ymax-m.ymin)/500.)+1
    
    hazdat = m.transform_scalar(data,lons,lats,nx,ny)
    
    im = m.imshow(hazdat, cmap=cmap, norm=norm, alpha=1.0)
    
    ##########################################################################################
    # get land & lake polygons for masking
    ##########################################################################################
    # mask non-AU polygons
    nonmask = [0, 1, 2, 3, 4, 6, 7, 11, 13, 16, 17] # polygon number
    nonmask = [0, 9, 15, 18, 19] # polygon number
    landpolys = []
    for pidx, polygon in enumerate(m.landpolygons):
        maskPoly = True
        
        # test poly number
        #poly = polygon.get_coords()
        #plt.text(poly[0,0], poly[0,1], str(pidx), fontsize=15)
        
        for nmidx in nonmask:
            if pidx == nmidx:
                maskPoly = False 
        if maskPoly == True:
            poly = polygon.get_coords()
            plt.fill(poly[:,0], poly[:,1], 'w')
    
    #mask_outside_polygon(polys[1][::-1], ax=None)
    polys = get_map_polygons(m)
    mask_outside_polygons(polys, '0.9', plt)
    
    # get lake ploygons
    polygons = []
    for polygon in m.lakepolygons:
        poly = polygon.get_coords()
        plt.fill(poly[:,0], poly[:,1], '0.9')
        polygons.append(poly)
    
    ##########################################################################################
    # format main axis
    ##########################################################################################
    
    # plt letter
    xlim = ax.get_xlim()
    xtxt = xlim[1] * 0.02
    ylim = ax.get_ylim()
    ytxt = ylim[1] * 0.96
    plt.text(xtxt, ytxt, pltlett[ii], fontsize=30, va='top', ha='left')
    
    xtxt = xlim[1] * 0.02
    ylim = ax.get_ylim()
    ytxt = ylim[1] * 0.02
    plt.text(xtxt, ytxt, modnames[ii], fontsize=23, va='bottom', ha='left')
    
'''
###########################################################################################
make colourbar
###########################################################################################
'''    

# set colourbar
plt.gcf().subplots_adjust(bottom=0.1)
cax = figure.add_axes([0.2,0.1,0.6,0.05]) # setup colorbar axes.
#norm = colors.Normalize(vmin=vmin, vmax=vmax)
cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')

cb.set_ticks(bounds)
labels = ['0'+str('%0.3f' % x).strip('0') for x in bounds]
labels[0] = '0.0'
cb.set_ticklabels(labels)
cb.ax.tick_params(labelsize=17)

# set title
titlestr = '1972-2021 Observed PGA (g)'
cb.set_label(titlestr, fontsize=24)
'''
# check to see if maps exists
if path.isdir('maps') == False:
    mkdir('maps')
'''    
# now save png file
plt.savefig(path.join('figures/composite_shakemaps.png'), dpi=600, format='png', bbox_inches='tight')

plt.show()
#plt.close('all')

#