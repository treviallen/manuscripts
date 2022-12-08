# -*- coding: utf-8 -*-
from sys import argv
from matplotlib import colors, colorbar #, cm
from os import path, mkdir, getcwd, system
#import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset as NetCDFFile
from numpy import arange, array, log10, mean, mgrid, ogrid, percentile, ma, isnan, nan, where, delete, floor, ceil
from mapping_tools import get_map_polygons, mask_outside_polygons, cpt2colormap # drawshapepoly, labelpolygon, 
from misc_tools import remove_last_cmap_colour
import shapefile
from scipy.constants import g
import matplotlib.gridspec as gridspec

##############################################################################
# set some default values here
##############################################################################
#mpl.use('Agg')
mpl.style.use('classic')
mpl.rcParams['pdf.fonttype'] = 42

# set map resolution
res = 'i' 

# get current working directory (and computer!)
cwd = getcwd()

# set figure size based on number of models
maprows = 2
mapcols = 2
figure = plt.figure(1,figsize=(15,15))

gs1 = gridspec.GridSpec(maprows, mapcols)
#gs1.update(wspace=0.025, hspace=0.05)

hspace = -0.09 * maprows - 0.05
gs1.update(wspace=0.025, hspace=hspace) # negative looks bad in "show", but ok in pngs

ncfiles = ['haz_maps/gaull90/gaull90_interp.0.05.grd',
           'haz_maps/gshap/gshap_interp.0.05.grd',
           'haz_maps/nshm12/nshm12_resampled.0.05.grd',
           'haz_maps/nsha18/nsha18_interp.0.05.grd']
#ncfiles = ['haz_maps/gshap/gshap_interp.0.05.grd'] # for testing
           
pltlett = ['a)', 'b)', 'c)', 'd)', 'e)']
modnames = ['Gaull et al (1990)', 'GSHAP (1999)', 'NSHM12', 'NSHA18']

##############################################################################
# get cmap
##############################################################################

# get colormap from cpt file
cptfile = '/Users/trev/Documents/DATA/GMT/cpt/BlueWhiteOrangeRed.cpt'
cptfile = '/Users/trev/Documents/DATA/GMT/cpt/blue-tan-d15.cpt'
#cptfile = '/Users/trev/Documents/DATA/GMT/es_landscape_90.cpt'
cptfile = '/Users/trev/Documents/DATA/GMT/cpt/cool-warm-d15.cpt'

ncolours = 17
cmap, zvals = cpt2colormap(cptfile, ncolours, rev=False)
cmap = remove_last_cmap_colour(cmap)


print('Making map...')

#logbounds = log(bounds)

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
    
    ##############################################################################    
    # read netcdf
    ##############################################################################
    
    print( 'Generating ratio grid...')
    # calculate log ratios on the fly
    system('rm sm_div_haz'+str(ii)+'.grd')

    system('gmt5 grdmath max_pga_grid_au_banda.grd ' + ncf + ' DIV ABS au_land_mask.grd MUL = sm_div_haz'+str(ii)+'.grd') #au_land_mask.grd MUL 
    
    # ste small values to zero
    system('gmt5 grdmath sm_div_haz'+str(ii)+'.grd 0.01 GE sm_div_haz'+str(ii)+'.grd MUL = sm_div_haz'+str(ii)+'.grd')
    
    #system('gmt5 grdsample sm_div_haz.grd -Gsm_div_haz.grd -R110/156/-46/-9 -I0.05')
    #system('gmt5 grd2xyz sm_div_haz.grd | gmt5 surface -Rsm_div_haz.grd -I0.05 -Gsm_div_haz.grd')
    
    print( 'Reading netCDF file...')
    nc = NetCDFFile('sm_div_haz'+str(ii)+'.grd')
    
    try:
       data = nc.variables['z'][:]
       lons = nc.variables['x'][:]
       lats = nc.variables['y'][:]
    except:
       data = nc.variables['z'][:]
       lons = nc.variables['lon'][:]
       lats = nc.variables['lat'][:]
    
    # transform to metres
    nx = int((m.xmax-m.xmin)/500.)+1
    ny = int((m.ymax-m.ymin)/500.)+1
    
    hazdat = m.transform_scalar(data,lons,lats,nx,ny)
    
    # set cmap
    bounds = array([1/500.0, 1/2.75, 1/2.5, 1/2.25, 1/2.0, 1/1.75, 1/1.5, 1/1.25, \
                1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 500.0]) #nsha = 1/50

    #cmap = plt.get_cmap('bwr', len(bounds)-1) # original submission used this!!!!!!!!!!!!!!!!!
    #cmap = plt.get_cmap('coolwarm', len(bounds)-1)

    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=ncolours)
    cmap.set_under('w', 1.0)
    #cmap.set_over('w', 1.0)
    cmap.set_bad('w', 1.0)
    
    im = m.imshow(hazdat, cmap=cmap, norm=norm, alpha=1.0)
    
    ##########################################################################################
    # get land & lake polygons for masking
    ##########################################################################################
    # mask non-AU polygons
    nonmask = [0, 1, 2, 3, 4, 6, 7, 11, 13, 16, 17] # polygon number
    nonmask = [0, 9, 14, 17, 18, 19] # polygon number
    landpolys = []
    for pidx, polygon in enumerate(m.landpolygons):
        maskPoly = True
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
    plt.text(xtxt, ytxt, modnames[ii], fontsize=27, va='bottom', ha='left')
    
'''
###########################################################################################
make colourbar
###########################################################################################
'''    

# set colourbar
plt.gcf().subplots_adjust(bottom=0.1)
cax = figure.add_axes([0.2,0.11,0.6,0.03]) # setup colorbar axes.
#norm = colors.Normalize(vmin=vmin, vmax=vmax)
cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal') #, extend='both')

cb.set_ticks(bounds[::2])
labels = [str('%0.2f' % x) for x in bounds[::2]]
labels[0] = '0.33'
labels[-1] = '3.00'
cb.set_ticklabels(labels)
cb.ax.tick_params(labelsize=17)

# set title
titlestr = 'Hazard Ratio (Observed / Predicted)'
cb.set_label(titlestr, fontsize=24)
'''
# check to see if maps exists
if path.isdir('maps') == False:
    mkdir('maps')
'''    
# now save png file
plt.savefig('figures/multi_au_ratio_maps.png', dpi=300, format='png', bbox_inches='tight')

plt.show()
#plt.close('all')

#