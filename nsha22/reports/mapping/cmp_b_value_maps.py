from numpy import array, arange, argsort, where, delete, hstack, sqrt, \
                  unique, mean, percentile, log10, ceil, floor, \
                  nan, isnan, around, diff, interp, exp, ones_like, loadtxt, \
                  interp, floor, mgrid, ogrid, ma
from os import path, sep, mkdir, getcwd, walk
from scipy.interpolate import griddata
from shapely.geometry import Point, Polygon
#from osgeo import ogr
from datetime import datetime
from sys import argv
import shapefile
import matplotlib.pyplot as plt
from matplotlib import colors, colorbar, style
from mpl_toolkits.basemap import Basemap
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import shapefile
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes

style.use('classic')

# import non-standard functions

from mapping_tools import get_field_data, get_field_index, drawoneshapepoly, \
                          drawshapepoly, labelpolygon, get_WGS84_area
from gmt_tools import cpt2colormap 
from misc_tools import remove_last_cmap_colour, get_log_xy_locs

#from misc_tools import listdir_extension
#from make_nsha_oq_inputs import write_oq_sourcefile
        
###############################################################################
# set defaults
###############################################################################
fig = plt.figure(1, figsize=(18, 9))
gs1 = gridspec.GridSpec(1, 2)
gs1.update(wspace=0.025, hspace=0.05)

cptfile = '//Users//trev//Documents//DATA//GMT//cpt//Paired_10.cpt'
shpfile18  = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/zones/2018_mw/Domains_multi_mc/shapefiles/Domains_NSHA18_MFD.shp'
shpfile23  = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/zones/2023_mw/Domains_multi_mc/shapefiles/Domains_NSHA23_MFD.shp'

shpfiles = [shpfile18, shpfile23]

"""    
ncolours = 14
cmap, zvals = cpt2colormap(cptfile, ncolours+1)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncolours)))
"""
# get colour index
ncolours=14
b_min = 0.7
b_max = 1.4
cmap = plt.get_cmap('viridis_r', ncolours)
norm = colors.BoundaryNorm(boundaries=arange(b_min, b_max+0.1, 0.1), ncolors=ncolours)
cs = (cmap(arange(ncolours-1)))

###############################################################################
# parse shapefile and make shapely objects
###############################################################################

def parse_shp_attributes(shpfile):
    print('Reading source shapefile...')
    sf = shapefile.Reader(shpfile)
    shapes = sf.shapes()
    polygons = []
    polygonsCopy = []
    for poly in shapes:
        polygons.append(Polygon(poly.points))
        polygonsCopy.append(Polygon(poly.points))
        
    # get input arrays from shapefile
    src_code = get_field_data(sf, 'CODE', 'str')
    src_name = get_field_data(sf, 'SRC_NAME', 'str')
    src_class = get_field_data(sf, 'CLASS', 'float')
    src_rte_adj = get_field_data(sf, 'RTE_ADJ_F', 'float')
    src_usd = get_field_data(sf, 'USD', 'float')
    src_lsd = get_field_data(sf, 'LSD', 'float')
    src_overwrite_lsd = get_field_data(sf, 'OW_LSD', 'float')
    src_mmin = get_field_data(sf, 'MIN_MAG', 'float')
    src_mmin_reg = get_field_data(sf, 'MIN_RMAG', 'float')
    src_mmax = get_field_data(sf, 'MMAX_BEST', 'float')
    src_mmax_u = get_field_data(sf, 'MMAX_UPPER', 'float')
    src_mmax_l = get_field_data(sf, 'MMAX_LOWER', 'float')
    src_bval = get_field_data(sf, 'BVAL_BEST', 'float')
    src_bval_u = get_field_data(sf, 'BVAL_UPPER', 'float')
    src_bval_l = get_field_data(sf, 'BVAL_LOWER', 'float')
    src_n0 = get_field_data(sf, 'N0_BEST', 'float')
    src_n0_u = get_field_data(sf, 'N0_UPPER', 'float')
    src_n0_l = get_field_data(sf, 'N0_LOWER', 'float')
    src_bval_fix = get_field_data(sf, 'BVAL_FIX', 'float')
    src_bval_fix_sd = get_field_data(sf, 'BVAL_FIX_S', 'float') # too many chars - does not recognise "D"
    src_mcomp = get_field_data(sf, 'MCOMP', 'str')
    src_ycomp = get_field_data(sf, 'YCOMP', 'str')
    src_shmax = get_field_data(sf, 'SHMAX', 'float')
    src_shm_sig = get_field_data(sf, 'SHMAX_SIG', 'float')
    src_ymax = get_field_data(sf, 'CAT_YMAX', 'float')
    src_cat = get_field_data(sf, 'CAT_FILE', 'str')
    sortind = argsort(src_code)
    
    return src_code, src_bval, src_n0, src_class, sf

###############################################################################
# map nsha domains b-values
###############################################################################
letters = ['(a) NSHA18', '(b) NSHA23']
for i, shpfile in enumerate(shpfiles):
    # parse data
    src_code, src_bval, src_n0, src_class, sf = parse_shp_attributes(shpfile)
    
    ax = fig.add_subplot(gs1[i])
    
    # set national-scale basemap
    llcrnrlat = -45
    urcrnrlat = -7
    llcrnrlon = 105
    urcrnrlon = 153
    
    lon_0 = mean([llcrnrlon, urcrnrlon])
    lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
    lat_2 = percentile([llcrnrlat, urcrnrlat], 75)
    
    m2 = Basemap(llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
                urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,
                projection='lcc',lat_1=lat_1,lat_2=lat_2,lon_0=lon_0,
                resolution='l',area_thresh=1000.)
    
    # annotate
    m2.drawcoastlines(linewidth=0.75,color='0.25')
    m2.drawcountries(linewidth=0.75,color='0.25')
    m2.drawstates(linewidth=0.75,color='0.25')
    
    # draw parallels and meridians.
    ll_space = 6
    m2.drawparallels(arange(-90.,90.,ll_space/2.0), labels=[0,0,0,0],fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.5)
    m2.drawmeridians(arange(0.,360.,ll_space), labels=[0,0,0,0], fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.5)
    
    cindex = []
    
    # loop thru b values
    for b, n0 in zip(src_bval, src_n0):
        if not isnan(n0):
            idx = interp(b, [b_min, b_max], [0, ncolours-1])
            cindex.append(int(round(idx)))
        
    # plt source zone boundary
    drawshapepoly(m2, plt, sf, cindex=cindex, cmap=cmap, ncolours=ncolours, fillshape=True, edgecolor='r', lw=1.)
    
    xlim = ax.get_xlim()
    xtxt = xlim[1] * 0.02
    ylim = ax.get_ylim()
    ytxt = ylim[1] * 0.02
    plt.text(xtxt, ytxt, letters[i], fontsize=20, va='bottom', ha='left')
    
    # overlay domains
    sf = shapefile.Reader('/Users/trev/Documents/Geoscience_Australia/NSHA2018/source_models/zones/shapefiles/Domains/Domains_Sep2011.shp')
    drawshapepoly(m2, plt, sf, lw=1., edgecolor='r', fillshape=False)
    
##########################################################################################
# get land & lake polygons for masking
##########################################################################################
'''
# mask non-AU polygons
nonmask = [0, 5, 8, 10, 13, 14, 15, 16, 17, 20, 21] #, 2, 3, 4, 6, 7, 11, 13, 16, 17] # polygon number
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
    plt.fill(poly[:,0], poly[:,1], 'w')
    polygons.append(poly)
'''
##########################################################################################
# add domains polys
##########################################################################################

# overlay domains
sf = shapefile.Reader('/Users/trev/Documents/Geoscience_Australia/NSHA2018/source_models/zones/shapefiles/Domains/Domains_Sep2011.shp')
drawshapepoly(m2, plt, sf, lw=1., edgecolor='r', fillshape=False)

###############################################################################
# finish maps
###############################################################################

# set colorbar
plt.gcf().subplots_adjust(bottom=0.1)
cax = fig.add_axes([0.33,0.11,0.34,0.035]) # setup colorbar axes.
norm = colors.Normalize(vmin=b_min, vmax=b_max)
cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')

# set cb labels
#linticks = array([0.01, 0.03, 0.1, 0.3 ])
ticks = arange(b_min, b_max+0.1, 0.1)
cb.set_ticks(ticks)
labels = [str('%0.1f' % x) for x in ticks]

#cb.set_ticklabels(labels, fontsize=10)
cb.ax.set_xticklabels(labels, fontsize=15)
cb.set_label('b-value', fontsize=17)
  
# save figure
plt.savefig('maps/cmp_2018_2023_b_val.png', fmt='png', dpi=300, bbox_inches='tight')
   
plt.show()                
