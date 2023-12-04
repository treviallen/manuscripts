# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 11:47:02 2015

convert gridded hazard to map

Usage:
    run plt_multi_haz_maps.py <path to param file> <prob>
    
    run plt_multi_haz_maps.py multi_haz_map_PGA_regional.param 10
    

@author: tallen
"""
from sys import argv
#from matplotlib.mlab import griddata
from scipy.interpolate import griddata
from matplotlib import colors, colorbar #, cm
from os import path, mkdir, getcwd, system
#import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from numpy import arange, array, log10, mean, mgrid, ogrid, percentile, ma, isnan, nan, where, delete, floor, ceil
from tools.mapping_tools import get_map_polygons, mask_outside_polygons, cpt2colormap # drawshapepoly, labelpolygon, 
import shapefile
from scipy.constants import g
import matplotlib.gridspec as gridspec
from gmt_tools import remove_last_cmap_colour
from netCDF4 import Dataset as NetCDFFile

#from gmt_tools import cpt2colormap
#from shapely.geometry import Point, Polygon

##############################################################################
# set some default values here
##############################################################################
mpl.use('Agg')
mpl.style.use('classic')
mpl.rcParams['pdf.fonttype'] = 42

drawshape = False # decides whether to overlay seismic sources   
pltGSHAP = False # don't use this here!

# set map resolution
res = 'i' 

# get current working directory (and computer!)
cwd = getcwd()

##############################################################################
# define inputs
##############################################################################

# set with paths to hazard grids - format = plotting name; file path
#paramfile = argv[1]

# which probability - acceptable values are: 2 (2%), 9 (9.5%) or 10 (10%)
pltProbability = 10

'''
e.g.: run plt_multi_haz_maps.py multi_haz_map_PGA_regional.param 10
'''

##############################################################################
# parse param file
##############################################################################
'''
lines = open(paramfile).readlines()
hazfiles = []
modnames = []
outfile = lines[0].strip()
for line in lines[1:]:
    modnames.append(line.strip().split(';')[0])
    hazfiles.append(line.strip().split(';')[1])
'''
# set figure size based on number of models
maprows = 1
figure = plt.figure(1,figsize=(24,12*maprows+1))

gs1 = gridspec.GridSpec(maprows, 2)
#gs1.update(wspace=0.025, hspace=0.05)

hspace = -0.09 * maprows - 0.05
gs1.update(wspace=0.025, hspace=hspace) # negative looks bad in "show", but ok in pngs


pltlett = ['a)', 'b)', 'c)', 'd)', 'e)']
##############################################################################
# loop thru files and parse hazard grids
##############################################################################
hazfile = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/complete_model/2023_final/results_maps_PGA_SC_B/hazard_map-mean_1.csv'

# set axes
#ax = figure.add_subplot(maprows, 2, ii+1)
ax = figure.add_subplot(gs1[0])

# parse hazard grid file 
lines = open(hazfile).readlines()

# get keys for model
if lines[0].startswith('#'):
    line = lines[1]
else:
    line = lines[0]

# get dictionary keys
keys = line.strip().split(',')[2:]

# make grid dictionary
grddict = []
gshap = False
#print('\nReading', modnames[ii])
for line in lines[2:]:
    dat = line.strip().split(',')
    
    tmpdict = {'lon':float(dat[0]), 'lat':float(dat[1])}
    
    # fill keys
    idx = 2
    for key in keys:
        if gshap == True:
            # convert to m/s**2
            tmpdict[key] = float(dat[idx]) / g
        else:
            tmpdict[key] = float(dat[idx])
        idx += 1
    
    # add to grid list
    grddict.append(tmpdict)
    
##############################################################################    
# get key index for plotting
##############################################################################

for i, key in enumerate(keys):
    keyProb = str(int(floor(100*float(key.split('-')[-1]))))
    if keyProb == str(pltProbability):
        mapidx = i

##############################################################################    
# now make maps
##############################################################################

#keys = ['PGA_10', 'PGA_02', 'SA02_10', 'SA02_02', 'SA10_10', 'SA10_02']
for i, key in enumerate([keys[mapidx]]): # just plot 1 for now!
    if i > 0:
        plt.clf()
        plt.cla()

    # get IM period
    period = key.split('-')[0]
    period = period.replace('(','')
    period = period.replace(')','')
    period = period.replace('.','')
    
    # get map probability of exceedance
    probFraction = str(float(key.split('-')[-1]))
    probability = str(100*float(key.split('-')[-1])).split('.')[0]+'%'
    #probability = str(100*float(key.split('-')[-1]))+'%'
    if probability == '9%':
        probability = '9.5%'
    print('Probability', probability)
    
    bbox = '108/152/-44/-8' # map boundary - lon1/lon2/lat1/lat2
    bbox = '107.0/153.0/-43.0/-8.0'
    
    bbox = bbox.split('/')
    minlon = float(bbox[0])
    maxlon = float(bbox[1])
    minlat = float(bbox[2])
    maxlat = float(bbox[3])
    mbuff = 1.
    
    # build data to plot
    hazvals = []
    latlist = []
    lonlist = []
    
    # add buffer to data
    for gridval in grddict:
        lonlist.append(gridval['lon'])
        latlist.append(gridval['lat'])
        if gridval[key] == 0.0:
            hazvals.append(0.0)
        else:
            hazvals.append(gridval[key])
            
    #idx = array(range(0, len(lonlist), 10)) # resample for quickly testing mapping
    idx = array(range(0, len(lonlist), 1))
    lonlist = array(lonlist)[idx]
    latlist = array(latlist)[idx]
    hazvals = array(hazvals)[idx]
    
    # delete zero hazvals
    idx =where(hazvals==0)[0]
    '''
    lonlist = delete(lonlist, idx)
    latlist = delete(latlist, idx)
    hazvals = delete(hazvals, idx)
    '''
    hazvals[idx] = 1E-20
    
    # get map bounds
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
    
    # draw parallels and meridians.
    if maxlon-minlon > 40:
        xlabel = 6.
    elif maxlon-minlon > 20:
        xlabel = 4.
    elif maxlon-minlon > 10:
        xlabel = 2.
    else:
        xlabel = 1.
        
    if maxlat-minlat > 40:
        ylabel = 6.
    elif maxlat-minlat > 20:
        ylabel = 4.
    elif maxlat-minlat > 10:
        ylabel = 2.
    else:
        ylabel = 1.
            
    m.drawparallels(arange(-90.,90.,ylabel), labels=[0,0,0,0],fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.5)
    m.drawmeridians(arange(0.,360.,xlabel), labels=[0,0,0,0], fontsize=10, dashes=[2, 2], color='0.5', linewidth=0.5)
    
    # first make regular cartesian grid
    print('Resampling data...')
    N = 500j
    extent = (minlon-mbuff, maxlon+mbuff, minlat-mbuff, maxlat+mbuff)
    xs,ys = mgrid[extent[0]:extent[1]:N, extent[2]:extent[3]:N]
    	
    #resampled = griddata(lonlist, latlist, hazvals, xs, ys, interp='linear')
    resampled = griddata((lonlist, latlist), hazvals, (xs, ys), method='linear')
    #resampled = griddata(lonlist, latlist, log10(hazvals), lonlist, latlist, interp='linear') # if this suddenly works, I have no idea why!
    
    # get 1D lats and lons for map transform
    lons = ogrid[extent[0]:extent[1]:N]
    lats = ogrid[extent[2]:extent[3]:N]
    
    # transform to map projection
    nx = int((m.xmax-m.xmin)/2000.)+1
    ny = int((m.ymax-m.ymin)/2000.)+1
    
    # differences in the way different machines deal with grids - weird!
    if cwd.startswith('/nas'):
        transhaz = m.transform_scalar(resampled.T,lons,lats,nx,ny)
    else:
        transhaz = m.transform_scalar(resampled.T,lons,lats,nx,ny)
    
    masked_array = ma.array(transhaz, mask=isnan(transhaz))
    #masked_array = masked_array.set_fill_value(0)
    
    # get colormap from cpt file
    cptfile = '/Users/trev/Documents/Geoscience_Australia/NSHA2018/postprocessing/maps/cw1-013_mod.cpt'
    #ncols = 9
    
    #cmap = cm.rainbow
    #print(period
    if period == 'PGA':
    
        if probability == '10%' or probability == '9.5%': # kluge to get on same scale
            ncolours = 12
            #vmin = -2.
            #vmax = -0.25 - 0.125 # so there is an odd number for which to split the cpt
            
        elif probability == '2%':
            ncolours = 12
            #vmin = -1.5
            #vmax = vmin + 0.25 * ncolours/2.
        T = 'PGA'
    
    elif period == 'SA005':
        
        if probability == '10%' or probability == '9.5%': # kluge to get on same scale
            ncolours = 12
            
        elif probability == '2%':
            ncolours = 10
        T = 'Sa(0.05)'
        
    elif period == 'SA01':
        
        if probability == '10%' or probability == '9.5%': # kluge to get on same scale
            ncolours = 12
            
        elif probability == '2%':
            ncolours = 10
        T = 'Sa(0.1)'
        
    elif period == 'SA02':
        if probability == '10%' or probability == '9.5%':
            ncolours = 13
            
        elif probability == '2%':
            ncolours = 12
            
        T = 'Sa(0.2)'
        
    elif period == 'SA10':
        if probability == '10%' or probability == '9.5%':
            ncolours = 13
        
        elif probability == '2%':
            ncolours = 14
            
        T = 'Sa(1.0)'
    
    #ncolours = 13
        
    try:
        cmap, zvals = cpt2colormap(cptfile, ncolours, rev=True)
    except:
        try:
            if pltGSHAP == 'True':
                nascptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/gshap_mpl.cpt'
                ncolours = 10
                cmap, zvals = cpt2colormap(nascptfile, ncolours, rev=False)
                cmap = remove_last_cmap_colour(cmap)
                
            else:
                nascptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/postprocessing/maps/'+ cptfile
                cmap, zvals = cpt2colormap(nascptfile, ncolours, rev=True)
            
            #cptfile = '/nas/gemd/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/postprocessing/maps/GMT_no_green.cpt'
            
        except:
            try:
                ncicptfile = '/short/w84/NSHA18/sandpit/tia547/NSHA2018/postprocessing/maps/'+ cptfile
                cmap, zvals = cpt2colormap(ncicptfile, ncolours, rev=True)
   
            except:
                ncicptfile = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/postprocessing/maps/'+ cptfile
                cmap, zvals = cpt2colormap(ncicptfile, ncolours, rev=True)
    
    print('Making map...')
    cmap.set_bad('w', 1.0)
    
    if probability == '10%' or probability == '9.5%':
        if pltGSHAP == 'True':
            bounds = array([0., 0.2, 0.4, 0.8, 1.6, 2.4, 3.2, 4.0, 4.8, 6.0])
            #ncolours = 9
            #norm = colors.Normalize(vmin=0,vmax=10)
        else:
            bounds = array([0, 0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.12, 0.16, 0.24])
            ncolours = 12
            norm = colors.BoundaryNorm(boundaries=bounds, ncolors=ncolours)
    else:
        bounds = array([0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.16, 0.20, 0.3, 0.5, 0.7, 1.0])
        #ncolours = 12
    norm = colors.BoundaryNorm(boundaries=bounds, ncolors=ncolours)
    
    #m.imshow(masked_array, cmap=cmap, extent=extent, vmin=vmin, vmax=vmax, zorder=0)
    m.imshow(masked_array, cmap=cmap, extent=extent, norm=norm, zorder=0)
    
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
    '''
    titlestr = ' '.join((modnames[i], T, probability, 'in 50-Year Mean Hazard on AS1170.4 Site Class '))    
    plt.title(titlestr+'$\mathregular{B_e}$')
    '''
    
    # plt letter
    xlim = ax.get_xlim()
    xtxt = xlim[1] * 0.02
    ylim = ax.get_ylim()
    ytxt = ylim[1] * 0.96
    plt.text(xtxt, ytxt, pltlett[0], fontsize=30, va='top', ha='left')
    '''
    xtxt = xlim[1] * 0.02
    ylim = ax.get_ylim()
    ytxt = ylim[1] * 0.02
    plt.text(xtxt, ytxt, modnames[ii], fontsize=27, va='bottom', ha='left')
    '''
    ##########################################################################################
    # plt subduction profile - comment out when not using
    ##########################################################################################
    '''
    proffile = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/shared/north_aus_profile.csv'
    proflons = []
    proflats = []
    lines = open(proffile).readlines()
    for line in lines:
        dat = line.strip().split(',')
        
        proflons.append(float(dat[0]))
        proflats.append(float(dat[1]))
        
    proflons = array(proflons)
    proflats = array(proflats)
    
    # plot profile on map
    x,y = m(proflons, proflats)
    m.plot(x,y,'k+', ms=7)
    '''
    # get map bbox
    if i == 0:
        map_bbox = ax.get_position().extents
    
'''
###########################################################################################
make colourbar
###########################################################################################
'''    

# set colourbar
plt.gcf().subplots_adjust(bottom=0.05)
cby = 0.10 + 0.02/maprows + 0.017*(maprows-1)
cbh = 0.00 + 0.027 / maprows
cax = figure.add_axes([0.14,cby,0.35,cbh]) # setup colorbar axes.
#norm = colors.Normalize(vmin=vmin, vmax=vmax)
cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')

cb.set_ticks(bounds)
labels = ['0'+str('%0.3f' % x).strip('0') for x in bounds]
labels[0] = '0.0'
if bounds[-1] >= 1.0:
    labels[-1] = '1.0'
cb.set_ticklabels(labels)
cb.ax.tick_params(labelsize=18)

# set title
titlestr = ' '.join((T, probability, 'in 50-Year Mean Hazard (g)'))
cb.set_label(titlestr, fontsize=22)
       
'''
###########################################################################################
plot ratio map
###########################################################################################
'''
ax = figure.add_subplot(gs1[1])
def parse_oq_hazard_grid(hazfile, pltprob):
    
    # parse csv files
    lines = open(hazfile).readlines()
    
    # get keys for model
    if lines[0].startswith('#'):
        line = lines[1]
    else:
        line = lines[0]
    
    # get dictionary keys
    keys = line.strip().split(',')#[2:]
    print(keys)
    	
    '''
    for i, key in enumerate(keys):
        keyProb = str(int(floor(100*float(key.split('-')[-1]))))
        if keyProb == pltProbability:
            mapidx = i
    '''    
    # make grid dictionary
    grddict = []
    
    #print('\nReading', modnames[ii])
    for line in lines[2:]:
        dat = line.strip().split(',')
        #print(dat)
        tmpdict = {'lon':float(dat[0]), 'lat':float(dat[1])}
        
        # fill keys
        #idx = 2
        for i, key in enumerate(keys):
            #print(key)
            #if i == mapidx:
            if i >= 2:
                tmpdict[key] = float(dat[i])
        
        # add to grid list
        grddict.append(tmpdict)
        
    return grddict
    
def dict2netcdf(gridDict, outGrid, zKey):
    from os import system
    
    # write temp file and make netcdf
    csvtxt = ''
    for gd in gridDict:
        csvtxt += ','.join((str(gd['lon']), str(gd['lat']), str(gd[zKey]))) + '\n'
    f = open('grid.csv', 'w')
    f.write(csvtxt)
    f.close()
    
    # write netcdf
    system('gmt5 surface grid.csv -G'+outGrid+' -R110/156/-46/-9 -I0.05')
    system('gmt5 grdmath '+outGrid+' 0.002 MAX = '+outGrid)

##############################################################################    
# make grids
##############################################################################

pltProbability = '10'
if pltProbability == '10':
    grid1 = 'nsha23_0.1_pga.grd'
    grid2 = 'nsha18_0.1_pga.grd'
    pltkey = 'PGA-0.1'
else:
    grid1 = 'nsha23_0.02_pga.grd'
    grid2 = 'nsha18_0.02_pga.grd'
    pltkey = 'PGA-0.02'

# read file 1
'''
#hazfile1 = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/complete_model/2023_final/results_maps_PGA_ta_pref/hazard_map-mean_1.csv'
hazfile1 = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/complete_model/2023_final/results_maps_PGA_SC_B/hazard_map-mean_1.csv'
nsha23_grddict = parse_oq_hazard_grid(hazfile1, pltProbability)
dict2netcdf(nsha23_grddict, grid1, pltkey)

# read file 2
hazfile2 = '/Users/trev/Documents/Geoscience_Australia/NSHA2018/source_models/complete_model/final/results_maps_PGA/hazard_map-mean_PGA.csv'
nsha18_grddict = parse_oq_hazard_grid(hazfile2, pltProbability)
dict2netcdf(nsha18_grddict, grid2, pltkey)
'''
system(' '.join(('gmt5 grdmath', grid1, grid2, ' DIV ABS au_land_mask.grd MUL = map_ratio.grd'))) #au_land_mask.grd MUL 

# set small values to zero
system('gmt5 grdmath map_ratio.grd 0.01 GE map_ratio.grd MUL = map_ratio.grd')

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

            
##############################################################################    
# now make maps
##############################################################################

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
                
m.drawparallels(arange(-90.,90.,6), labels=[0,0,0,0],fontsize=14, dashes=[2, 2], color='0.5', linewidth=0.5)
m.drawmeridians(arange(0.,360.,6), labels=[0,0,0,0], fontsize=14, dashes=[2, 2], color='0.5', linewidth=0.5)

##############################################################################    
# read netcdf
##############################################################################
bounds = array([1/500.0, 1/2.75, 1/2.5, 1/2.25, 1/2.0, 1/1.75, 1/1.5, 1/1.25, \
            1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 500.0]) #nsha = 1/50
            
norm = colors.BoundaryNorm(boundaries=bounds, ncolors=ncolours)

print( 'Reading netCDF file...')
nc = NetCDFFile('map_ratio.grd')

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

plt.text(xtxt, ytxt, pltlett[1], fontsize=30, va='top', ha='left')

'''
###########################################################################################
make colourbar
###########################################################################################
'''    

plt.gcf().subplots_adjust(bottom=0.05)
cby = 0.10 + 0.02/maprows + 0.017*(maprows-1)
cbh = 0.00 + 0.027 / maprows
cax = figure.add_axes([0.535,cby,0.35,cbh]) # setup colorbar axes.
#norm = colors.Normalize(vmin=vmin, vmax=vmax)
cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')

cb.set_ticks(bounds[::2])
labels = [str('%0.2f' % x) for x in bounds[::2]]
labels[0] = '0.33'
labels[-1] = '3.00'
cb.set_ticklabels(labels)
cb.ax.tick_params(labelsize=18)

# set title
titlestr = 'Hazard Ratio (NSHA23 / NSHA18)'
cb.set_label(titlestr, fontsize=22)
 
# now save png file
plt.savefig(path.join('map_nsha23_w_ratio.png'), \
            dpi=200, format='png', bbox_inches='tight')

# save pdf file
'''
plt.savefig(path.join('maps', 'hazard_map_'+modelName.replace(' ','_')+'.'+key+'.pdf'), \
            dpi=300, format='pdf', bbox_inches='tight')
'''
plt.show()
plt.close('all')

#