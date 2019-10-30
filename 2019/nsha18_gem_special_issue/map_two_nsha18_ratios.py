# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 11:47:02 2015

convert gridded hazard to map

Usage:
    python map_nsha18.py <path to csv file>
    

@author: tallen
"""
from sys import argv
#from matplotlib.mlab import griddata
from scipy.interpolate import griddata
from matplotlib import colors, colorbar #, cm
from os import path, mkdir, getcwd
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.basemap import Basemap
from numpy import arange, array, log10, mean, mgrid, ogrid, percentile, ma, isnan, logspace, nan, linspace, ceil
from mapping_tools import get_map_polygons, mask_outside_polygons, cpt2colormap # drawshapepoly, labelpolygon, 
from misc_tools import remove_last_cmap_colour
import matplotlib.gridspec as gridspec
#from gmt_tools import cpt2colormap
#from shapely.geometry import Point, Polygon

##############################################################################
# set some default values here
##############################################################################
mpl.use('Agg')
mpl.style.use('classic')
mpl.rcParams['pdf.fonttype'] = 42
gs1 = gridspec.GridSpec(1, 2)
hspace = 0.15
gs1.update(wspace=0.03, hspace=hspace) # negative looks bad in "show", but ok in pngs


drawshape = False # decides whether to overlay seismic sources

bbox = '107.0/153.0/-45.0/-7.0' # map boundary - lon1/lon2/lat1/lat2
bbox = bbox.split('/')
minlon = float(bbox[0])
maxlon = float(bbox[1])
minlat = float(bbox[2])
maxlat = float(bbox[3])
mbuff = 2.

pltlett = ['a)', 'b)', 'c)', 'd)', 'e)']

# set map resolution
res = 'i' 

# get current working directory (and computer!)
cwd = getcwd()

addContours = False # for testing only

##############################################################################
# parse hazard map file
##############################################################################

# set map file to plot
gridfiles2 = ['/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/source_models/complete_model/final/results_maps_PGA_cat_test/hazard_map-mean_1.csv', \
             '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/source_models/complete_model/final/results_maps_PGA/hazard_map-mean_1.csv']

# set map file to plot
gridfiles1 = ['/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/source_models/complete_model/final_mx/results_maps_PGA/hazard_map-mean_1.csv', \
              '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/source_models/complete_model/final/results_maps_PGA_gaull/hazard_map-mean_1.csv']
              #]

# get map name for plotting
modelName = 'NSHA18_ratios'

# set lims
clim = 4. # valid inputs 3 or 4

##############################################################################
# parse hazard map files and get ratio
##############################################################################

# get model name from input file
#model = path.split(gridfile1)[-1].split('_')[2].split('.')[0] # this will likely need modifying depending on filename format

def parse_grd_file(gridfile):
    # parse hazard grid file 
    lines = open(gridfile).readlines()
    
    # get keys for model
    if lines[0].startswith('#'):
        line = lines[1]
    else:
        line = lines[0]
    
    # get dictionary keys
    keys = line.strip().split(',')[2:]
    
    # make grid dictionary
    grddict = []
    
    print('\nReading', modelName)
    for line in lines[2:]:
        dat = line.strip().split(',')
        tmpdict = {'lon':float(dat[0]), 'lat':float(dat[1])}
        
        # fill keys
        idx = 2
        for key in keys:
            tmpdict[key] = float(dat[idx])
            idx += 1
        
        # add to grid list
        grddict.append(tmpdict)
        
    return grddict, keys
    
##############################################################################    
# now begin loops
##############################################################################
fig = plt.figure(1, figsize=(20, 12))

j = 0
for gridfile1, gridfile2 in zip(gridfiles1, gridfiles2):
    
    grddict1, keys  = parse_grd_file(gridfile1)
    grddict2, keys  = parse_grd_file(gridfile2)
    
    # make ratio grddict - just 10% in 50 for now!!!
    print('\nJust 10% in 50 for now!!!\n')
    grddict = []
    for gd1, gd2 in zip(grddict1, grddict2):
        if gd1['lon'] == gd2['lon'] and gd1['lat'] == gd2['lat']:
            if gd2[keys[0]] == 0.0:
                tdict = {'lon': gd1['lon'], 'lat': gd1['lat'],
                         'logratio': nan}
            else:
                 tdict = {'lon': gd1['lon'], 'lat': gd1['lat'],
                         'logratio': log10(gd1[keys[0]] / gd2[keys[0]])}
                         
            grddict.append(tdict)
                   
    ##############################################################################    
    # now make maps
    ##############################################################################
    
    #keys = ['PGA_10', 'PGA_02', 'SA02_10', 'SA02_02', 'SA10_10', 'SA10_02']
    for i, key in enumerate([keys[0]]): # just plot 1 for now!
        
        # get IM period
        period = key.split('-')[0]
        T = 'PGA'
        
        # get map probability of exceedance
        probability = str(100*float(key.split('-')[-1])).split('.')[0]+'%'
        
        #figure = plt.figure(i,figsize=(19,12))
        
        ax = fig.add_subplot(gs1[j])
        
        '''
        # get shpfile for masking hazard values
        shpfile = solfile.split('.')[0]
        
        inshape = '../Grids/2005_grids/canada_2005grid_released.shp'
        sf = shapefile.Reader(inshape)
        sf = sf.shapes()
        maskpoly = Polygon(sf[0].points)
        '''
        
        # build data to plot
        logRatio = []
        latlist = []
        lonlist = []
        
        # add buffer to data
        for gridval in grddict:
            lonlist.append(gridval['lon'])
            latlist.append(gridval['lat'])
            if isnan(gridval['logratio']):
                logRatio.append(nan)
            else:
                logRatio.append(gridval['logratio'])
                
            
        #idx = array(range(0, len(lonlist), 100)) # resample for quickly testing mapping
        idx = array(range(0, len(lonlist), 1))
        lonlist = array(lonlist)[idx]
        latlist = array(latlist)[idx]
        logRatio = array(logRatio)[idx]
        
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
        resampled = griddata(lonlist, latlist, logRatio, xs, ys, interp='linear')
        
        
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
        #Optimus Prime
        cptfile = 'BlueWhiteOrangeRed.cpt'
        #cptfile = 'BlueYellowRed.cpt'
        
        if clim == 4.:
            ncolours = 17
            nticks = 9
        elif clim == 3.:
            ncolours = 17
            nticks = 9
        
        # get colour lims in log10
        vmin = log10(1./clim)
        vmax = log10(clim)
            
        try:
            cmap, zvals = cpt2colormap(cptfile, ncolours, rev=False)
        except:
            try:
                nascptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/postprocessing/maps/'+ cptfile
                capfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/shared/capitals_names.csv'
                cmap, zvals = cpt2colormap(nascptfile, ncolours, rev=False)
                #cmap, zvals = cpt2colormap(nascptfile, ncolours, rev=True)
                cmap = remove_last_cmap_colour(cmap)
            except:
                ncicptfile = '/short/w84/NSHA18/sandpit/tia547/NSHA2018/postprocessing/maps/'+ cptfile
                cmap, zvals = cpt2colormap(ncicptfile, ncolours, rev=False)
                cmap = remove_last_cmap_colour(cmap)
    
        
        print('Making map...')
        cmap.set_bad('w', 1.0)
        m.imshow(masked_array, cmap=cmap, extent=extent, vmin=vmin, vmax=vmax, zorder=0)
        
        ##########################################################################################
        # get land & lake polygons for masking
        ##########################################################################################
        # mask non-AU polygons
        nonmask = [0, 2, 3, 4, 6, 7, 11, 13, 16, 17] # polygon number
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
        
        '''
        ###########################################################################################
        annotate cities
        ###########################################################################################
        '''
        
        import matplotlib.patheffects as PathEffects
        pe = [PathEffects.withStroke(linewidth=2.5, foreground="w")]
                  
        llat = []
        llon = []
        locs = []
        textoffset = []
        
        # read data
        lines = open(capfile).readlines()
        for line in lines:
            llon.append(float(line.strip().split(',')[0]))
            llat.append(float(line.strip().split(',')[1]))
            locs.append(line.strip().split(',')[2])
            textoffset.append(float(line.strip().split(',')[3]))
        
        # plot locs on map
        x, y = m(array(llon), array(llat))
        plt.plot(x, y, 's', markerfacecolor='None', markeredgecolor='k', markeredgewidth=1., markersize=8)
        
        # label cities
        for i, loc in enumerate(locs):
            if textoffset[i] == 0.:
                x, y = m(llon[i]-0.35, llat[i]+0.12)
                plt.text(x, y, loc, size=15, ha='right', va='bottom', weight='normal')
            elif textoffset[i] == 1.:
                x, y = m(llon[i]+0.35, llat[i]+0.12)
                #plt.text(x, y, loc, size=15, ha='left', va='bottom', weight='light', path_effects=pe)
                plt.text(x, y, loc, size=15, ha='left', va='bottom', weight='normal')
            elif textoffset[i] == 2.:
                x, y = m(llon[i]+0.3, llat[i]-0.3)
                plt.text(x, y, loc, size=15, ha='left', va='top', weight='normal')
            elif textoffset[i] == 3.:
                x, y = m(llon[i]-0.3, llat[i]-0.2)
                plt.text(x, y, loc, size=15, ha='right', va='top', weight='normal')
                
        # plt letter
        xlim = ax.get_xlim()
        xtxt = xlim[1] * 0.02
        ylim = ax.get_ylim()
        ytxt = ylim[1] * 0.02
        plt.text(xtxt, ytxt, pltlett[j], fontsize=24, va='bottom', ha='left')
        
        ##########################################################################################
        # plot contours
        ##########################################################################################
        if addContours == True:
            x, y = m(xs, ys)
            levels = linspace(vmin, vmax, ncolours)
                
            '''
            elif probability == '2%':
                #levels = arange(0.05, 0.31, 0.05)
                levels = bounds[1:]
            '''
            print(levels)
            csm = plt.contour(x, y, resampled, levels, colors='0.2', lw=0.3)    
            #csm_lo = plt.contour(x, y, resampled, levels_lo, colors='0.2', lw=0.3)
            
            plt.clabel(csm, inline=1, fontsize=10, fmt='%0.2f')
            #plt.clabel(csm_lo, inline=1, fontsize=10, fmt='%0.3f')
        
        j += 1
        
        """
        ##########################################################################################
        # add GA logo
        ##########################################################################################
        # load logo
        map_bbox = ax.get_position().extents
        try:
            im = plt.imread('../GAlogo.png')
        except:
            # cover all bases
            try:
                im = plt.imread('/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/postprocessing/GAlogo.png')
            except:
                try:
                    im = plt.imread('/short/w84/NSHA18/sandpit/tia547/NSHA2018/postprocessing/GAlogo.png')
                except:
                    im = plt.imread('/Users/tallen/Documents/Geoscience_Australia/NSHA2018/postprocessing/GAlogo.png')
        
        # set bbox for logo
        imoff = 0.02
        #logo_bbox = mpl.transforms.Bbox(array([[map_bbox[0]+imoff,map_bbox[1]+imoff],[0.15,0.15]]))
        #logo_bbox = [map_bbox[0]+0.11,map_bbox[1]-0.005,0.15,0.15]
        logo_bbox = [map_bbox[0]+0.14,map_bbox[1]-0.075,0.25,0.25]
        newax = figure.add_axes(logo_bbox) #, zorder=-1)
        newax.imshow(im)
        newax.axis('off')
        
        '''
        ##########################################################################################
        # add CC-by
        ##########################################################################################
        '''
        # load logo
        try:
            im = plt.imread('../ccby_narrow.png')
        except:
            # covering all bases again
            try:
                im = plt.imread('/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/postprocessing/ccby_narrow.png')
                
            except:
                try:
                    im = plt.imread('/short/w84/NSHA18/sandpit/tia547/NSHA2018/postprocessing/ccby_narrow.png')
                    
                except:
                    im = plt.imread('/Users/tallen/Documents/Geoscience_Australia/NSHA2018/postprocessing/ccby_narrow.png')
                    
    
        # set bbox for logo
        imoff = 0.02
        logo_bbox = [map_bbox[0]+0.11,map_bbox[1]-0.005,0.2,0.2]
        logo_bbox = [0.66,map_bbox[1]-0.03,0.1,0.1]
        newax = figure.add_axes(logo_bbox) #, zorder=-1)
        newax.imshow(im)
        newax.axis('off') 
        """
        '''
        ##########################################################################################
        # superimpose area source shapefile
        ##########################################################################################
        '''
        '''
        shpfile =['..//final_inputs//SWCan_T3EclC_area1.shp', \
                  '..//final_inputs//WArctic_area.shp']
        if drawshape == True:
            for shp in shpfile:
                sf = shapefile.Reader(shp)
                drawshapepoly(m, plt, sf, col='k', lw=1.5, polyline=True)
                labelpolygon(m, plt, sf, 'CODE', fsize=14)    
        '''
        
    
'''
###########################################################################################
make colourbar
###########################################################################################
'''    

# set colourbar
plt.gcf().subplots_adjust(bottom=0.1)
cax = fig.add_axes([0.34,0.14,0.33,0.025]) # setup colorbar axes.
norm = colors.Normalize(vmin=vmin, vmax=vmax)
cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')

# set cb labels
#nticks = ceil(ncolours/2.)
linticks = logspace(vmin, vmax, nticks)
logticks = log10(linticks)
cb.set_ticks(logticks)
labels = [str('%0.2f' % 10**x) for x in logticks]
cb.set_ticklabels(labels)
cb.ax.tick_params(labelsize=16)

# set title
titlestr = ' '.join(('Ratio of', T, probability, 'in 50-Year Mean Hazard'))
cb.set_label(titlestr, fontsize=20)

# check to see if maps exists
if path.isdir('map_ratio') == False:
    mkdir('map_ratio')
    
# now save png file
print('Saving', 'map_ratio_'+modelName.replace(' ','_')+'.'+key+'.png')
plt.savefig('map_ratio_'+modelName.replace(' ','_')+'.'+key+'.png', \
            dpi=800, format='png', bbox_inches='tight')

# save pdf file
'''
plt.savefig(path.join('maps', 'hazard_map_'+modelName.replace(' ','_')+'.'+key+'.pdf'), \
            dpi=300, format='pdf', bbox_inches='tight')
'''
plt.show()
    