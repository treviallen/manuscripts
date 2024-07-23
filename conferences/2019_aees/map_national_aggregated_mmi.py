def return_dyfi_data(jsonfile):
    import json


    with open(jsonfile) as f:
        data = json.load(f)

    dyfi_dict = []
    for feature in data['features']:
        tmp = {'geomerty':feature['geometry']['coordinates'][0],
               'centroid':feature['properties']['center']['coordinates'],
               'intensity':feature['properties']['intensityFine'],
               'nresp':feature['properties']['nresp']}
        
        # append to list
        dyfi_dict.append(tmp)
        
    return dyfi_dict

###############################################################################
# start main code here
###############################################################################

from io_catalogues import parse_ga_event_query
from misc_tools import listdir_extension
from numpy import array, arange, isnan, nan, unique, mean, percentile, arange, where, hstack
from os import path, getcwd
from mapping_tools import get_map_polygons, mask_outside_polygons, annotate_cities, distance

###############################################################################
# parse catalogue
###############################################################################

gacsv = 'earthquakes_export_edit.csv'
gacsv = '202406_earthquakes_export.csv'

events = parse_ga_event_query(gacsv)

###############################################################################
# grab 10-km aggregated
###############################################################################
years = arange(2018,2025)
json_folder_suffix = '_bulk_felt_reports_geojson'

# start looping by year
for y, year in enumerate(years):
    print('Making map for year: '+str(year))
    # now just keep 10-km agg
    geojson_10 = []
    geojson_years = []
    
    # loop trhu min-max years
    yearloop = 2018
    while yearloop <= year:
        print('    '+str(yearloop))
        # getting aggregated files
        all_geojson = listdir_extension(str(yearloop)+json_folder_suffix, 'geojson')
        
        for jsonfile in all_geojson:
            #if jsonfile.endswith('10km_filtered.geojson'):
            if jsonfile.endswith('5km.geojson'):
                geojson_10.append(jsonfile)
                geojson_years.append(yearloop)
        
        yearloop += 1
                
    ###############################################################################
    # loop through events and match json files
    ###############################################################################
    
    for i, event in enumerate(events):
        events[i]['dyfi'] = nan # initialise
        
        for jsonfile, geojson_year in zip(geojson_10, geojson_years):
            if jsonfile.startswith(event['event_id']):
                dyfi_dict = return_dyfi_data(path.join(str(geojson_year)+json_folder_suffix, jsonfile))
                
                events[i]['dyfi'] = dyfi_dict
    
    ###############################################################################
    # loop through events get unique centroids
    ###############################################################################
    
    cent_lonlist = []
    cent_latlist = []
    eqla = []
    eqlo = []
    eqmag = []
    cent_mmi = []
    cent_nresp = []
    event_ids = []
    total_nresp = 0
    ev_count = 0
    centroids = []
    
    for i, event in enumerate(events):
        if isinstance(event['dyfi'],list):
           # loop through dyfi recs
           for rec in event['dyfi']:
               cent_lonlist.append(rec['centroid'][0])
               cent_latlist.append(rec['centroid'][1])
               cent_mmi.append(rec['intensity'])
               cent_nresp.append(rec['nresp'])
               centroids.append(array(rec['centroid']))
               event_ids.append(event['event_id'])
               eqla.append(event['lat'])
               eqlo.append(event['lon'])
               eqmag.append(event['mag'])
               
               # add nresp
               total_nresp += rec['nresp']
               len_recs = len(event['dyfi'])
           
           if len_recs > 2:
               ev_count += 1
    
    cent_lonlist = array(cent_lonlist)
    cent_latlist = array(cent_latlist)
    cent_mmi = array(cent_mmi)
    cent_nresp = array(cent_nresp)
    event_ids = array(event_ids)
    eqlo = array(eqlo)
    eqla = array(eqla)
    eqmag = array(eqmag)
    
    # get source-site distance
    mmi_dist = []
    for ela, elo, cla, clo in zip(eqla, eqlo, cent_latlist, cent_lonlist):
        mmi_dist.append(distance(ela, elo, cla, clo)[0])
        
    mmi_dist = array(mmi_dist)
    
    '''
    # get unique values
    unique_lons = unique(array(cent_lonlist))
    unique_lats = unique(array(cent_latlist))
    unique_cent = []
    
    i = 0
    for clo, cla in zip(cent_lonlist, cent_latlist):
        cent_match = False
        if i == 0:
            unique_cent.append([clo, cla])
        else:
            for ucent in unique_cent:
                if ucent[0] == clo and ucent[1] == cla:
                    cent_match = True
                    
            if cent_match == False:
                if clo > 105 and clo < 165 and cla > -50 and cla < -5:
                    unique_cent.append([clo, cla])
            
        i += 1
    '''
    '''
    # get max mmi vals
    max_mmi_10k = []
    for ucent in unique_cent:
        max_mmi = 0.
        for cent in centroids:
            if cent == ucent:
               if cm > max_mmi:
                   max_mmi = cm
        max_mmi_10k.append(max_mmi)
    '''
    
    ###############################################################################
    # start mapping
    ###############################################################################
    
    from mpl_toolkits.basemap import Basemap
    from matplotlib import colors, colorbar
    from gmt_tools import cpt2colormap, remove_last_cmap_colour
    from mapping_tools import distance, reckon
    import matplotlib.pyplot as plt
    import matplotlib as mpl     
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    mpl.style.use('classic')
    
    ##########################################################################################
    # get cmap
    ##########################################################################################
    # get colormap
    if getcwd().startswith('/nas'):
        cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/mi_pop.cpt'
    else:
        cptfile = '//Users//trev//Documents//DATA//GMT//cpt//mi_pop.cpt'
    ncols = 10
    cmap, zvals = cpt2colormap(cptfile, ncols+1, rev=False)
    cmap = remove_last_cmap_colour(cmap)
    cs = (cmap(arange(ncols)))
    
    ##########################################################################################
    # make map
    ##########################################################################################
    
    bbox = '110.5/156.5/-44.5/-9.5'
    
    bbox = bbox.split('/')
    minlon = float(bbox[0])
    maxlon = float(bbox[1])
    minlat = float(bbox[2])
    maxlat = float(bbox[3])
    mbuff = 1.
    
    # make bounds on the fly - adjust to make square-ish
    urcrnrlat = maxlat
    llcrnrlat = minlat
    urcrnrlon = maxlon
    llcrnrlon = minlon
    
    # set up figure
    fig = plt.figure(y+1, figsize=(18,10))
    plt.tick_params(labelsize=16)
    ax = fig.add_subplot(111)
    
    # set projection
    lon_0 = mean([llcrnrlon, urcrnrlon])
    lat_1 = percentile([llcrnrlat, urcrnrlat], 25)
    lat_2 = percentile([llcrnrlat, urcrnrlat], 75)
    
    m = Basemap(projection='merc',\
                llcrnrlon=llcrnrlon,llcrnrlat=llcrnrlat, \
                urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,\
                rsphere=6371200.,resolution='i',area_thresh=1000)
    
    # draw coastlines, state and country boundaries, edge of map.
    m.drawcoastlines()
    m.drawstates()
    m.drawcountries()
    
    m.fillcontinents(color='0.9',lake_color='lightskyblue')
    m.drawparallels(arange(-90.,90.,6), labels=[1,0,0,0],fontsize=16, dashes=[2, 2], color='0.99', linewidth=0.0)
    m.drawmeridians(arange(0.,360.,10), labels=[0,0,0,1], fontsize=16, dashes=[2, 2], color='0.99', linewidth=0.0)
    
    ##########################################################################################
    # get land & lake polygons for masking
    ##########################################################################################
    
    polys = get_map_polygons(m)
    
    mask_outside_polygons(polys, 'lightskyblue', plt)
    
    # get lake ploygons and mask
    polygons = []
    for polygon in m.lakepolygons:
        poly = polygon.get_coords()
        plt.fill(poly[:,0], poly[:,1], 'lightskyblue')
        polygons.append(poly)
    
    ###############################################################################
    # make 50-km grid the ugly way
    ###############################################################################
    
    cell_size = 30.
    min_nresp = 2
    
    grd_lons = []
    grd_lats = []
    grd_geom = []
    grd_mmi = []
    grd_count = []
    lon_grd = minlon
    
    while lon_grd < maxlon:
        lat_grd = minlat
        nxt_lon, tmp_lat = reckon(lat_grd, lon_grd, cell_size*0.95, 90)
        
        while lat_grd < maxlat:
           tmp_lon, nxt_lat = reckon(lat_grd, lon_grd, cell_size, 0)
           
           # fill grids
           geom = [[lon_grd, lat_grd], [lon_grd, nxt_lat], [nxt_lon, nxt_lat], [nxt_lon, lat_grd], [lon_grd, lat_grd]]
           grd_geom.append(geom)
           
           # check if mmi data exists
           idx1 = where((cent_lonlist >= lon_grd) & (cent_lonlist < nxt_lon) & (cent_mmi >= 2.) \
                       & (cent_latlist >= lat_grd) & (cent_latlist < nxt_lat) & (cent_nresp >= min_nresp))[0]
           
           idx2 = where((cent_lonlist >= lon_grd) & (cent_lonlist < nxt_lon) & \
                        (cent_latlist >= lat_grd) & (cent_latlist < nxt_lat) & \
                        (mmi_dist <= cell_size) & (cent_mmi < 4.5) & (cent_nresp >= 1))[0]
           
           idx3 = where((cent_lonlist >= lon_grd) & (cent_lonlist < nxt_lon) & \
                        (cent_latlist >= lat_grd) & (cent_latlist < nxt_lat) & \
                        (mmi_dist <= cell_size) & (eqmag > 4.) & (cent_nresp >= 1))[0]
           
           idx2 = hstack((idx2, idx3))
           
           idx = unique(hstack((idx1, idx2)))
                       
           if len(idx) > 0:
               grd_mmi.append(max(cent_mmi[idx]))
               grd_count.append(len(idx))
           else:
               grd_mmi.append(nan)
               grd_count.append(nan)
               
           lat_grd = nxt_lat
           
        lon_grd = nxt_lon
    
    for geom, mi in zip(grd_geom, grd_mmi):
        if not isnan(mi):
        
            # map grid
            pltx = array(geom)[:,0]
            plty = array(geom)[:,1]
            
            x, y = m(pltx, plty) 
            
            colidx = int(round(mi))-1
            c= tuple(cs[colidx][:-1])
            plt.fill(x, y, fc=c, ec='0.45', lw=0.25, zorder=100)
            
           # plt.plot(x, y, '-', c='0.5', lw=0.5)   
    '''
    # map centroids
    ulos = array(unique_cent)[:,0]
    ulas = array(unique_cent)[:,1]
    x, y = m(ulos, ulas)
    plt.plot(x, y, 'b+')
    '''
    ##########################################################################################
    # annotate
    ##########################################################################################
    
    # add n responses txt
    lonrng = m.urcrnrlon - m.llcrnrlon
    latrng = m.urcrnrlat - m.llcrnrlat
    plttxt = 'Number of Responses = '+str(total_nresp)
    x, y = m(m.llcrnrlon+0.02*lonrng, m.llcrnrlat+0.98*latrng)
    
    props = dict(boxstyle='round', facecolor='w', alpha=1)
    plt.text(x, y, plttxt, size=16, ha='left', va='top', bbox=props, zorder=10000)
    
    # add cities
    capfile = 'city_names.csv'
    
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
    plt.plot(x, y, 's', markerfacecolor='None', markeredgecolor='k', markeredgewidth=0.5, markersize=8, zorder=10000)
    
    # label cities
    for i, loc in enumerate(locs):
        if textoffset[i] == 0.:
            x, y = m(llon[i]-0.35, llat[i]+0.12)
            plt.text(x, y, loc, size=15, ha='right', va='bottom', weight='normal', zorder=10000, path_effects=pe)
        elif textoffset[i] == 1.:
            x, y = m(llon[i]+0.35, llat[i]+0.12)
            #plt.text(x, y, loc, size=15, ha='left', va='bottom', weight='light', path_effects=pe)
            plt.text(x, y, loc, size=15, ha='left', va='bottom', weight='normal', zorder=10000, path_effects=pe)
        elif textoffset[i] == 2.:
            x, y = m(llon[i]+0.3, llat[i]-0.3)
            plt.text(x, y, loc, size=15, ha='left', va='top', weight='normal', zorder=10000, path_effects=pe)
        elif textoffset[i] == 3.:
            x, y = m(llon[i]-0.3, llat[i]-0.2)
            plt.text(x, y, loc, size=15, ha='right', va='top', weight='normal', zorder=10000, path_effects=pe)
    
    ##########################################################################################
    # add colorbar
    ##########################################################################################
    
    axins = inset_axes(ax,
                       width="50%",  # width = 50% of parent_bbox width
                       height="4%",  # height : 5%   
                       loc=3,
                       bbox_to_anchor=(0.01,0.062,1,1), bbox_transform=ax.transAxes)
                        
    norm = mpl.colors.Normalize(vmin=0.5, vmax=10.5)#myb
    cb = colorbar.ColorbarBase(axins, cmap=cmap, norm=norm, orientation='horizontal')
    
    # set cb labels
    ticks = range(1,11)
    rom_num = ['I', 'II', 'III', 'IV', 'V', 'VI','VII','VIII','IX','X']
    cb.set_ticks(ticks)
    cb.set_ticklabels(rom_num)
    
    titlestr = 'Macroseismic Intensity'
    cb.set_label(titlestr, fontsize=16)
    
    
    plt.savefig(str(year)+'_map_nat_agg_mmi.png', fmt='png', dpi=300, bbox_inches='tight')
plt.show()