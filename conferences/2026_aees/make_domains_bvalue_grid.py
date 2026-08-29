
"""
# get neotectonic domain number and Mmax from zone centroid
def get_simple_neotectonic_domain_params(plon, plat, refShpFile):
    import shapefile
    from shapely.geometry import Point, Polygon
    from mapping_tools import get_field_data #, get_shp_centroid
    
    # load domains shp
    #print(refShpFile)
    domshp = shapefile.Reader(refShpFile)
    
    # get domains
    #neo_doms = get_field_data(dsf, 'DOMAIN', 'float')
    neo_mmax = get_field_data(domshp, 'MMAX_BEST', 'float')
    neo_bval = get_field_data(domshp, 'BVAL_BEST', 'float')
    neo_bval_l = get_field_data(domshp, 'BVAL_LOWER', 'float')    
    neo_mmax = get_field_data(domshp, 'MMAX_BEST', 'float')  
    
    # get bval sigma
    bval_sig = neo_bval_l - neo_bval
    
    # get domain polygons
    dom_shapes = domshp.shapes()
    domain = []
    min_rmag = []
    mmax = []
    bval = []
    bval_sig = []
    
    # loop through domains and find point in poly
    matchidx = -99
    point = Point(plon, plat)
    for i in range(0, len(dom_shapes)):
        # make sure trts match
        dom_poly = Polygon(dom_shapes[i].points)
            
        # check if target centroid in domains poly
        if point.within(dom_poly):
            matchidx = i
            
    #min_rmag = neo_min_reg[matchidx]
    mmax = neo_mmax[matchidx]
    trt = neo_trt[matchidx]
    bval = neo_bval[matchidx]
    bval_sig = bval_sig[matchidx]
    
    return mmax, bval, bval_sig
"""
from numpy import arange, nan
import shapefile
from os import path
from get_bvalue import get_simple_neotectonic_domain_params
    
#refShpFile = path.join('shapefiles','Domains_NSHA23_MFD.shp')
refShpFile = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/zones/2026_testing/domains_2026amt/shapefiles/Domains_2026AMT_MFD.shp'
domshp = shapefile.Reader(refShpFile)

#mmax, bval, bval_sig = get_simple_neotectonic_domain_params(148, -30, domshp)    

spacing = 0.1
grid_lims = {"xmin": 108.,
             "xmax": 156.,
             "ymin": -48.,
             "ymax": -10.,
             "zmin": 0.,
             "zmax": 1000.,
             "xspc": spacing,
             "yspc": spacing,
             "zspc": 1000.} 

# set 0.1 degree grid

xgrd = arange((grid_lims['xmin'] + spacing/2), grid_lims['xmax'], spacing)
ygrd = arange((grid_lims['ymin'] + spacing/2), grid_lims['ymax'], spacing)

# loop through and set neotectonic domains b-value
txt = 'longitude,latitude,bvalue\n'
for x in xgrd:
    for y in ygrd:
        bval = nan
        try:
            mmax, bval, bval_sig = get_simple_neotectonic_domain_params(x, y, domshp)
        except:
            dummy = 0
            #print('Point not matched: '+str(x)+', '+str(y))
        txt += ','.join((str('%0.4f' % x), str('%0.4f' % y), str('%0.4f' % bval))) + '\n'
        
f = open('neodomains_bval_grid.xyz', 'w')
f.write(txt)
f.close()

'''
# Now run:
> pygmt #changing conda environment
> gmt xyz2grd neodomains_bval_grid.xyz -Gneodomains_bval.grd -I0.1d -R108/156/-48/-10
> gmt surface neodomains_bval_grid.xyz -Gneodomains_bval.grd -I30c -R108/156/-48/-10
'''
  