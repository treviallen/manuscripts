from catalogue.parsers import parse_NSHA2018_catalogue
from io_catalogues import parse_iscgem, parse_ga_event_query
from mag_tools import get_au_ml_zone, nsha18_ml2mw, nsha18_mb2mw
from misc_tools import dictlist2array
from numpy import array, where, unique
from datetime import datetime
from os import path, mkdir, getcwd

def check_banda_zone(eqlos, eqlas):
    '''
    eqlos & eqlas are a list of lats/lons
    
    returns list of regions
    '''
    
    import shapefile
    from shapely.geometry import Point, Polygon
    from mapping_tools import get_field_data
    from os import getcwd
    from numpy import array
    
    shpfile = 'shapefiles/isc_gem_selection.shp'
    
    sf = shapefile.Reader(shpfile)
    shapes = sf.shapes()
    polygons = []
    for poly in shapes:
        polygons.append(Polygon(poly.points))
        
    # loop thru events and polygons
    banda_region = []
    for lo, la in zip(eqlos, eqlas):
        in_zone = 0
            
        for poly in polygons:
            pt = Point(lo, la)
            if pt.within(poly):
                in_zone = 1
                
        banda_region.append(in_zone)
            
    return array(banda_region)
################################################################################
# parse nsha catalogue
################################################################################

if getcwd().startswith('/nas'):
    catfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/catalogue/data/NSHA18CAT.MW.V0.1.csv'
else:
    catfile = '/Users/trev/Documents/Catalogues/isc-gem/isc-gem-cat.csv'

# parse catalogue
cat = parse_iscgem(catfile)

# get data arrays
iscgem_mw = dictlist2array(cat, 'mw')
iscgem_auth = dictlist2array(cat, 'mo_auth')
iscgem_lats = dictlist2array(cat, 'lat')
iscgem_lons = dictlist2array(cat, 'lon')
iscgem_deps = dictlist2array(cat, 'dep')
iscgem_dt = dictlist2array(cat, 'datetime')

# check if in Banda Polygon
iscgem_banda_reg = check_banda_zone(iscgem_lons, iscgem_lats)

################################################################################
# parse updated GA catalogue
################################################################################

ga_csv = 'earthquakes_export_2018-2022_banda.csv'

# parse  csv
cat = parse_ga_event_query(ga_csv)
neac_mag = dictlist2array(cat, 'mag')
neac_magType = dictlist2array(cat, 'magType')
neac_lons = dictlist2array(cat, 'lon')
neac_lats = dictlist2array(cat, 'lat')
neac_deps = dictlist2array(cat, 'dep')
neac_dt = dictlist2array(cat, 'datetime')

# get ml zone
print('Getting Banda zone...')
neac_banda_reg = check_banda_zone(neac_lons, neac_lats)

################################################################################
# build arrays to append
################################################################################

new_lons = []
new_lats = []
new_deps = []
new_mags = []
new_mw = []
new_magTypes = []
new_magreg = []
new_dt = []
new_reg = []

for i in range(0, len(neac_mag)):
    # check if in date & mag range
    if neac_dt[i] > iscgem_dt[-1] and neac_mag[i] >= 2.0:
        # check if in mag polys
        new_dt.append(neac_dt[i])
        new_lons.append(neac_lons[i])
        new_lats.append(neac_lats[i])
        new_deps.append(neac_deps[i])
        new_reg.append(neac_banda_reg[i])
        
        # get pref mags
        new_mags.append(neac_mag[i])
        new_magTypes.append(neac_magType[i])
        
        # convert to mw if necessary
        if neac_magType[i] == 'ML':
            new_mw.append(iscgem_ml2mw(neac_mag[i]))
        elif neac_magType[i] == 'mb':
            new_mw.append(iscgem_mb2mw(neac_mag[i]))
        else:
            # assume mw of some type
            new_mw.append(neac_mag[i])
        
print(unique(array(new_magTypes)))

################################################################################
# make shakemap inputs
################################################################################
def make_event_folder(dt):
    from os import path, mkdir
    
    # get event ID & make folders
    evid = dt.strftime('%Y%m%d%H%M')
    evpath1 = path.join('banda_data', evid)
    evpath2 = path.join('banda_data', evid, 'current')
    if not path.exists(evpath2):
        #mkdir('banda_data')
        mkdir(evpath1)
        mkdir(evpath2)
        
    return evpath2

def make_shakemap_xml(dt, lat, lon, mag, dep, evpath):
    # read template
    lines = open('event_xml_template.txt').readlines()

    # now make shakemap txt
    newlines = ''
    for line in lines[0:-1]:
        newlines += line
    
    newline = ' '.join(('<earthquake id="'+dt.strftime('%Y%m%d%H%M')+'"',
                        'lat="'+str(lat)+'"',
                        'lon="'+str(lon)+'"',
                        'mag="'+str(mag)+'"',
                        'time="'+dt.strftime('%Y-%m-%dT%H:%M:%SZ')+'"',
                        'timezone="GMT"',
                        'depth="'+str(dep)+'"',
                        'locstring="" description="" type="RS" created="1056046245" netid="au" />'))
    
    newlines += newline
                        
    # export xml
    xmlpath = path.join(evpath, 'event.xml')
    f = open(xmlpath, 'w')
    f.write(newlines)
    f.close()
    
'''
<earthquake id="201206191053" lat="-38.259" lon="146.290" mag="5.17" year="2012" 
month="06" day="19" hour="10" minute="53" second="29" timezone="GMT" depth="17.2" 
locstring="Moe, Victoria" description="Moe Mean" type="RS" created="1056046245" network="au" />
'''      

sm_dt   = []
sm_lats = []
sm_lons = []
sm_deps = []
sm_mags = []

# make data foler if not exist
if not path.exists('data'):
    mkdir('data')

# loop through NSHA cat and make inputs    
for i in range(0, len(iscgem_dt)):
    if iscgem_dt[i] >= datetime(1972, 1, 1) and iscgem_dt[i] < datetime(2019, 1, 1) \
       and iscgem_mw[i] >= 6.5 and iscgem_banda_reg[i] == 1:
        
        # get event ID & make folders
        evpath = make_event_folder(iscgem_dt[i])
            
        # now make shakemap txt - skip small TC events
        make_shakemap_xml(iscgem_dt[i], iscgem_lats[i], iscgem_lons[i], iscgem_mw[i], iscgem_deps[i], evpath)
        
        sm_dt.append(iscgem_dt[i])
        sm_lats.append(iscgem_lats[i])
        sm_lons.append(iscgem_lons[i])
        sm_deps.append(iscgem_deps[i])
        sm_mags.append(iscgem_mw[i])
    
    # M 6.0 in arafura sea
    if iscgem_dt[i] >= datetime(2000, 12, 23) and iscgem_dt[i] < datetime(2000, 12, 24) \
       and iscgem_mw[i] >= 5.8 and iscgem_banda_reg[i] == 1:
        
        # get event ID & make folders
        evpath = make_event_folder(iscgem_dt[i])
            
        # now make shakemap txt - skip small TC events
        make_shakemap_xml(iscgem_dt[i], iscgem_lats[i], iscgem_lons[i], iscgem_mw[i], iscgem_deps[i], evpath)
        
        sm_dt.append(iscgem_dt[i])
        sm_lats.append(iscgem_lats[i])
        sm_lons.append(iscgem_lons[i])
        sm_deps.append(iscgem_deps[i])
        sm_mags.append(iscgem_mw[i])
                            
# loop through new NEAC events    
for i in range(0, len(new_dt)):
    if new_mw[i] >= 6.5 and new_reg[i] == 1:
        
        # get event ID & make folders
        evpath = make_event_folder(new_dt[i])
            
        # now make shakemap txt
        make_shakemap_xml(new_dt[i], new_lats[i], new_lons[i], new_mw[i], new_deps[i], evpath)
        
        sm_dt.append(new_dt[i])
        sm_lats.append(new_lats[i])
        sm_lons.append(new_lons[i])
        sm_deps.append(new_deps[i])
        sm_mags.append(new_mw[i])
        
# write list of shakemaps
smtxt = 'DATETIME,LON,LAT,DEP,MW\n'
for i in range(0, len(sm_dt)):
    smtxt += ','.join((sm_dt[i].strftime('%Y-%m-%dT%H:%M:%SZ'), str(sm_lons[i]), str(sm_lats[i]), \
                       str(sm_deps[i]), str(sm_mags[i]))) + '\n'
                       
# write to file
csvfile = 'shakemap_banda_event_list.csv'
f = open(csvfile, 'w')
f.write(smtxt)
f.close()
