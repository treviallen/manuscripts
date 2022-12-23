from catalogue.parsers import parse_NSHA2018_catalogue
from io_catalogues import parse_ga_event_query
from mag_tools import get_au_ml_zone, nsha18_ml2mw, nsha18_mb2mw
from misc_tools import dictlist2array
from numpy import array, where, unique, isnan
from datetime import datetime
from os import path, mkdir, getcwd

################################################################################
# parse nsha catalogue
################################################################################

if getcwd().startswith('/nas'):
    catfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/catalogue/data/NSHA18CAT.MW.V0.1.csv'
else:
    catfile = '/Users/trev/Documents/Geoscience_Australia/NSHA2018/catalogue/data/NSHA18CAT.MW.V0.1.csv'

# parse catalogue
cat = parse_NSHA2018_catalogue(catfile)

# get data arrays
nsha18_mw = dictlist2array(cat, 'prefmag')
mw_src = dictlist2array(cat, 'mw_src')
ml_region = dictlist2array(cat, 'ml_region')
nsha18_auth = dictlist2array(cat, 'auth')
nsha18_lats = dictlist2array(cat, 'lat')
nsha18_lons = dictlist2array(cat, 'lon')
nsha18_deps = dictlist2array(cat, 'dep')
nsha18_dt = dictlist2array(cat, 'datetime')

################################################################################
# parse updated GA catalogue
################################################################################

ga_csv = 'earthquakes_export_1972-2022.csv'

# parse  csv
cat = parse_ga_event_query(ga_csv)
neac_mag = dictlist2array(cat, 'mag')
neac_magType = dictlist2array(cat, 'magType')
neac_lons = dictlist2array(cat, 'lon')
neac_lats = dictlist2array(cat, 'lat')
neac_deps = dictlist2array(cat, 'dep')
neac_dt = dictlist2array(cat, 'datetime')

# get ml zone
print('Getting MLa zone...')
neac_ml_zone = get_au_ml_zone(neac_lons, neac_lats)

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

for i in range(0, len(neac_mag)):
    # check if in date & mag range
    if neac_dt[i] > nsha18_dt[-1] and neac_mag[i] >= 2.0 and neac_ml_zone[i] != '':
        # check if in mag polys
        new_dt.append(neac_dt[i])
        new_lons.append(neac_lons[i])
        new_lats.append(neac_lats[i])
        new_deps.append(neac_deps[i])
        
        # get pref mags
        new_mags.append(neac_mag[i])
        new_magTypes.append(neac_magType[i])
        
        # convert to mw if necessary
        if neac_magType[i] == 'ML':
            new_mw.append(nsha18_ml2mw(neac_mag[i]))
        elif neac_magType[i] == 'mb':
            new_mw.append(nsha18_mb2mw(neac_mag[i]))
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
    evpath1 = path.join('data', evid)
    evpath2 = path.join('data', evid, 'current')
    if not path.exists(evpath2):
        mkdir(evpath1)
        mkdir(evpath2)
        
    return evpath1, evpath2

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
sm_path = []

# make data foler if not exist
if not path.exists('data'):
    mkdir('data')

# loop through NSHA cat and make inputs    
for i in range(0, len(nsha18_dt)):
    if nsha18_dt[i] >= datetime(1972, 1, 1) and nsha18_dt[i] < datetime(2022, 1, 1) \
       and nsha18_mw[i] >= 4.25 and ml_region[i].startswith('Other') == False:
        
        # get event ID & make folders
        evpath1, evpath = make_event_folder(nsha18_dt[i])
            
        # now make shakemap txt - skip small TC events
        if nsha18_lats[i] > -20.206 and nsha18_lats[i] < -19.657 \
           and nsha18_lons[i] > 133.475 and nsha18_lons[i] < 134.112 and nsha18_mw[i] < 4.8:
            print('Skipping TC event: ' + str(nsha18_mw[i]))
        else:
            if isnan(nsha18_deps[i]) or nsha18_deps[i] > 35.:
                nsha18_deps[i] = 10.
                
            make_shakemap_xml(nsha18_dt[i], nsha18_lats[i], nsha18_lons[i], nsha18_mw[i], nsha18_deps[i], evpath)
        
        sm_dt.append(nsha18_dt[i])
        sm_lats.append(nsha18_lats[i])
        sm_lons.append(nsha18_lons[i])
        sm_deps.append(nsha18_deps[i])
        sm_mags.append(nsha18_mw[i])
        sm_path.append(evpath1)
                            
# loop through new NEAC events    
for i in range(0, len(new_dt)):
    if new_mw[i] >= 4.25:
        
        # get event ID & make folders
        evpath1, evpath = make_event_folder(new_dt[i])
            
        # now make shakemap txt
        if new_lats[i] > -20.206 and new_lats[i] < -19.657 \
           and new_lons[i] > 133.475 and new_lons[i] < 134.112 and new_mw[i] < 4.8:
            print('Skipping TC event: ' + str(new_mw[i]))
        else:
            make_shakemap_xml(new_dt[i], new_lats[i], new_lons[i], new_mw[i], new_deps[i], evpath)
        
        sm_dt.append(new_dt[i])
        sm_lats.append(new_lats[i])
        sm_lons.append(new_lons[i])
        sm_deps.append(new_deps[i])
        sm_mags.append(new_mw[i])
        sm_path.append(evpath1)
        
# write list of shakemaps
smtxt = 'DATETIME,LON,LAT,DEP,MW,PATH,FAULT_FILE\n'
for i in range(0, len(sm_dt)):
    smtxt += ','.join((sm_dt[i].strftime('%Y-%m-%dT%H:%M:%SZ'), str(sm_lons[i]), str(sm_lats[i]), \
                       str(sm_deps[i]), str(sm_mags[i]), sm_path[i])) + '\n'
                       
# write to file
csvfile = 'shakemap_event_list.csv'
f = open(csvfile, 'w')
f.write(smtxt)
f.close()
