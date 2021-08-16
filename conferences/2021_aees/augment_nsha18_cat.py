from catalogue.parsers import parse_NSHA2018_catalogue
from io_catalogues import parse_ga_event_query
from mag_tools import get_au_ml_zone, nsha18_ml2mw, nsha18_mb2mw
from misc_tools import dictlist2array
from numpy import array, where, unique

################################################################################
# parse nsha catalogue
################################################################################

catfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/catalogue/data/NSHA18CAT.MW.V0.1.csv'

# parse catalogue
cat = parse_NSHA2018_catalogue(catfile)

# get data arrays
mw_pref = dictlist2array(cat, 'prefmag')
mw_src = dictlist2array(cat, 'mw_src')
ml_region = dictlist2array(cat, 'ml_region')
nsah18_auth = dictlist2array(cat, 'auth')
nsah18_lats = dictlist2array(cat, 'lat')
nsah18_lons = dictlist2array(cat, 'lon')
nsah18_deps = dictlist2array(cat, 'dep')
nsah18_dt = dictlist2array(cat, 'datetime')

################################################################################
# parse updated GA catalogue
################################################################################

ga_csv = 'earthquakes_export_from_1840_sorted.csv'

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
new_mags_mw = []
new_magTypes = []
new_magreg = []
new_dt = []

for i in range(0, len(neac_mag)):
    # check if in date & mag range
    if neac_dt[i] > nsah18_dt[-1] and neac_mag[i] >= 2.0 and neac_ml_zone[i] != '':
        # check if in mag polys
        new_dt.append(neac_dt[i])
        new_lons.append(neac_lons[i])
        new_lats.append(neac_lats[i])
        new_deps.append(neac_deps[i])
        new_dt.append(neac_dt[i])
        
        # get pref mags
        new_mags.append(neac_mag[i])
        new_magTypes.append(neac_magType[i])
        
        # convert to mw if necessary
        if neac_magType[i] == 'ML':
            new_mags_mw.append(nsha18_ml2mw(neac_mag[i]))
        elif neac_magType[i] == 'mb':
            new_mags_mw.append(nsha18_mb2mw(neac_mag[i]))
        else:
            # assume mw of some type
            new_mags_mw.append(neac_mag[i])
        
print(unique(array(new_magTypes)))