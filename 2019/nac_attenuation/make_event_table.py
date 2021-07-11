import pickle
import shapefile
from shapely.geometry import Point, Polygon
from mapping_tools import get_field_data
from misc_tools import dictlist2array
from numpy import unique

# Function to print sum
def checkKey(dict, key):
    if key in dict.keys():
        return True
    else:
        return False

################################################################################
# save/load pickle
################################################################################

stdict = pickle.load(open("stdict_ampfact.pkl", "rb"))


################################################################################
# save/load shp
################################################################################

shpfile = 'shapefiles/2021_nac_gmm_zones.shp'
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))
    
zone_code = get_field_data(sf, 'CODE', 'str')
zone_group = get_field_data(sf, 'ZONE_GROUP', 'str')

mmin = 5.25

tab_stdict = []
for poly, zcode, zgroup in zip(polygons, zone_code, zone_group):
    for sd in stdict:
        if checkKey(sd, 'eqlo') == True:
            pt = Point(sd['eqlo'], sd['eqla'])
            if pt.within(poly) and sd['mag'] >= mmin and sd['rhyp'] <= 1750:
                tab_stdict.append(sd)

################################################################################
# save/load shp
################################################################################


evids = unique(dictlist2array(tab_stdict, 'ev'))

tabtxt = 'ORIGIN_TIME,LON,LAT,DEPTH,MW,N_RECS\n'
tabDict = []
for ev in evids:
    cnt = 0
    for sd in tab_stdict:
        if ev == sd['ev']:
            cnt += 1
            eqla = sd['eqla']
            eqlo = sd['eqlo']
            dep = sd['dep']
            mag = sd['mag']
            
    tabDict.append({'ev':ev, 'eqlo':eqlo, 'eqla':eqla, 'dep':dep, 'mw':mag, 'cnt':cnt})
    tabtxt += ','.join((ev, str('%0.3f' % eqlo), str('%0.3f' % eqla), str('%0.0f' % dep), \
                        str('%0.1f' % mag), str(cnt))) + '\n'

f = open('submitted/event_table.csv', 'w')
f.write(tabtxt)
f.close()
