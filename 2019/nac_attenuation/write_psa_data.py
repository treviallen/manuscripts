import pickle
import shapefile
from shapely.geometry import Point, Polygon
from mapping_tools import get_field_data, distance
from misc_tools import dictlist2array
from numpy import array, unique, log, exp, where, argsort, interp, isnan

# Function to print sum
def checkKey(dict, key):
    if key in dict.keys():
        return True
    else:
        return False

################################################################################
# save/load pickle
################################################################################
pklfile = 'stdict.pkl'
pklfile = 'stdict_ampfact.pkl'
stdict = pickle.load(open(pklfile, "rb"))


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

ts = tab_stdict
################################################################################
# save/load shp
################################################################################

evids = dictlist2array(tab_stdict, 'ev')
uevids = unique(dictlist2array(tab_stdict, 'ev')) 
rhyps = dictlist2array(tab_stdict, 'rhyp')

T = array([10.0, 7.5, 6.0, 5.0, 4.0, 3.0, 2.5, 2.0, 1.5, 1.25, 1.0, 0.75, 0.6, 0.5, 0.45, 0.4, 0.3, 0.25, 0.20, 0.15, 0.125, 0.10])[::-1] # use NGA-E periods
Tstr = ''
for t in T:
   Tstr += ',T' + str('%0.2f' % t)


tabtxt = 'ORIGIN_TIME,LON,LAT,DEPTH,MW,STA,NET,INST_TY,VS30,RHYP,AZIM,PGV,PGA'+Tstr+'\n'
tabDict = []
for ev in uevids:
    # get idx for event
    idx = where(evids == ev)[0]
    
    # sort by distance
    didx = argsort(rhyps[idx])
    
    #print(rhyps[idx][didx])
    
    # start writing txt
    for i in idx[didx]:
       rec = ts[i]
       # interpolate SA data
       sa_interp = exp(interp(log(T), log(rec['per']), log(rec['geom'])))
       
       sa_str = ''
       for sai in sa_interp:
           sa_str += ','+str('%0.4e' % sai)
       
       inst_ty = rec['chstr'][0][0:2]+'*'
       	
       if isnan(rec['azim']):
           azim = distance(rec['eqla'], rec['eqlo'], rec['stla'], rec['stlo'])[1]
           rec['azim'] = azim
       
       tabtxt += ','.join((rec['ev'], str('%0.3f' % rec['eqlo']), str('%0.3f' % rec['eqla']), \
                           str('%0.0f' % rec['dep']), str('%0.1f' % rec['mag']), rec['sta'], rec['net'], \
                           inst_ty, str('%0.0f' % rec['vs30']), str('%0.0f' % rec['rhyp']), \
                           str('%0.0f' % rec['azim']), str('%0.4e' % rec['pgv']), \
                           str('%0.4e' % rec['pga']))) + sa_str + '\n'
    
if pklfile.startswith('stdict_ampfact'):
    f = open('submitted/base_amp_model_flatfile.csv', 'w')
else:
    f = open('submitted/base_model_flatfile.csv', 'w')
f.write(tabtxt)
f.close()
