import pickle
import shapefile
from shapely.geometry import Polygon, Point
import matplotlib.pyplot as plt
import matplotlib as mpl
from mapping_tools import get_field_data

mpl.style.use('classic')

def return_mag_dist(recs, polygons):
    maxdist = 1750.
    for j, poly in enumerate(polygons):
        
        if zone_code[j] == 'BS':
            mag = []
            rhyp = []
                    
            for rec in recs: 
                if rec['rhyp'] <= maxdist:
                    pt = Point(rec['eqlo'], rec['eqla'])
                    if pt.within(poly):
                        mag.append(rec['mag'])
                        rhyp.append(rec['rhyp'])
    
    return mag, rhyp

##########################################################################################
# load shapefile
##########################################################################################

shpfile = 'shapefiles/nac_gmm_zones.shp'
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
zone_code = get_field_data(sf, 'CODE', 'str')
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))

##########################################################################################
# parse amp SA files
##########################################################################################
# reads sa data files and returns period (T) and acceleration (SA) vectors
print('Loading pkl file...')
recs = pickle.load(open("stdict_ampfact.pkl", "rb" ))

amag, arhyp = return_mag_dist(recs, polygons)

##########################################################################################
# parse base SA files
##########################################################################################
# reads sa data files and returns period (T) and acceleration (SA) vectors
print('Loading pkl file...')
recs = pickle.load(open("stdict.pkl", "rb" ))

bmag, brhyp = return_mag_dist(recs, polygons)

##########################################################################################
# plot
##########################################################################################

fig = plt.figure(1, figsize=(8,8))
plt.plot(arhyp, amag, 'o', mec='orangered', mfc='none', mew=1.5, ms=7, label='Site Only')
plt.plot(brhyp, bmag, 'o', mfc='dodgerblue', mec='none', ms=7, label='Base Model')

plt.xlabel('Hypocentral Distance (km)', fontsize=18)
plt.ylabel('Moment Magnitude', fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid(which='both')

plt.legend(loc=2, fontsize=16, numpoints=3)

plt.savefig('figures/2021_banda_mag_dist.png', fmt='png', dpi=300, bbox_inches='tight')
plt.show()
