import shapefile
from shapely.geometry import Polygon
from numpy import array, ones_like, nan, hstack

##########################################################################################
# parse mag zones
##########################################################################################

mlRegFile = '/Users/trev/Documents/Geoscience_Australia/NSHA2018/catalogue/magnitude/ml/australia_ml_regions.txt'
lines = open(mlRegFile).readlines()

# get WA polygons
walon = []
walat = []
for line in lines:
    if line.startswith('SWAustralia'):
        
        walon.append(float(line.split()[3]))
        walat.append(float(line.split()[2]))

walat = array(walat).reshape(len(walon),1) 
walon = array(walon).reshape(len(walon),1) 

# get SEA polygons
ealon = []
ealat = []
for line in lines:
    if line.startswith('SEAustralia'):
        
        ealon.append(float(line.split()[3]))
        ealat.append(float(line.split()[2]))
        
ealat = array(ealat).reshape(len(ealat),1) 
ealon = array(ealon).reshape(len(ealat),1) 

# get SA polygons
salon = []
salat = []
for line in lines:
    if line.startswith('SAustralia'):
        
        salon.append(float(line.split()[3]))
        salat.append(float(line.split()[2]))
        
salat = array(salat).reshape(len(salat),1) 
salon = array(salon).reshape(len(salat),1) 

##########################################################################################
# make polygons
##########################################################################################

polygons = []
#polygons.append(Polygon(hstack((walon, walat))))

# write shapefile
w = shapefile.Writer('australia_ml_regions', shapeType=5)

w.field('ML_REGION','C','5')

w.poly([hstack((walon, walat))])
w.record('WCA')
w.poly([hstack((ealon, ealat))])
w.record('EA')
w.poly([hstack((salon, salat))])
w.record('SA')

w.close()
 
# write prj file
prjfile = 'australia_ml_regions.prj'
f = open(prjfile, 'w')
f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
f.close()