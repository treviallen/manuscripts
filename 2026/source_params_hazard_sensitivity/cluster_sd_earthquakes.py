print('!!!!! Use conda activate py311 !!!!!')

from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from numpy import array, log10, where, zeros_like, unique
import numpy as np
from shapely.geometry import Polygon, Point, MultiPoint, MultiPolygon, LineString
from shapely.ops import voronoi_diagram
from shapely import wkt
import shapefile

# load Brune data
csvfile = 'brune_stats.csv'
#data = loadtxt(csvfile, delimiter=',', skiprows=1)

print(csvfile)
lines = open(csvfile).readlines()
lats = []
lons = []
mags = []
qual = []
stressdrops = []

temptextlist = []
for line in lines[1:]:
    dat = line.strip().split(',')
    lons.append(float(dat[2]))
    lats.append(float(dat[3]))
    mags.append(float(dat[8]))
    qual.append(float(dat[-1]))
    stressdrops.append(float(dat[10]))
    
    if qual[-1] == 1:
        temptextlist.append(line)
    
lons = array(lons)
lats = array(lats)
logsd = log10(array(stressdrops))
#logsd = array(stressdrops)
sd = array(stressdrops)
qual = array(qual)
idx = where(qual == 1)[0]

################################################################################
# do clustering
print('Starting Cluster Analysis ...')
data = list(zip(logsd[idx], lons[idx], lats[idx]))
#data = list(zip(sd[idx], lons[idx], lats[idx]))
inertias = []
'''
n = 16
for i in range(1,n):
    print('    N cluster = '+str(i))
    kmeans = KMeans(n_clusters=i)
    kmeans.fit(data)
    inertias.append(kmeans.inertia_)

plt.plot(range(1,n), inertias, marker='o')
plt.title('Elbow method')
plt.xlabel('Number of clusters')
plt.ylabel('Inertia')
plt.show() 
'''

################################################################################
# standardise data

scaler = StandardScaler()
scaled_data = scaler.fit_transform(data)
#print(scaled_data)

# scale sd data
scaled_data[:,0] /= 5.

################################################################################
# plot clusters
n_clusters = 10
#kmeans = KMeans(n_clusters=9, random_state=1)
kmeans = KMeans(n_clusters=n_clusters, random_state=1)
kmeans.fit(data)
kmeans.fit(scaled_data)
plt.scatter(lons[idx], lats[idx], c=kmeans.labels_, zorder=10)

################################################################################
# export polygons
'''
https://spatial-dev.guru/2024/04/28/automated-polygon-splitting-using-voronoi-diagrams-and-clustering/
'''
xmin = 111.5
xmax = 155
ymin = -46
ymax = -5

# WKT for main polygon for which random points needs to be generated
wkt_str = 'POLYGON ((111.5 -45, 111.5 -10, 155 -10, 155 -45, 111.5 -45))'
 
polygon = wkt.loads(wkt_str)

# make points 
points = []
for lon, lat in zip(lons[idx], lats[idx]):
    point = Point(lon, lat)
    points.append(point)
    
# Converting point geometry to plain coordinates
points = [[point.x,point.y] for point in points]

# Get cluster labels
labels = kmeans.labels_
 
# Create empty lists to store points for each cluster
cluster_points = [[] for _ in range(n_clusters)]
 
# Separate points into clusters based on labels
for point_idx, label in enumerate(labels):
    cluster_points[label].append(points[point_idx])
 
# Print number of points in each cluster
for i, points in enumerate(cluster_points):
    print(f"Cluster {i+1} has {len(points)} points.")
 
# Get convex_hulls of clusters
convex_hulls = [MultiPoint(x).convex_hull for x in cluster_points]
 
# Get Centroid of convex_hulls
convex_hulls_centroids = [x.centroid for x in convex_hulls]
 
# Get voronoi for given centroids
voronoi_polygons = list(voronoi_diagram(MultiPoint([centroid for centroid in convex_hulls_centroids])).geoms)
 
# Perform Intersection with main polygon for alignment
aligned_polygon = MultiPolygon([x.buffer(0).intersection(polygon) for x in voronoi_polygons])

polys = []
for voronoi_polygon in voronoi_polygons:
    polys.append(voronoi_polygon.buffer(0).intersection(polygon))
 
print(aligned_polygon)

################################################################################
# write to shapefile

# set shapefile to write to 
outshp = 'cluster_polygons.shp'
w = shapefile.Writer(outshp)
w.field('LABELS','F', 2, 0)

ulabels = unique(labels)

# loop through fitted_contours
for poly, label in zip(polys, ulabels):
    xx, yy = poly.exterior.coords.xy
    points = []
    for x, y in zip(xx, yy):
        points.append([x, y])
    
    w.record(int(label))
    w.line(array([points]))
    
# write projection file
print(outshp)
prjfile = 'cluster_polygons.prj'
f = open(prjfile, 'w')
f.write('GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]]')
f.close()

################################################################################
# write clusters
text = lines[0].strip()+',CLUSTER\n'

for i, line in enumerate(temptextlist):
    text += line.strip() + ',' + str(kmeans.labels_[i]) + '\n'
    
f = open('brune_stats_cluster.csv', 'w')
f.write(text)
f.close()

plt.show() 

