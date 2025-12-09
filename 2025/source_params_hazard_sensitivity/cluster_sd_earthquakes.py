print('!!!!! Use conda activate py311 !!!!!')

from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from numpy import array, log10, where
import numpy as np


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
# plot clusters
kmeans = KMeans(n_clusters=9, random_state=1)
kmeans = KMeans(n_clusters=12, random_state=1)
kmeans.fit(data)
plt.scatter(lons[idx], lats[idx], c=kmeans.labels_)
#plt.show() 

################################################################################
# export polygons

Z_kmeans = kmeans.predict(np.c_[lons.ravel(), lats.ravel()])
Z_kmeans = Z_kmeans.reshape(lons[idx].shape)
plt.contourf(lons[idx], lats[idx], Z_kmeans, cmap=viridis, alpha=0.5)

################################################################################
# write clusters
text = lines[0].strip()+',CLUSTER\n'

for i, line in enumerate(temptextlist):
    text += line.strip() + ',' + str(kmeans.labels_[i]) + '\n'
    
f = open('brune_stats_cluster.csv', 'w')
f.write(text)
f.close()
