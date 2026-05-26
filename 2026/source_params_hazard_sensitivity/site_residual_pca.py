import pickle
import pandas as pd 
from numpy import unique, array, arange, log, log10, exp, mean, nanmean, ndarray, sqrt, \
                  nanmedian, hstack, pi, nan, isnan, interp, where, zeros_like, polyfit, \
                  hstack, nanstd, loadtxt, nanmedian, vstack
from datetime import timedelta
from obspy import UTCDateTime
from sys import argv
from get_mag_dist_terms import get_distance_term, get_magnitude_term, get_kappa_term, get_regional_term
from gmt_tools import cpt2colormap, remove_last_cmap_colour, remove_first_cmap_colour
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
from matplotlib import colorbar
import matplotlib as mpl
mpl.style.use('classic')
plt.rcParams['pdf.fonttype'] = 42

###############################################################################
# load datasets 
###############################################################################

recs = pickle.load(open('../../2023/au_stress_drop/fft_data.pkl', 'rb' ))
#recs = pickle.load(open('fft_data_mag_match.pkl', 'rb' ))
chan = recs[0]['channels'][0]
freqs = recs[0][chan]['freqs']
fidx = arange(25,len(freqs),5)
plt_freqs = freqs[fidx]

fidx = set(arange(25,len(freqs),5))

# load atten coeffs
print('../../2023/au_stress_drop/atten_coeffs_1.3_5km.pkl')
coeff_pkl = argv[1]
coeffs = pickle.load(open(coeff_pkl, 'rb' ))
print(len(coeffs))

clen = len(coeffs)

# load station sets
lines = open('../../2023/au_stress_drop/station_sets.csv').readlines()
sta_sets = []
for line in lines:
    sta_sets.append(set(line.strip().split(',')))

# get stas to ignore
ignore_stas = open('../../2023/au_stress_drop/sta_ignore_pca.txt').readlines()
ignore_stas = set([x.strip() for x in ignore_stas])

######################################################################################
# make dataset
######################################################################################
keep_nets = set(['AU', 'IU', 'S1', 'II', 'G', 'MEL', 'ME', '2O', 'AD', 'SR', 'UM', 'AB', 'VI', 'GM', 'M8', 'DU', 'WG', '4N', \
                 '1P', '1P', '2P', '6F', '7K', '7G', 'G', '7B', '4N', '7D', '', 'OZ', 'OA', 'WG', 'XX', 'AM', 'YW', '3B', '1K', \
                 '1Q', '3O', '7F', '6K', '5G', '5C'])

def get_inter_event_terms(recs, magType=1):
    '''
    magType: 1 = Brune mags
             2 = ML-MW cluster correction (MW - dMW) where dMW from mean MW - MW(conv)
             3 = Use converted magnitude
    '''
    resDict = []
    for rec in recs:
        try:
            prefmw = rec['mag']
            channel = rec['channels'][0]
            chan = recs[0]['channels'][0]
            yres = []                                    
            
            # if qual == 1, "mag" is Mwb
            if rec['net'] in keep_nets and rec['qual'] == 1:
                if not rec['sta'] in ignore_stas: 
            
                    for f, c in enumerate(coeffs):
                    
                        if f in fidx:
                            #print("Coeffs Freq = " +str('%0.3f' % c['freq']))
                                
                            freq = recs[0][chan]['freqs'][f]
                            #print("Reg Freq = " +str('%0.3f' % freq))
                            
                            if not freq == c['freq']:
                               print('\n!!!! Frequency of coefficients inconsistent with data index !!!!\n')
                               crash
                            
                            ###############################################################################
                            # parse coefs and get model prediction
                            ###############################################################################
                    
                         
                        
                            # filter by instrument type
                            addData = True
                            if channel.startswith('SH') or channel.startswith('EH'):
                                if rec[channel]['freqs'][f] < 0.9 and rec['pazfile'].endswith('s6000-2hz.paz'):
                                    addData = False
                                elif rec[channel]['freqs'][f] < 0.4:
                                    addData = False
                            
                            # filer by sample-rate
                            if rec[channel]['freqs'][f] > (0.4 * rec[channel]['sample_rate']):
                                addData = False
                            
                            if rec[channel]['sn_ratio'][f] >= 4. and addData == True:
                                
                                magterm = get_magnitude_term(prefmw, c)
                                #print(magterm)
                                
                                # get distance term
                                distterm = get_distance_term(rec['rhyp'], c)
                                
                                #	get distance independent kappa
                                kapterm = get_kappa_term(rec['sta'], c['freq'])
                                
                                #	get regional term
                                regterm = get_regional_term(rec['rhyp'], c, rec['eqdom'])
                                
                                # get total correction
                                ypred = magterm + distterm + kapterm + regterm
                                #print(ypred)
                            
                                yobs = log10(rec[channel]['p-swave_spec'][f])
                                yres.append(yobs - ypred)
                            
                            else:
                                yres.append(nan)
                    
                    resData = {'yres':array(yres), 'mag':rec['mag'], 'ev':rec['ev'], 'sta':rec['sta'], \
                               'rhyp':rec['rhyp'], 'evdt':rec['evdt']}
                    
                    resDict.append(resData)
        
        except:
            dummy = 0
                    
    return resDict


print('Getting raw residual data ...')
'''
resDict = get_inter_event_terms(recs)
pklfile = open('site_residual_data.pkl', 'wb')
pickle.dump(resDict, pklfile, protocol=-1)
pklfile.close()
'''
# load residuals
resDict = pickle.load(open('site_residual_data.pkl', 'rb' ))

# make data array
data = []
stas = []
for rd in resDict:
    # find occurrence of nan values
    nanskip = False
    idx = where(isnan(rd['yres']))[0]
    if len(idx) > 10:
        nanskip == True
    
    if len(rd['yres']) > 0 and nanskip == False and rd['rhyp'] <= 500.:
        stas.append(rd['sta'])
        if len(data) == 0:
            data = rd['yres']
        else:
            data = vstack((data, rd['yres']))

######################################################################################
# now get station average
######################################################################################
nandat = arange(0,25) * nan

ustas = unique(stas)
sta_means = []
for us in ustas:
    # add dummy data that will be ignored
    sta_data = nandat
    
    # check station sets
    sta_set = set([us])
    
    for ss in sta_sets:
        if us in ss:
            sta_set = ss
    
    for d, sta in zip(data, stas):
        if sta in sta_set:
            sta_data = vstack((sta_data, d))
    
    print(us, sta_data.shape)            
    # get station mean
    sta_mean = nanmean(sta_data, axis=0)
    if len(sta_means) == 0:            
        sta_means = sta_mean
    else:
        sta_means = vstack((sta_means, sta_mean))
######################################################################################
# now do PCA
######################################################################################
# set target site classes
n_clusters = 6
targets = arange(1,n_clusters+1)
target_names = ['SC'+str(x) for x in targets]

#fix nan values
from sklearn.impute import SimpleImputer
imputer = SimpleImputer(missing_values = nan, strategy='mean')
imputer = imputer.fit(sta_means)
data_imp =imputer.transform(sta_means)

# Apply PCA with two components (for 2D visualization)
pca = PCA(n_components=2)
x_scaled = StandardScaler().fit_transform(data_imp)

pca_features = pca.fit_transform(x_scaled)

################################################################################
# do clustering
################################################################################

print('Starting Cluster Analysis ...')
inertias = []

n = 16
for i in range(1,n):
    print('    N cluster = '+str(i))
    kmeans = KMeans(n_clusters=i)
    kmeans.fit(pca_features)
    inertias.append(kmeans.inertia_)

plt.clf()
fig = plt.figure(1, figsize=(5,5))
plt.plot(range(1,n), inertias, marker='o')
plt.title('Elbow method')
plt.xlabel('Number of clusters')
plt.ylabel('Inertia')
plt.savefig('figures/elbow.png', bbox_inches='tight')
plt.show() 
plt.clf()

################################################################################
# do actual cluster
################################################################################
# Plot the results
plt.clf()

fig = plt.figure(2, figsize=(6,6))
plt.cla()
cptfile = '//Users//trev//Documents//DATA//GMT//cpt//keshet.cpt'

cmap, zvals = cpt2colormap(cptfile, n_clusters+1, rev=True)
cmap = remove_first_cmap_colour(cmap)

cols = (cmap(arange(n_clusters)))

# plot clusters
#n_clusters = 10
kmeans = KMeans(n_clusters=n_clusters, random_state=0, max_iter=1000)
kmeans.fit(pca_features)
clusters = kmeans.labels_

for i in range(0, n_clusters):
    idx = where(kmeans.labels_ == i)[0]
    #print(kmeans.labels_[idx])
    plt.plot(pca_features[:, 0][idx], pca_features[:, 1][idx], 'o', c=cols[i], ms=6, label=target_names[i])
    print('Cluster '+str(i))
  
'''
plt.scatter(pca_features[:, 0], pca_features[:, 1], c=kmeans.labels_, s=30, cmap='viridis', \
	          vmin=0, vmax=n_clusters-1, edgecolor='k')
'''
# add polygons
#https://stackoverflow.com/questions/60913394/can-there-be-overlap-in-k-means-clusters

#plt.title('PCA of Iris Dataset')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.legend(fontsize=8, loc=2, numpoints=1, ncol=2)
plt.savefig('figures/pca_scatter_plot.png',dpi=300,bbox_inches='tight')	

plt.show()

################################################################################
# write clusters
################################################################################

clust_txt = 'STA,CLUSTER\n'

for us, cluster in zip(ustas, clusters):
    clust_txt += ','.join((us, str(cluster+1))) + '\n'
    
f = open('site_pca_cluster.csv', 'w')
f.write(clust_txt)
f.close()


################################################################################
# now plot spectra
################################################################################
fig = plt.figure(3, figsize=(8,3))
plt.cla()
# loop through clusters
class_means = []
for t, tn in zip(targets, target_names):
    class_data = nandat
    for c, sm in zip(clusters, sta_means):
        if c == t-1:
            plt.semilogx(plt_freqs, sm, '-', c=cols[t-1], lw=0.25, alpha=0.25)
            class_data = vstack((class_data, sm))
    
    class_mean = nanmean(class_data, axis=0)
    if len(class_means) == 0:            
        class_means = class_mean
    else:
        class_means = vstack((class_means, class_mean))

# plot means
i = 0
for target_name, class_mean in zip(target_names, class_means):
    plt.semilogx(plt_freqs, class_mean, '-', c=cols[i], lw=2, label=target_name)
    i += 1

plt.legend(loc=2, ncol=3, fontsize=10)

plt.xlim([0.1, 25])
plt.ylim([-4,4])
plt.xlabel('Frequency (Hz)')
plt.ylabel('ln Residual')    
plt.grid(which='both')

plt.savefig('figures/clustered_site_conditions.png',dpi=300,bbox_inches='tight')	

plt.show()