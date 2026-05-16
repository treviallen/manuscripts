#import matplotlib as mpl
#mpl.style.use('classic')
import matplotlib.pyplot as plt
from numpy import array, arange, nan, zeros_like, where, nanmean, nanstd, loadtxt, unique, log10
import matplotlib as mpl
mpl.style.use('classic')
plt.rcParams['pdf.fonttype'] = 42
import warnings
warnings.filterwarnings("ignore")

################################################################################

# load Brune data
csvfile = 'brune_stats_cluster.csv'
#data = loadtxt(csvfile, delimiter=',', skiprows=1)

print(csvfile)
lines = open(csvfile).readlines()
lats = []
lons = []
mw = []
qual = []
stressdrops = []
brunedate = []
cluster = []

for line in lines[1:]:
    dat = line.strip().split(',')
    brunedate.append(dat[0])
    lons.append(float(dat[2]))
    lats.append(float(dat[3]))
    mw.append(float(dat[8]))
    qual.append(float(dat[-5]))
    cluster.append(float(dat[-1]))
    stressdrops.append(float(dat[10]))
    
lons = array(lons)
lats = array(lats)
stressdrops = array(stressdrops)
qual = array(qual)
mw = array(mw)
cluster = array(cluster)
n = len(unique(cluster))

################################################################################
# load ml data
'''
EVENT,EVID,ML_2800,ML_2800_STD,ML_2080,ML_2080_STD,MW,SD,REG,NRECS
1989-12-27T23:26:57.000000Z,,5.528,0.328,5.446,0.328,nan,nan,E,10
'''
# load Brune data
csvfile = '../../2023/au_stress_drop/ml_mw_stats.csv'

print(csvfile)
lines = open(csvfile).readlines()
ml = []
mldate = []

for line in lines[1:]:
    dat = line.strip().split(',')
    mldate.append(dat[0])
    ml.append(float(dat[2]))

ml = array(ml)

################################################################################
# lmatch events
ml_match = zeros_like(mw) * nan

for i, bd in enumerate(brunedate):
    for j, mld in enumerate(mldate):
        if bd == mld:
            print(bd)
            ml_match[i] = ml[j]

################################################################################
# load NSHA23 ML-MW 2800 coeffs

xplt = arange(3, 6.81, 0.01)
s0, s1, s2 = loadtxt('../../2023/au_stress_drop/mw-ml_coeffs_2800.csv', delimiter=',', skiprows=1)
yplt_nsha23 = s0 * xplt**2 + s1 * xplt + s2

################################################################################
# loop thru clusters and get average residual
cluster_txt = 'CLUSTER,MEAN +- STD\n'

# get mean for all events
nsha23_mw = s0 * ml_match**2 + s1 * ml_match + s2

# get mw residual
mw_res = mw - nsha23_mw

# get mean & std
mw_res_mean = nanmean(mw_res)
mw_res_std = nanstd(mw_res)
print(-1, mw_res_mean, mw_res_std) # all clusters

cluster_txt += ','.join(('-1', str(mw_res_mean), str(mw_res_std))) + '\n'

plt.cla()
plt.clf()
fig = plt.figure(1, figsize=(6, 15))
plt.rc('xtick',labelsize=8)
plt.rc('ytick',labelsize=8)

bins = arange(-1.3, 2.5, 0.2)
plt.subplot(4,3,1)
plt.hist(log10(stressdrops), bins=bins, facecolor='0.8') #, width=0.8)
plt.title('All Data', fontsize=10)
plt.ylabel('Count')
 
for i in range(0,n):
    idx = where((cluster == i) & (qual == 1))[0]
    
    # get predicted mw from rrelationship
    #print(ml_match[idx])
    nsha23_mw = s0 * ml_match[idx]**2 + s1 * ml_match[idx] + s2
    
    # get mw residual
    mw_res = mw[idx] - nsha23_mw
    
    # get mean & std
    mw_res_mean = nanmean(mw_res)
    mw_res_std = nanstd(mw_res)
    print(i, mw_res_mean, mw_res_std)
    
    cluster_txt += ','.join((str(i+1), str('%0.2f' %  mw_res_mean)+ ' +- '+str('%0.2f' % mw_res_std))) + '\n'
    
    # now plot hist
    plt.subplot(4,3,i+2)
    plt.hist(log10(stressdrops[idx]), bins=bins, facecolor='0.8') #, width=0.8)
    plt.title('Cluster '+str(i+1), fontsize=10)
    
    if i == 2 or i == 5 or i == 8:
        plt.ylabel('Count')
    if i >= 8:
        plt.xlabel('log SD')
    
f = open('ml_mw_bias_cluster.csv', 'w')
f.write(cluster_txt)
f.close()

# save fig
plt.tight_layout()
plt.savefig('figures/cluster_sd_histograms.png', format='png', dpi=300, bbox_inches='tight')       
plt.show()




