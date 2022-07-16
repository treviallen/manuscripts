from numpy import arange, around, array, delete, log10, linspace, random, unique, where
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('classic')
import pickle

lines = open('sea_stn_mags.csv').readlines()[1:]

rhyps = []
mags = []
dts = []

for line in lines:
    data = line.strip().split(',')
    rhyps.append(float(data[2]))
    mags.append(float(data[-2]))
    dts.append(float(data[0]))
    
##############################################################################
# get distance distributions
##############################################################################
fig = plt.figure(1, figsize=(15,10))
plt.suptitle('log Distance Distribution')

mrng = [1.5, 2.25, 2.75, 3.25, 3.75, 4.25, 6.]
bins = arange(0.45, 3, 0.1)

txt = 'minmag,mu,std\n'

# loop thru mrng
distPDF = []
for i in range(0, len(mrng)-1):
    log_rhyps = []
    for mag, rhyp in zip(mags, rhyps):
        if mag >= mrng[i] and mag < mrng[i+1]:
            log_rhyps.append(log10(rhyp))
    
    log_rhyps = array(log_rhyps)        
    
    # plot histogram 
    ax = plt.subplot(2,3,i+1)
    plt.hist(log_rhyps, bins=bins, color='0.75', density=True)
    plt.xlim((0.5, 3.))
    
    # Fit a normal distribution to the data:
    ymin, ymax = plt.ylim()
    mu, std = norm.fit(log_rhyps)
    
    # Plot the PDF.
    xmin, xmax = plt.xlim()
    x = linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)
    plt.plot(x, p, 'k', linewidth=2)
    title = 'M'+str(mrng[i])+'-'+str(mrng[i+1])+": mu = %.2f,  std = %.2f" % (mu, std)
    plt.title(title)
    plt.xlabel('log Ryhp (km)')
    plt.ylabel('PDF')
    
    txt += ','.join((str(mrng[i]), str(mu), str(std))) + '\n'
    
    distPDF.append({'mmax': mrng[i], 'mmin': mrng[i+1], 'mu':mu, 'std':std})

# write to file
f = open('mu_std_distance.csv', 'w')
f.write(txt)
f.close()

# save pickle
pickle.dump(distPDF, open( "distPDF.pkl", "wb" ))

plt.savefig('mag-dependent_log_Rhyp_dist.png', fmt='png', bbox_inches='tight')   
plt.show()
    
##############################################################################
# get N stations distributions
##############################################################################
fig = plt.figure(2, figsize=(15,10))
uevs, nstas = unique(array(dts), return_counts=True)

# get evmag
evmags = []
for ue in uevs:
    idx = where(ue == dts)[0]
    evmags.append(mags[idx[0]])
	
txt = 'minmag,mu,std\n'

bins = arange(4, 30, 2)
numPDF = []
# loop thru mrng
for i in range(0, len(mrng)-1):
    sta_dat = []
    for uev, nsta, evmag in zip(uevs, nstas, evmags):
        if evmag >= mrng[i] and evmag < mrng[i+1]:
            sta_dat.append(nsta)
    sta_dat = array(sta_dat)
           
    # plot histogram 
    ax = plt.subplot(2,3,i+1)
    plt.hist(sta_dat, bins=bins, color='0.75', density=True)
    #plt.xlim((0.5, 3.))
    
    # Fit a normal distribution to the data:
    ymin, ymax = plt.ylim()
    mu, std = norm.fit(sta_dat)
    
    # Plot the PDF.
    xmin, xmax = plt.xlim()
    x = linspace(xmin, xmax, 100)
    p = norm.pdf(x, mu, std)
    plt.plot(x, p, 'k', linewidth=2)
    title = 'M'+str(mrng[i])+'-'+str(mrng[i+1])+": mu = %.2f,  std = %.2f" % (mu, std)
    plt.title(title)
    plt.xlabel('# Stations')
    plt.ylabel('PDF')
    
    plt.suptitle('N Stas')
    
    txt += ','.join((str(mrng[i]), str(mu), str(std))) + '\n'
    
    numPDF.append({'mmax': mrng[i], 'mmin': mrng[i+1], 'mu':mu, 'std':std})

# write to file
f = open('mu_std_nsta.csv', 'w')
f.write(txt)
f.close()

#save pickle
pickle.dump(numPDF, open( "numPDF.pkl", "wb" ))
plt.savefig('mag-dependent_station_number.png', fmt='png', bbox_inches='tight')   
plt.show()

##############################################################################
# for given mag, sample distance distribution
##############################################################################

mx = 2.25
#x = linspace(0, 40, 40)

# get n stasions
for npdf in numPDF: 
    if npdf['mmin'] >= mx and npdf['mmax'] < mx:
        print(npdf)
        sample_nsta = around(norm.rvs(npdf['mu'], npdf['std'], size=1))
        
# now get distances for Nstas
for dpdf in distPDF: 
    if dpdf['mmin'] >= mx and dpdf['mmax'] < mx:
        print(dpdf)
        sample_dists = 10**(norm.rvs(dpdf['mu'], dpdf['std'], size=int(sample_nsta)))
        
print(sample_dists)
print('\n')

# strip stas > 800 km
idx = sample_dists > 800.
sample_dists = delete(sample_dists, idx)
print(sample_dists)