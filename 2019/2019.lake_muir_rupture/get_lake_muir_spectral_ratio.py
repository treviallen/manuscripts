#import obspy
import matplotlib.pyplot as plt
import numpy as np
from numpy import array,logspace,zeros_like,log10,nanmedian,where,isnan,mean,hstack,arange, nan, vstack, exp, log
from obspy import UTCDateTime
from datetime import datetime, timedelta
from data_fmt_tools import return_sta_data 
from mapping_tools import distance
from scipy import stats
from misc_tools import get_binned_stats, get_binned_stats_mean, remove_last_cmap_colour
from scipy.stats import linregress
import scipy.odr.odrpack as odrpack
from obspy import read
import pickle



#def read_datafile(file_name):
    # the skiprows keyword is for heading, but I don't know if trailing lines
    # can be specified
   # data = np.loadtxt(file_name, delimiter=str, skiprows=0)
    #return data

rec_pkl = 'wa_fft_data.pkl'
records_all = pickle.load(open(rec_pkl, "rb" ))

###############################################################################
# strip like records - not sure why they're here!
###############################################################################
records = []
for i, rec1 in enumerate(records_all):
    addRec = True
    if len(records) > 0:
        for j, rec2 in enumerate(records):
            if rec1['sta'] == rec2['sta'] and rec1['channels'][0] == rec2['channels'][0] \
                and rec1['pick_file']['datetime'] == rec2['pick_file']['datetime']:
                    addRec = False
                
    if addRec == True:
        records.append(rec1)       

###############################################################################
# match stations
###############################################################################
stacked_ratios = []
sn_thresh = 10.
for rec1 in records:
    if str(rec1['pick_file']['datetime']).startswith('2018-11-08'):
        for rec2 in records:
            if rec1['sta'] == rec2['sta'] and rec1['channels'][0] == rec2['channels'][0] \
                and str(rec2['pick_file']['datetime']).startswith('2018-10-12'):
                    print rec1['pick_file']['datetime'], rec1['sta']
                    
                    # get data
                    dat1 = rec1[rec1['channels'][0]]
                    dat2 = rec2[rec2['channels'][0]]
                    
                    # set low S/N recs to nan
                    idx = where(dat1['sn_ratio'] < sn_thresh)[0]
                    dat1['swave_spec'][idx] = nan
                    
                    idx = where(dat2['sn_ratio'] < sn_thresh)[0]
                    dat2['swave_spec'][idx] = nan
                    
                    # get spectral ratio
                    spec_rat = dat1['swave_spec'] / dat2['swave_spec']
                    
                    # add to data array
                    if rec1['pick_file']['rhyp'] >= 00. and rec1['pick_file']['rhyp'] < 1500.:
                        if len(stacked_ratios) == 0:
                            stacked_ratios = spec_rat
                        else:
                            stacked_ratios = vstack((stacked_ratios, spec_rat))
                
freqs = records[0][records[0]['channels'][0]]['freqs']

# get mean ratios
median_rat = []
for i in range(0, len(freqs)):
    median_rat.append(exp(nanmedian(log(stacked_ratios[:,i]))))

            
###############################################################################
#plot data
###############################################################################
fig = plt.figure(1, figsize=(12,12))

for rat in stacked_ratios:
    plt.loglog(freqs, rat, '-', lw=0.75, c='dodgerblue')

plt.loglog(freqs, median_rat, '-', lw=2., c='k')
plt.grid(which='both')
plt.show()

###############################################################################
# plot wood_anderson_attenuation
###############################################################################

"""plt.figure(2, figsize=(10,10))
plt.plot(distdict, ampdict, '+', color='0.6')
plt.xlim([log10(3), 3])
#plt.ylim([-4.5, -0.25])

#bin data
bins = arange(0., 3.1, 0.1)
#bins = log10(hstack((arange(5., 101., 5.), arange(110., 201., 10.), \
#                     arange(220., 301., 20.), arange(350., 601., 50.))))
#bins = log10(arange(10., 800, 20))
medbin, stdbin, meanx, outbins, nperbin = get_binned_stats(bins, distdict, ampdict)
idx = where(~isnan(medbin))[0]
plt.errorbar(meanx, medbin, yerr=stdbin, color='r', fmt='s', mec='r')

# try separating Gippsland vs non_Gippsland
#gidx = where(evla < -37.8)[0]
#medbin, stdbin, meanx, outbins, nperbin = get_binned_stats(bins, log10(rhyp[gidx]), lognormamp[gidx])
#idx = where(~isnan(medbin))[0]
#plt.errorbar(meanx, medbin, yerr=stdbin, color='b', fmt='^', mec='b')

plt.xlabel('log Hypocentral Distance', fontsize=16)
plt.ylabel('Normalised log Wood-Anderson (mm)', fontsize=16) 
plt.title('Geometrical Attenuation Regression')
plt.grid(which='major')

###############################################################################
# plot regression
###############################################################################
# fit all data 
data = odrpack.RealData(meanx, medbin) # regress binned data
#data = odrpack.RealData(log10(rhyp), lognormamp) # regress all data

R1 = 70.
R2 = 900.
# get dists LE 90
didx = where(10**meanx <= R1)[0]
tri1 = linregress(meanx[didx], medbin[didx])
#tri1 = linregress(log10(rhyp[didx]), lognormamp[didx]) # regress all data

def fit_b2(c, x):
    # c0 defined above
    
    from numpy import zeros_like, where
    
    ans = zeros_like(x)
    
    idx1 = where(x <= log10(R1))[0]
    idx2 = where((x > log10(R1)) & (x <= log10(R2)))[0]
        
    ans[idx1] = c0 + b1 * x[idx1]
    
    ans[idx2] = c0 + b1 * log10(R1) \
                + c[0] * log10(10**x[idx2] / R1)
 
    return ans

b1 = tri1[0]
#2 = out.beta[0] # 90-150 km
c0 = tri1[1]

didx = where(10**meanx <= R2)[0]

data = odrpack.RealData(meanx[didx], medbin[didx]) # regress binned data
#data = odrpack.RealData(log10(rhyp[didx]), lognormamp[didx]) # regress all data

tri2 = odrpack.Model(fit_b2)
odr = odrpack.ODR(data, tri2, beta0=[0.0])
odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq
out = odr.run()

b2 = out.beta[0] # 90-150 km
#b3 = out.beta[1] # GT 150 km
    
#tri3 = odrpack.Model(fit_b3)
#odr = odrpack.ODR(data, tri3, beta0=[-1.7])
#odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq
#out = odr.run()
#out = odr.run()
#
#b3 = out.beta[0] # 90-150 km
#b2 = out.beta[0] # 90-150 km
##b3 = out.beta[1] # GT 150 km


# plot model
x1 = log10([5, R1])
y1 = b1*x1 + c0
plt.plot(x1, y1, '-', color='limegreen', lw=2.0)

x2 = log10([R1, R2])
y2 = c0 + b1 * log10(R1) \
     + b2 * log10(10**x2 / R1)
plt.plot(x2, y2, '-', color='limegreen', lw=2.0)


plt.savefig('norm_logamps.png', bbox_inches='tight')

plt.show()"""

