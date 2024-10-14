import pickle
from numpy import unique, array, arange, log, log10, exp, mean, nanmean, ndarray, sqrt, \
                  nanmedian, hstack, pi, nan, isnan, interp, where, zeros_like, polyfit, logspace
from misc_tools import get_binned_stats, dictlist2array, get_mpl2_colourlist
from mag_tools import nsha18_mb2mw, nsha18_ml2mw
from get_mag_dist_terms import get_distance_term, get_magnitude_term, get_kappa_term
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.dates as mdates
from scipy.stats import linregress
import scipy.odr.odrpack as odrpack
from obspy import UTCDateTime
from response import stationlist2dict
mpl.style.use('classic')
plt.rcParams['pdf.fonttype'] = 42
import warnings
warnings.filterwarnings("ignore")

fidx = 30 # 0.13
#fidx = 50 # 0.3Hz
fidx = 75 # 1Hz
fidx = 90 # 2Hz
#fidx = 115 # 6.3 Hz
#fidx = 100 # 3.2 Hz
#fidx = 110 # 5
###############################################################################
# load datasets 
###############################################################################
recs = pickle.load(open('fft_data.pkl', 'rb' ))


# set most recent Brune mags
"""
lines = open('brune_stats.csv').readlines()[1:]
brunedat = []
for line in lines:
    dat = line.strip().split(',')
    tmp = {'datetime':UTCDateTime(dat[0]), 'mw':float(dat[8]), 'qual':int(float(dat[-1]))}
    
    brunedat.append(tmp)
    
def get_brune_deets(rec_ev):
    bruneStats = {'qual':0}
    for evnum, ev in enumerate(brunedat): 
        if rec_ev == ev['datetime']:
            #print(rec_ev)
            bruneStats = ev
               
    return bruneStats

# convert mags to MW
for i, rec in enumerate(recs):
    bruneMag = False
    bruneStats = get_brune_deets(rec['evdt'])
    
    if bruneStats['qual'] > 0:
        recs[i]['mag'] = bruneStats['mw']
        bruneMag = True
    
    if bruneMag == False:
        if rec['magType'].lower().startswith('mb'):
            recs[i]['mag'] = nsha18_mb2mw(rec['mag'])
        elif rec['magType'].lower().startswith('ml'):
            # additional fix for use of W-A 2800 magnification pre-Antelope
            if UTCDateTime(rec['ev']) < UTCDateTime(2008, 1, 1):
                recs[i]['mag'] -= 0.07
            
            # now fix ML
        recs[i]['mag'] = nsha18_ml2mw(rec['mag'])
"""
# dump with updated mags
'''
pklfile = open('fft_data.pkl', 'wb')      
pickle.dump(recs, pklfile, protocol=-1)
pklfile.close()                           
'''

# load coeffs
coeffs = pickle.load(open('atten_coeffs.pkl', 'rb' ))
#coeffs = pickle.load(open('atten_coeffs_1f.pkl', 'rb' ))

c = coeffs[fidx]
#c = coeffs[0]
print("Coeffs Freq = " +str('%0.3f' % c['freq']))

# load station sets
lines = open('station_sets.csv').readlines()
sta_sets = []
for line in lines:
    sta_sets.append(set(line.strip().split(',')))

###############################################################################
# set datasets
###############################################################################

events = unique(dictlist2array(recs, 'ev'))
mags = dictlist2array(recs, 'mag')
stations = unique(dictlist2array(recs, 'sta'))

chan = recs[0]['channels'][0]
freq = recs[0][chan]['freqs'][fidx]
print("Reg Freq = " +str('%0.3f' % freq))

if not freq == c['freq']:
   print('\n!!!! Frequency of coefficients inconsistent with data index !!!!\n')
   crash

stalist = stationlist2dict()
stalist_start = dictlist2array(stalist, 'start')
stalist_sta = dictlist2array(stalist, 'sta')


mdist_lookup_mags = arange(3.5,7.1,0.5)
mdist_lookup_dists = array([550, 1200, 1700, 2000, 2200, 2200, 2200, 2200])


###############################################################################
# parse coefs and get model prediction
###############################################################################
rhyps = []
yres = []    
for i, rec in enumerate(recs):
    try:
        channel = rec['channels'][0]
        
        idx = where(rec['mag'] >= mdist_lookup_mags)[0]
        if len(idx) == 0:
            mag_dist = mdist_lookup_dists[0]
        else:
            mag_dist = mdist_lookup_dists[idx[-1]]
            
        if rec[channel]['sn_ratio'][fidx] >= 4. and rec['mag'] <= mag_dist:
            rhyps.append(rec['rhyp'])
            
            # get mag term
            magterm = get_magnitude_term(rec['mag'], c)
            
            # get dist term
            distterm = get_distance_term(rec['rhyp'], c)
             
            #	get distance independent kappa
            kapterm = get_kappa_term(rec['sta'], c['freq'])
            
            # get total correction
            ypred = magterm + distterm + kapterm
            
            yobs = log10(rec[channel]['swave_spec'][fidx])
            yres.append(yobs - ypred)
            recs[i]['yres'] = yobs - ypred
            recs[i]['yobs'] = yobs
            recs[i]['ypred'] = ypred
            recs[i]['magterm'] = magterm
            recs[i]['distterm'] = distterm
            recs[i]['kapterm'] = kapterm
            
        else:
            yres.append(nan)
            rhyps.append(rec['rhyp'])
            recs[i]['yres'] = nan
    
    except:
        print('No data')
        recs[i]['yres'] = nan
    
fig = plt.figure(1, figsize=(18,5))

plt.semilogx(rhyps, yres, '+', c='0.5')
plt.semilogx([5, 2200], [0, 0], 'k--')
plt.xlim([5, 2200])
plt.ylim([-3, 3])

# get binned data
logbins = arange(0.1, 3.5, 0.1)
logresbin, stdbin, medx, binstrp, nperbin = get_binned_stats(logbins, log10(rhyps), yres)
plt.semilogx(10**medx, logresbin, 'rs', ms=7)

plt.show()
crash
###############################################################################
# get stns res
###############################################################################
dateRng = [UTCDateTime(1989,12,1).datetime, UTCDateTime(2024,8,30).datetime]
fig = plt.figure(1, figsize=(19,11))
i = 1
ii = 1
for sta in stations:
    # check station sets
    sta_set = set([sta])
    
    for ss in sta_sets:
        if sta in ss:
            sta_set = ss
               
    staRes = []
    staDate = []
    for rec in recs:
        if rec['sta'] in sta_set:
            staRes.append(rec['yres'])
            staDate.append(UTCDateTime(rec['ev']).datetime)
            
    if len(staRes) >= 3:
        ax = plt.subplot(4,2,i)
        plt.plot(staDate, staRes, 'r+', lw=1)
        plt.plot(dateRng, [0, 0], 'k--', lw=0.75)
        plt.ylim([-2, 2])
        plt.xlim(dateRng)
        plt.title(sta)
        
        # now plot instrument changes
        start_times = []
        for sl_sta, sl_start in zip(stalist_sta, stalist_start):
            if sl_sta == sta:
                start_times.append(sl_start)
                
        start_times = unique(array(start_times))
        for st in start_times:
            plt.plot([st, st], [-2, 2], 'g--', lw=0.75)
        
        if i == 1 or i == 3 or i == 5 or i == 7:
           plt.ylabel('log Obs - log Pred')
        #plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        #plt.gcf().autofmt_xdate()
        
        i += 1
    
    if i > 8:
        plt.savefig('sta_res/sta_res_'+str(ii)+'.png', fmt='png', bbox_inches='tight')
        i = 1
        ii += 1
        fig = plt.figure(ii, figsize=(19,11))
        
plt.savefig('sta_res/sta_res_'+str(ii)+'.png', fmt='png', bbox_inches='tight')
plt.show()
    