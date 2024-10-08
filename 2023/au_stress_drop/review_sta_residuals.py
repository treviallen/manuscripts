import pickle
from numpy import unique, array, arange, log, log10, exp, mean, nanmean, ndarray, sqrt, \
                  nanmedian, hstack, pi, nan, isnan, interp, where, zeros_like, polyfit
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

fidx = 75 # 1Hz
fidx = 90 # 2Hz
###############################################################################
# load datasets 
###############################################################################

recs = pickle.load(open('fft_data.pkl', 'rb' ))

# convert mags to MW
for i, rec in enumerate(recs):
    if rec['magType'].lower().startswith('mb'):
        recs[i]['mag'] = nsha18_mb2mw(rec['mag'])
    elif rec['magType'].lower().startswith('ml'):
        # additional fix for use of W-A 2800 magnification pre-Antelope
        if UTCDateTime(rec['ev']) < UTCDateTime(2008, 1, 1):
            recs[i]['mag'] -= 0.07
        
        # now fix ML
        recs[i]['mag'] = nsha18_ml2mw(rec['mag'])

# load atten coeffs
coeffs = pickle.load(open('atten_coeffs.pkl', 'rb' ))

c = coeffs[fidx]
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

###############################################################################
# parse coefs and get model prediction
###############################################################################
rhyps = []
yres = []    
for i, rec in enumerate(recs):
    try:
        channel = rec['channels'][0]
            
        if rec[channel]['sn_ratio'][fidx] >= 4.:
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

plt.show()
#crash
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
    