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

fidx = 75
###############################################################################
# load pickle 
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

###############################################################################
# set datasets
###############################################################################

events = unique(dictlist2array(recs, 'ev'))
mags = dictlist2array(recs, 'mag')
stations = unique(dictlist2array(recs, 'sta'))
networks = unique(dictlist2array(recs, 'net'))

#fidx = 35

stalist = stationlist2dict()
stalist_start = dictlist2array(stalist, 'start')
stalist_sta = dictlist2array(stalist, 'sta')

#fidx = 75
chan = recs[0]['channels'][0]
freq = recs[0][chan]['freqs'][fidx]
print("Reg Freq = " +str('%0.3f' % freq))

if not freq == c['freq']:
   print('\n!!!! Frequency of coefficients inconsistent with data index !!!!\n')
   crash

###############################################################################
# parse coefs and get model prediction
###############################################################################

rhyps = []
yres = []    
for i, rec in enumerate(recs):
    try:
        channel = rec['channels'][0]
            
        if rec[channel]['sn_ratio'][fidx] >= 5.:
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

###############################################################################
# get stns res
###############################################################################
dateRng = [UTCDateTime(1989,12,1).datetime, UTCDateTime(2023,1,1).datetime]
fig = plt.figure(1, figsize=(19,11))
i = 1
ii = 1
residuals = []
for ev in events:
    
    staRes = []
    stas = []
    rhyps = []
    nets = []
    for rec in recs:
        if rec['ev'] == ev:
            staRes.append(rec['yres'])
            stas.append(rec['sta'])
            rhyps.append(rec['rhyp'])
            mag = rec['mag']
            nets.append(rec['net'])
            place = rec['place']
            if len(place) > 50:
                place = place[0:50]
    
    # fill residual array
    staRes = array(staRes)
    nets = array(nets)
    rdict = {'ev':ev, 'stas':stas, 'stares': staRes, 'rhyps':rhyps, 'mag':mag, 'nets':nets}
    	
    # add AU mean
    nidx = where((nets=='AU') | (nets=='IU') | (nets=='II'))[0]
    if len(nidx) > 0:
        meanAU = nanmean(staRes[nidx])
    else:
        meanAU = 0.
    aurelres = staRes - meanAU
    rdict['aurelres'] = aurelres
    residuals.append(rdict)
    
    if len(rhyps) >= 0:
        ax = plt.subplot(4,2,i)
        plt.plot(rhyps, staRes, 'r+', lw=1)
        plt.plot([0, 2300], [0, 0], 'k--', lw=0.75)
        plt.ylim([-2., 2.])
        plt.xlim([0, 2200])
        plt.title(' - '.join((ev, 'M'+str('%0.2f' % mag), place)), fontsize=10)
        
        # now plot sta name
        for rhyp, res, sta, net in zip(rhyps, staRes, stas, nets):
            if not isnan(res):
                stanet = '.'.join((sta, net))
                plt.text(rhyp+10, res+0.05, stanet, ha='left', va='top', fontsize=7)
        
        if i == 1 or i == 3 or i == 5 or i == 7:
           plt.ylabel('log Obs - log Pred')
        
        i += 1
    
    if i > 8:
        plt.savefig('ev_res/ev_res_'+str(ii)+'.png', fmt='png', bbox_inches='tight')
        i = 1
        ii += 1
        fig = plt.figure(ii, figsize=(19,11))
        
plt.savefig('ev_res/ev_res_'+str(ii)+'.png', fmt='png', bbox_inches='tight')
plt.show()

###############################################################################
# write file
###############################################################################

txt = 'ORIGIN_TIME,STA,NET,MAG,RHYP,STA_RES,AU_REL_RES\n'
for rdict in residuals:
    for i, sta in enumerate(rdict['stas']):
        if not isnan(rdict['stares'][i]):
            txt += ','.join((str(rdict['ev']), sta, rdict['nets'][i], str('%0.2f' % rdict['mag']), str(rdict['rhyps'][i]), \
                             str(rdict['stares'][i]), str(rdict['aurelres'][i]))) + '\n'
                         
f = open('station_residuals.csv', 'w')
f.write(txt)
f.close()