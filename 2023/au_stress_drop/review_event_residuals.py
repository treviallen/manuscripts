import pickle
from numpy import unique, array, arange, log, log10, exp, mean, nanmean, ndarray, sqrt, \
                  nanmedian, hstack, pi, nan, isnan, interp, where, zeros_like, polyfit
from misc_tools import get_binned_stats, dictlist2array, get_mpl2_colourlist
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

###############################################################################
# load pickle 
###############################################################################

recs = pickle.load(open('fft_data.pkl', 'rb' ))

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

###############################################################################
# parse coefs and get model prediction
###############################################################################

lines = open('basic_atten_coeffs_tmp.txt').readlines()
m0 = float(lines[2].strip().split('\t')[1])
m1 = float(lines[3].strip().split('\t')[1])
m2 = float(lines[4].strip().split('\t')[1])
r1 = float(lines[5].strip().split('\t')[1])
r2 = float(lines[6].strip().split('\t')[1])
#r3 = float(lines[7].strip().split('\t')[1])
r4 = float(lines[8].strip().split('\t')[1])
r5 = float(lines[9].strip().split('\t')[1])
rref = float(lines[10].strip().split('\t')[1])
fidx = int(lines[12].strip().split('\t')[1])

chan = recs[0]['channels'][0]
freq = recs[0][chan]['freqs'][fidx]
print("Reg Freq = " +str('%0.3f' % freq))

rhyps = []
yres = []    
for i, rec in enumerate(recs):
    try:
        channel = rec['channels'][0]
            
        if rec[channel]['sn_ratio'][fidx] >= 4.:
            rhyps.append(rec['rhyp'])
            
            if rec['rhyp'] <= rref:
                D = sqrt(rec['rhyp']**2 + r2**2)
                y = m0 + m1 * (rec['mag']-4.)**2. + m2 * (rec['mag']-4.) + r1 * log10(D)
            else:
                D = sqrt(rref**2 + r2**2)
                y = m0 + m1 * (rec['mag']-4.)**2. + m2 * (rec['mag']-4.) + r1 * log10(D) \
                    + r4 * log10(rec['rhyp'] / rref) + r5 * (rec['rhyp'] - rref)
            
            yobs = log10(rec[channel]['swave_spec'][fidx])
            yres.append(yobs - y)
            recs[i]['yres'] = yobs - y
            
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
dateRng = [UTCDateTime(1989,12,1).datetime, UTCDateTime(2022,3,1).datetime]
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
        plt.title(ev)
        
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