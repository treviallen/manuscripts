import pickle
from numpy import unique, array, arange, log, log10, logspace, exp, mean, nanmean, ndarray, \
                  nanmedian, hstack, pi, nan, isnan, interp, polyfit, where, zeros_like, polyfit, sqrt, \
                  floor, zeros
from misc_tools import get_binned_stats, dictlist2array, get_mpl2_colourlist, savitzky_golay, dictlist2array, get_log_xy_locs
from mag_tools import nsha18_mb2mw, nsha18_ml2mw
from get_mag_dist_terms_swan import get_distance_term
from scipy.odr import Data, Model, ODR, models
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import linregress
import scipy.odr.odrpack as odrpack
from obspy import UTCDateTime
from sys import argv
mpl.style.use('classic')
plt.rcParams['pdf.fonttype'] = 42
import warnings
warnings.filterwarnings("ignore")
            
###############################################################################
# grunt defs
###############################################################################
#print('\n ASSIGN MW AND RE-RUN \n')
def normalise_data(p, recs, sn_ratio, events):
    
    #print('!!!!!!! UNCOMMENT FILTER BY SAMPLE RATE !!!!!!!')
    
    log_norm_amps = []
    stas = []

    chan = recs[0]['channels'][-1]
    
    print("Reg Freq = " +str('%0.3f' % freq))
    
    i = 1
    for e in events:
        cnt = 0
        
        mrhyps = []
        mamps  = []
        mmags = []
        mstas = []
    
        # get all records for each sta
        for rec in recs:
            if len(rec['channels']) > 0 and rec['ev'] == e:
                if rec['net'] in keep_nets:
                    if not rec['sta'] in ignore_stas:
                        channel = rec['channels'][0]
                        
                        # filter by instrument type
                        addData = True
                        if rec[channel]['freqs'][fidx[p]] < 0.4:
                            if  channel.startswith('SH') or channel.startswith('EH'):
                                addData = False
                        
                        # filer by sample-rate
                        if rec[channel]['freqs'][fidx[p]] > 0.45 * rec[channel]['sample_rate']:
                            addData = False
                        
                        # ignore dodgy CMSA data
                        if rec['sta'] == 'CMSA' and rec[channel]['freqs'][fidx[p]] < 0.5:
                            addData = False
                            
                        if rec[channel]['sn_ratio'][fidx[p]] >= sn_ratio and addData == True:
                            mrhyps.append(rec['rhyp'])
                            mmags.append(rec['mag'])
                            mstas.append(rec['sta'])
                            mamps.append(rec[channel]['swave_spec'][fidx[p]])
                        
        
        mrhyps = array(mrhyps)
        mamps = array(mamps)
        mmags = array(mmags)
        mstas = array(mstas)
        
        if len(mrhyps) > 0:
            
            i += 1
            
            # get binned data
            logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log10(mrhyps), log10(mamps))
            
            # normalise data @ 500 km
            nidx = where((binstrp > 2.69) & (binstrp < 2.71))[0] # 500 km
            nidx = where((binstrp >= 2.25) & (binstrp < 2.35))[0] # 200 km
            #nidx = where((binstrp > 1.99) & (binstrp < 2.01))[0] # 100 km
            
            if len(nidx) > 0:
                #print (m, nidx)
                namps = log10(mamps) - logmedamp[nidx]
                
                if len(log_norm_amps) == 0:
                    log_norm_amps = namps
                    logamps = log10(mamps)
                    norm_rhyps = mrhyps
                    mags = mmags
                    stas = mstas
                else:
                    log_norm_amps = hstack((log_norm_amps, namps))
                    logamps = hstack((logamps, log10(mamps)))
                    norm_rhyps = hstack((norm_rhyps, mrhyps))
                    mags = hstack((mags, mmags))
                    stas = hstack((stas, mstas))
    
    return log_norm_amps, norm_rhyps, mags, logamps, stas
    
###############################################################################
# load pickle 
###############################################################################

recs = pickle.load(open('fft_data.pkl', 'rb' ))

# remove bad recs
keep_nets = set(['AU', 'IU', 'S1', 'G', 'MEL', 'ME', '20', 'AD', 'SR', 'UM', 'AB', 'VI', 'GM' \
                 '1P', '1P', '2P', '6F', '7K', '7G', 'G', '7B', '4N', '7D', '', 'OZ', 'OA', 'WG','XX'])

# get stas to ignore
ignore_stas = open('sta_ignore.test').readlines()
ignore_stas = set([x.strip() for x in ignore_stas])

###############################################################################
# parse preliminary Mw and assign as mag
###############################################################################

lines = open('brune_stats.csv').readlines()[1:]

brune_ev = []
brune_mw = []
brune_flag = [] # if trust Mw

for line in lines:
    dat = line.strip().split(',')
    brune_ev.append(dat[0])
    brune_mw.append(dat[6])
    brune_flag.append(0) # for now!

####################################################################################
# start main
####################################################################################

events = unique(dictlist2array(recs, 'ev'))
mags = dictlist2array(recs, 'mag')
magTypes = dictlist2array(recs, 'magType')
rhyp = dictlist2array(recs, 'rhyp')
datetimes = dictlist2array(recs, 'ev')
omag = mags

# reset mag to brune mw
for i, event in enumerate(events):
    for j, bev in enumerate(brune_ev):
        if brune_flag == 1:
            mags[i] = brune_mw[j]
            #magTypes[i] = 

stations = unique(dictlist2array(recs, 'sta'))

# convert mags to MW
for i, mt in enumerate(magTypes):
    if mt.lower().startswith('mb'):
        mags[i] = nsha18_mb2mw(mags[i])
    elif mt.lower().startswith('ml'):
        # additional fix for use of W-A 2800 magnification pre-Antelope
        if UTCDateTime(datetimes[i]) < UTCDateTime(2008, 1, 1):
            mags[i] -= 0.07
        
        # now fix ML
        mags[i] = nsha18_ml2mw(mags[i])

mrng = arange(3.5, 6.9, 0.1)
minDist = 10**0.5
maxDist = 2000
#maxDist = 1300
minRegDist = 100
maxRegDist = 1000

bins = arange(log10(minDist), log10(maxDist), 0.1)
    
# loop through freqs
rec = recs[0]
chan = rec['channels'][-1]
freqs = rec[chan]['freqs']

sn_ratio = 4
    
# load atten coeffs
coeffs = pickle.load(open('atten_coeffs.pkl', 'rb' ))

fidx = [76, 117]
fidx = [76, 99]

xplt = arange(1,2200,1)
for p, freq in enumerate(freqs[fidx]):
    fig = plt.figure(1, figsize=(12,4))
    ax = plt.subplot(1,2,p+1)
    
    
    # set freq-depeendent sn_ratio
    sn_ratio = 4
    
    # get normalised amplitudes
    log_norm_amps, norm_rhyps, mags, logamps, stas = normalise_data(p, recs, sn_ratio, events)
    plt.loglog(norm_rhyps, 10**log_norm_amps, '+', c='0.6', lw=0.5, ms=6)
    
    log_norm_rhyps = log10(norm_rhyps)
    logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log_norm_rhyps, log_norm_amps)
    
    plt.loglog(10**medx, 10**logmedamp, 'rs', ms=7)
    
    # calc distance specific terms
    c = coeffs[fidx[p]]	
    
    # get distance term
    distterms = []
    for x in xplt:
        distterm = get_distance_term(x, c)
        distterms.append(distterm)
    
    distterms = array(distterms)
    
    yplt = 10**(c['nc1s'] + distterms)
    plt.loglog(xplt, yplt, 'k-', lw=2.0)    
        
    plt.grid(which='major')
    plt.xlim([2, 500])
    plt.ylim([0.001, 1000])
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    xloc = get_log_xy_locs(xlims, 0.04)
    yloc = get_log_xy_locs(ylims, 0.04)
    props = dict(boxstyle='round', facecolor='w', alpha=1)
    plt.text(xloc, yloc, 'f = '+str('%0.1f' % freq)+' Hz', va='bottom', ha ='left', fontsize=12, bbox=props)
    if p == 0:
        plt.ylabel('Normailised Fourier Amplitude')
    plt.xlabel('Hypocentral Distance (km)')
    
plt.savefig('norm_geom_spread_paper.png', fmt='png', dpi=150, bbox_inches='tight')
plt.show()
