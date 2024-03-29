import pickle
from numpy import unique, array, arange, log, log10, logspace, exp, mean, nanmean, ndarray, \
                  nanmedian, hstack, pi, nan, isnan, interp, polyfit, where, zeros_like, polyfit, sqrt, \
                  floor, zeros
from misc_tools import get_binned_stats, dictlist2array, get_mpl2_colourlist, savitzky_golay, dictlist2array
from mag_tools import nsha18_mb2mw, nsha18_ml2mw
#from js_codes import get_average_Q_list, extract_Q_freq
#from mapping_tools import distance
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

# = True of plotting/testing; = False for full regression
pltTrue = argv[1]

# get near-source fit
nref = 5
def fit_near_source_saturation(c, x):
    from numpy import sqrt
    
    D = sqrt((10**x)**2 + nref**2) # were c1 ~= 5-10?
    #D = sqrt((10**x)**2 + nref**2)
    
    ans = c[0] * log10(D) + c[1]
    
    return ans

            
###############################################################################
# grunt defs
###############################################################################
print('\n ASSIGN MW AND RE-RUN \n')
def normalise_data(p, recs, sn_ratio):
    
    #print('!!!!!!! UNCOMMENT FILTER BY SAMPLE RATE !!!!!!!')
    
    log_norm_amps = []
    stas = []

    chan = recs[0]['channels'][-1]
    
    print("Reg Freq = " +str('%0.3f' % freq))
    
    i = 1
    for m in mrng:
        cnt = 0
        
        mrhyps = []
        mamps  = []
        mmags = []
        mstas = []
    
        # get all records for each sta
        for rec in recs:
            if len(rec['channels']) > 0 and rec['mag'] >= m-0.05 and rec['mag'] < m+0.05:
                if rec['net'] in keep_nets:
                    if not rec['sta'] in ignore_stas:
                        channel = rec['channels'][0]
                        
                        # filter by instrument type
                        addData = True
                        if rec[channel]['freqs'][fidx[p]] < 0.8:
                            if  channel.startswith('SH') or channel.startswith('EH'):
                                addData = False
                        
                        
                        # filer by sample-rate
                        if rec[channel]['freqs'][fidx[p]] > 0.45 * rec[channel]['sample_rate']:
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
            nidx = where((binstrp > 2.49) & (binstrp < 2.51))[0] # 250 km
            nidx = where((binstrp > 1.99) & (binstrp < 2.01))[0] # 100 km
            
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
    
def fit_near_source_atten(plt, norm_rhyps, log_norm_amps, r1, stas):
    if pltTrue == True:
        plt.loglog(norm_rhyps, 10**log_norm_amps, '+', c='0.6', lw=0.5, ms=6)
        plt.xlim([4, 2250])
        
        # plot dodgy stas
        for sta, log_norm_amp, norm_rhyp in zip(stas, log_norm_amps, norm_rhyps):
            if log_norm_amp > 0.6 and norm_rhyp < 160.:
                plt.text(norm_rhyp, 10**log_norm_amp, sta, fontsize=6)
            elif log_norm_amp < -1 and norm_rhyp < 100.:
                plt.text(norm_rhyp, 10**log_norm_amp, sta, fontsize=6)

    # get binned data
    log_norm_rhyps = log10(norm_rhyps)
    logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log_norm_rhyps, log_norm_amps)
    
    if pltTrue == True:
        #idx = where(logmedamp < 1.88)[0]
        plt.loglog(10**medx, 10**logmedamp, 'rs', ms=7)
    
    
    # fit all data - keep me!
    didx = where((medx >= log10(minr)) & (medx < log10(r1)) & (nperbin > 2) & (isnan(medx) == False))[0] 
    data = odrpack.RealData(medx[didx], logmedamp[didx])
    
    didx = where((log_norm_rhyps >= log10(minr)) & (log_norm_rhyps <= log10(r1)))[0]
    data = odrpack.RealData(log10(norm_rhyps[didx]), log_norm_amps[didx])
    
    afit = odrpack.Model(fit_near_source_saturation)
    #odr = odrpack.ODR(data, afit, beta0=[-1., 2, 1.])
    odr = odrpack.ODR(data, afit, beta0=[-1.5, 2])
    
    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
    out = odr.run()
    nc = out.beta
    
    return nc
    
def refit_near_source_atten(plt, norm_rhyps, log_norm_amps, r1, n0):
    
    # get binned data
    log_norm_rhyps = log10(norm_rhyps)
    logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log_norm_rhyps, log_norm_amps)
    
    
    def refit_near_source_saturation(c, x):
        from numpy import sqrt
        
        D = sqrt((10**x)**2 + nref**2) # were c1 ~= 5-10?
        
        ans = n0 * log10(D) + c[0]
        
        return ans
        
    # fit all data - keep me!
    didx = where((medx >= log10(minr)) & (medx < log10(r1)) & (nperbin > 2) & (isnan(medx) == False))[0] 
    data = odrpack.RealData(medx[didx], logmedamp[didx])
    
    #didx = where((log_norm_rhyps >= log10(minr)) & (log_norm_rhyps <= log10(r1)))[0]
    #data = odrpack.RealData(log10(norm_rhyps[didx]), log_norm_amps[didx])
    
    afit = odrpack.Model(refit_near_source_saturation)
    #odr = odrpack.ODR(data, afit, beta0=[-1., 2, 1.])
    odr = odrpack.ODR(data, afit, beta0=[2])
    
    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
    out = odr.run()
    nc_intercept = out.beta
    
    return nc_intercept

def fit_mid_field_atten(plt, norm_rhyps, log_norm_amps, n0, r1, r2, freq):
    #print(r2)
    # fit all data
    log_norm_rhyps = log10(norm_rhyps)
    
    logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log_norm_rhyps, log_norm_amps)
    
    # fit all data - keep me!
    didx = where((medx >= log10(r1)) & (medx < log10(r2)) & (isnan(medx) == False) & (logmedamp < 2))[0]
    data = odrpack.RealData(medx[didx], logmedamp[didx])
    
    #print(medx[didx])
    #print(logmedamp[didx])
    
    # if few datapoints, use full dataset
    #if len(didx) <= 2:
    #didx = where((log_norm_rhyps >= log10(r1)) & (log_norm_rhyps < log10(r2)))[0] # & (nperbin > 2))[0] #log10(r1))[0]
    #data = odrpack.RealData(log_norm_rhyps[didx], log_norm_amps[didx])

    # fit all as free params
    afit = odrpack.Model(fit_mid_field_first)
    odr = odrpack.ODR(data, afit, beta0=[-0.1, -0.001, 0.2])
    
    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
    out = odr.run()
    mcf = out.beta
    mc = zeros(3)
    #print('mcf',mcf)
    #nc = zeros(2)        
    
    mc[0] = mcf[0]
    mc[1] = mcf[1]
    mc[2] = mcf[2]
    
    # do some checks here
    xtest=arange(r1,r2,1)
    ytest = mc[0] * log10(xtest / r1) + mc[1] * (xtest - r1)
    
    mididx = int(floor(len(xtest)/2))
    
    
    # if conditions met, refit with log-linear 
    '''
    if ytest[mididx] > ytest[0] or ytest[-1] > ytest[0] or ytest[1] > ytest[0] or ytest[mididx] < ytest[-1] or mc[1] > 0.0:
       
       #print('refitting mid')
       afit = odrpack.Model(refit_mid_field_first)
       odr = odrpack.ODR(data, afit, beta0=[-0.2, 0.2])
       
       odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
       out = odr.run()
       remc = out.beta
       
       mc = zeros(3)
       mc[0] = remc[0]
       mc[1] = 0.
       mc[2] = remc[1]
    '''   
    # recalculate n1
    D1 = sqrt(r1**2 + nref**2)
    
    n1_hold = mc[0] * log10(r1 / r1) + mc[1] * (r1 - r1) + mc[2] - n0 * log10(D1)
    
    # plot
    xrng_mf = arange(log10(r1), log10(r2), 0.02)
    
    yrng_mf = n0 * log10(D1) + n1_hold + mc[0] * log10(10**xrng_mf / r1) + mc[1] * (10**xrng_mf - r1)
        
    if pltTrue == True:
       plt.loglog(10**xrng_mf, 10**yrng_mf, 'm-', lw=2)
    
    return mc
    
def refit_mid_field_atten(plt, norm_rhyps, log_norm_amps, n0, m0s, r1, r2, freq):
    
    log_norm_rhyps = log10(norm_rhyps)
    
    logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log_norm_rhyps, log_norm_amps)
    
    # fit all data - keep me!
    didx = where((medx >= log10(r1)) & (medx < log10(r2)) & (isnan(medx) == False))[0]
    data = odrpack.RealData(medx[didx], logmedamp[didx])
    
    # if few datapoints, use full dataset
    #if len(didx) <= 2:
    #didx = where((log_norm_rhyps >= log10(r1)) & (log_norm_rhyps < log10(r2)))[0] # & (nperbin > 2))[0] #log10(r1))[0]
    #data = odrpack.RealData(log_norm_rhyps[didx], log_norm_amps[didx])

    # fit all as free params
    afit = odrpack.Model(fit_mid_field_first_mc0_fix)
    odr = odrpack.ODR(data, afit, beta0=[-0.001, 0.2])
    
    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
    out = odr.run()
    mcf = out.beta
    mc = zeros(3)
    print('mcf refit ',mcf)
    
    mc[0] = m0s # smoothed m0
    mc[1] = mcf[0]
    mc[2] = mcf[1]
    
    # recalculate n1
    #n1 = mc[0] * log10(r1 / r1) + mc[1] * (r1 - r1) + mc[2] - n0 * log10(Dy)
        
    xtest=arange(r1,r2,1)
    ytest = mc[0] * log10(xtest / r1) + mc[1] * (xtest - r1)
    
    mididx = int(floor(len(xtest)/2))
    
    # if conditions met, refit with log-linear 
    if ytest[mididx] > ytest[0] or ytest[-1] > ytest[0] or ytest[1] > ytest[0] or ytest[mididx] < ytest[-1] or mc[1] > 0.0:
       
       print('refitting mid 2')
       afit = odrpack.Model(refit_mid_field_first_mc0_fix)
       odr = odrpack.ODR(data, afit, beta0=[0.2])
       
       odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
       out = odr.run()
       remc = out.beta
       
       mc[0] = m0s # smoothed m0
       mc[1] = 0. # assume linear - set to zero
       mc[2] = remc[0]
       
    # recalculate n1
    D1 = sqrt(r1**2 + nref**2)
    
    n1 = mc[0] * log10(r1 / r1) + mc[1] * (r1 - r1) + mc[2] - n0 * log10(D1)

    # plot
    xrng_mf = arange(log10(r1), log10(r2), 0.02)
    
    yrng_mf = n0 * log10(D1) + n1 + mc[0] * log10(10**xrng_mf / r1) + mc[1] * (10**xrng_mf - r1)
        
    if pltTrue == True:
       plt.loglog(10**xrng_mf, 10**yrng_mf, 'g-', lw=2)
    
    return mc, n1
    
def fit_far_field_atten(plt, mc, norm_rhyps, log_norm_amps, n0, n1, r1, r2, r3, freq):
    # fit all data
    log_norm_rhyps = log10(norm_rhyps)
    
    logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log_norm_rhyps, log_norm_amps)
    
    # fit all data - keep me!
    didx = where((medx >= log10(r2)) & (medx < log10(r3)) & (nperbin > 2) & (isnan(medx) == False))[0] 
    data = odrpack.RealData(medx[didx], logmedamp[didx])
    
    didx = where((log_norm_rhyps >= log10(r2)) & (log_norm_rhyps < log10(r3)))[0] # & (nperbin > 2))[0] #log10(r1))[0]
    data = odrpack.RealData(log_norm_rhyps[didx], log_norm_amps[didx])
    
    afit = odrpack.Model(fit_far_field)
    odr = odrpack.ODR(data, afit, beta0=[-1., -0.001])
    
    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
    out = odr.run()
    fc = out.beta
    
    # do some checks here
    xtest=arange(r2,r3,1)
    ytest = fc[0] * log10(xtest / r2) + fc[1] * (xtest - r2)
    
    mididx = int(floor(len(xtest)/2))
    
    # allow to fit data!
    '''
    # if conditions met, refit with log-linear 
    if ytest[mididx] > ytest[0] or ytest[-1] > ytest[0] or ytest[1] > ytest[0] or ytest[mididx] < ytest[-1]:
    
        print('refitting far')
        afit = odrpack.Model(refit_far_field)
        odr = odrpack.ODR(data, afit, beta0=[-0.5])
        
        odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
        out = odr.run()
        refc = out.beta
        
        fc[0] = refc[0]
        fc[1] = 0.
    '''
    # plot
    xrng_ff = arange(log10(r2), log10(r3), 0.02)
    D1 = sqrt(r1**2 + nref**2)
    
    yrng_ff = n0 * log10(D1) + n1 \
              + mc[0] * log10(r2 / r1) + mc[1] * (r2 - r1) \
              + fc[0] * log10(10**xrng_ff / r2) + fc[1] * (10**xrng_ff - r2)
    
    if pltTrue == True:
        plt.loglog(10**xrng_ff, 10**yrng_ff, 'b-', lw=2)
    
    return fc

def get_distance_residuals(n0, n1, mc, fc):
     D1 = sqrt((10**x)**2 + nref**2)
     ans = n0 * log10(D1) + n1
     
     # set mid-field
     idx = where((10**x > r1) & (10**x <= r2))[0]
     D1 = sqrt(r1**2 + nref**2)
     ans[idx] = n0 * log10(D1) + n1 + mc[0] * log10(10**x[idx] / r1) + mc[1] * (10**x[idx] - r1)
     
     idx = where(10**x > r2)[0]
     ans[idx] = n0 * log10(D1) + n1 \
                + mc[0] * log10(r2 / r1) + mc[1] * (r2 - r1) \
                + fc[0] * log10(10**x[idx] / r2) + fc[1] * (10**x[idx] - r2)
   
###############################################################################
# load pickle 
###############################################################################

recs = pickle.load(open('fft_data.pkl', 'rb' ))

# remove bad recs
keep_nets = set(['AU', 'IU', 'S1', 'G', 'MEL', 'ME', '20', 'AD', 'SR', 'UM', 'AB', 'VI', 'GM' \
                 '1P', '1P', '2P', '6F', '7K', '7G', 'G', '7B', '4N', '7D', '', 'OZ'])

# get stas to ignore
ignore_stas = open('sta_ignore.txt').readlines()
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
            magTypes[i] = 

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

mrng = arange(3.7, 6.9, 0.1)
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

if pltTrue == 'False':
    fidx=arange(0,len(freqs),1) # for regressing all coeffs
    pltTrue = False

else:
    fig = plt.figure(1, figsize=(18,11))
    fidx=arange(0,140,12)+2 # for testing
    #fidx=arange(0,12,1)+0 # for testing
    pltTrue = True

    
coeffs = []
nc0_array = []

minr = 5.
r1 = 50 # max dist for near source
r2 = 180
r3 = maxDist
    
# based on initial regression analysis, set mc0
'''
mc0_fix = -0.28
mc0_fix = -0.2
'''
for p, freq in enumerate(freqs[fidx]):
    if pltTrue == True:
        fig = plt.figure(1, figsize=(18,11))
        plt.subplot(3,4,p+1)
        plt.ylabel(str(freq), fontsize=7)
    
    # set freq-depeendent sn_ratio
    sn_ratio = 4
    
    # get normalised amplitudes
    log_norm_amps, norm_rhyps, mags, logamps, stas = normalise_data(p, recs, sn_ratio)
    
    ###############################################################################
    # plt near-field GR
    ###############################################################################
    
    nc = fit_near_source_atten(plt, norm_rhyps, log_norm_amps, r1, stas)
    #print(nc)
    nc = -1.2

    coeffs.append({'nc0': nc, 'r1': r1, 'r2': r2, 'nref':nref, 'freq':freq})
    
    nc0_array.append(nc)
    n1 = nan # not defined till later
    
###############################################################################
# smooth nc coeff
###############################################################################
if pltTrue == True:
    sg_window = 3
    sg_poly = 1
else:
    sg_window = 21
    sg_poly = 3
    
smooth_nc0 = savitzky_golay(array(nc0_array[0:]), sg_window, sg_poly) # near slope

for i, c in enumerate(coeffs):
    if i < 0:
        coeffs[i]['nc0s'] = coeffs[i]['nc0']
    else:
        coeffs[i]['nc0s'] = smooth_nc0[i]


###############################################################################
# fit mid-atten first
###############################################################################

# using smoothed n0 values, refit
mc0_array = []
for p, freq in enumerate(freqs[fidx]):
    if pltTrue == True:
        plt.subplot(3,4,p+1)
    
    # set n0 - use smoothed vals
    n0 = coeffs[p]['nc0s']
    #mc0_fix = coeffs[p]['mc0s']
    #print('n0 '+str(n0)) 
    
    # set freq-depeendent sn_ratio
    if freq >= 0.08 and freq <= 0.5:
        sn_ratio = 4
    else:
        sn_ratio = 4
    
    # get normalised amplitudes
    log_norm_amps, norm_rhyps, mags, logamps, stas = normalise_data(p, recs, sn_ratio)
    
    ###############################################################################
    # plt mid-field GR
    ###############################################################################
    
    # get mid-field fit
    #nref2 = 100.
    '''
    def fit_mid_field(c, x):
        from numpy import sqrt, log10
        
        #D = sqrt((10**x)**2 + nc[2]**2) # were c1 ~= 5-10?
        D1 = sqrt((10**x)**2 + nref**2)
        ans = n0 * log10(D1) + n1
        
        idx = where((10**x > r1) & (10**x <= r2))[0]
        D1 = sqrt(r1**2 + nref**2)
        
        ans[idx] = n0 * log10(D1) + n1 + c[0] * log10(10**x[idx] / r1) + c[1] * (10**x[idx] - r1)
            
        return ans
                
    def fit_mid_field_mc0_fix(c, x):
        from numpy import sqrt, log10
        
        #D = sqrt((10**x)**2 + nc[2]**2) # were c1 ~= 5-10?
        D1 = sqrt((10**x)**2 + nref**2)
        ans = n0 * log10(D1) + n1
        
        idx = where((10**x > r1) & (10**x <= r2))[0]
        D1 = sqrt(r1**2 + nref**2)
        
        ans[idx] = n0 * log10(D1) + n1 + mc0_fix * log10(10**x[idx] / r1) + c[0] * (10**x[idx] - r1)
            
        return ans
    '''    
    def fit_mid_field_first(c, x):
        from numpy import sqrt, log10
        
        ans = c[0] * log10(10**x / r1) + c[1] * (10**x - r1) + c[2]
            
        return ans 

    '''
    def fit_mid_field_first_mc0_fix(c, x):
        from numpy import sqrt, log10
        print('fit_mid_field_first_mc0_fix')
        
        ans = mc0_fix * log10(10**x / r1) + c[0] * (10**x - r1) + c[1]
        return ans
    '''    
    def refit_mid_field_first(c, x):
        from numpy import sqrt, log10
        
        ans = c[0] * log10(10**x / r1) +  c[1]
        return ans
    
    mc1 = fit_mid_field_atten(plt, norm_rhyps, log_norm_amps, n0, r1, r2, freq) # temp mc for smoothing
    
    mc0_array.append(mc1[0])


###############################################################################
# smooth mc0 and re-fit mid atten    
###############################################################################

mc0_array = array(mc0_array)
idx = where(mc0_array < -0.0)[0]
smooth_mc0 = savitzky_golay(mc0_array[idx], sg_window, sg_poly) # mid-slope

#interpolate to missing frequencies
interp_mc0 = interp(freqs[fidx], freqs[fidx][idx], smooth_mc0)

# set to coeffs
for i, c in enumerate(coeffs):
    coeffs[i]['mc0s'] = interp_mc0[i]
    coeffs[i]['mc0'] = mc0_array[i]

# using smoothed m0 values, refit
for p, freq in enumerate(freqs[fidx]):
    if pltTrue == True:
        plt.subplot(3,4,p+1)
        
    # set n0 - use smoothed vals
    n0 = coeffs[p]['nc0s']
    
    # get normalised amplitudes
    log_norm_amps, norm_rhyps, mags, logamps, stas = normalise_data(p, recs, sn_ratio)
            
    ###############################################################################
    # re-fit mid-field with fixed freq-dependent mc0
    ###############################################################################
    mc0_fix = coeffs[p]['mc0s']
    #print('1',mc0_fix)
    def fit_mid_field_first_mc0_fix(c, x):
        from numpy import sqrt, log10
        #print('fit_mid_field_first_mc0_fix')
        #print('mc0_fix',mc0_fix)
        ans = mc0_fix * log10(10**x / r1) + c[0] * (10**x - r1) + c[1]
        #print('2',mc0_fix)
        return ans
        
    def refit_mid_field_first_mc0_fix(c, x):
        from numpy import sqrt, log10
        
        ans = mc0_fix * log10(10**x / r1) + c[0]
        return ans
        
    mc, n1 = refit_mid_field_atten(plt, norm_rhyps, log_norm_amps, n0, mc0_fix, r1, r2, freq)
    
    ###############################################################################
    # plt near-field GR with updated n1
    ###############################################################################

    # plot
    xrng_nf = log10(arange(1, r1+1))
    
    #D = sqrt((10**xrng_nf)**2 + nc[2]**2)
    D1 = sqrt((10**xrng_nf)**2 + nref**2)
    #print(n1)
    yrng_nf = n0 * log10(D1) + n1
    
    if pltTrue == True:
        plt.loglog(10**xrng_nf, 10**yrng_nf, 'k-', lw=2)
        
    ###############################################################################
    # plt far-field GR
    ###############################################################################
    
    # get far-field fit
    #nref2 = 100.
    def fit_far_field(c, x):
        from numpy import sqrt, log10
        
        # set near-field
        D1 = sqrt((10**x)**2 + nref**2)
        ans = n0 * log10(D1) + n1
        
        # set mid-field
        idx = where((10**x > r1) & (10**x <= r2))[0]
        D1 = sqrt(r1**2 + nref**2)
        ans[idx] = n0 * log10(D1) + n1 + mc[0] * log10(10**x[idx] / r1) + mc[1] * (10**x[idx] - r1)
        
        idx = where(10**x > r2)[0]
        ans[idx] = n0 * log10(D1) + n1 \
                   + mc[0] * log10(r2 / r1) + mc[1] * (r2 - r1)  \
                   + c[0] * log10(10**x[idx] / r2) + c[1] * (10**x[idx] - r2)
            
        return ans
        
    '''
    def fit_far_field_mid_first(c, x):
        from numpy import sqrt, log10
        print('fit_far_field_mid_first')
        
        # get y offset at r1
        Dy = sqrt(r1**2 + nref**2)
        yoff = n0 * log10(Dy)
        
        D1 = sqrt((10**x)**2 + nref**2)
        
        ans = n0 * log10(D1) - yoff + mc[2]
        
        # set mid-field
        idx = where((10**x > r1) & (10**x <= r2))[0]
        D1 = sqrt(r1**2 + nref**2)
        ans[idx] = n0 * log10(D1) + mc[0] * log10(10**x[idx] / r1) + mc[1] * (10**x[idx] - r1) + mc[2]
        
        idx = where(10**x > r2)[0]
        ans[idx] = n0 * log10(D1) \
                   + mc[0] * log10(r2 / r1) + mc[1] * (r2 - r1) + mc[2] \
                   + c[0] * log10(10**x[idx] / r2) + c[1] * (10**x[idx] - r2)
        
        return ans
    '''    
    def refit_far_field(c, x):
        from numpy import sqrt, log10
        
        #D = sqrt((10**x)**2 + nc[2]**2) # were c1 ~= 5-10?
        D1 = sqrt((10**x)**2 + nref**2)
        ans = n0 * log10(D1) + n1
        
        # set mid-field
        idx = where((10**x > r1) & (10**x <= r2))[0]
        D1 = sqrt(r1**2 + nref**2)
        ans[idx] = n0 * log10(D1) + n1 + mc[0] * log10(10**x[idx] / r1) + mc[1] * (10**x[idx] - r1)
        
        idx = where(10**x > r2)[0]
        ans[idx] = n0 * log10(D1) + n1 \
                   + mc[0] * log10(r2 / r1) + mc[1] * (r2 - r1) \
                   + c[0] * log10(10**x[idx] / r2)
            
        return ans
        
    fc = fit_far_field_atten(plt, mc, norm_rhyps, log_norm_amps, n0, n1, r1, r2, r3, freq)
    #print(fc)
    
    '''
    
    '''
    
    ###############################################################################
    # fit mag intercept - mag scaling only needs to be rough as is for testing 
    # instrument response
    ###############################################################################
    
    #fig = plt.figure(2, figsize=(10,10))
    def fit_mag_intercept(c, x):
        from numpy import sqrt, log10
        
        D1 = sqrt(x**2 + nref**2) # were c1 ~= 5-10?
        
        ans = n0 * log10(D1) + n0
        
        idx = x >= r1
        D1 = sqrt(r1**2 + nc[2]**2) # were c1 ~= 5-10?
        #D = sqrt(r1**2 + nref**2)
        
        ans[idx] = n0 * log10(D1) + c[0] + fc[0] * log10(x[idx] / r1) + fc[1] * (x[idx] - r1)
        
        return ans
    
    # get residual
    D1 = sqrt(norm_rhyps**2 + nref**2)
    predamps = n0 * log10(D1) + n1
    
    # set mid-field
    idx = where((norm_rhyps > r1) & (norm_rhyps <= r2))[0]
    D1 = sqrt(r1**2 + nref**2)
    predamps[idx] = n0 * log10(D1) + n1 + mc[0] * log10(norm_rhyps[idx] / r1) \
                    + mc[1] * (norm_rhyps[idx] - r1)
    
    # set far-field
    idx = where(norm_rhyps > r2)[0]
    predamps[idx] = n0 * log10(D1) + n1 \
                    + mc[0] * log10(r2 / r1) + mc[1] * (r2 - r1) \
                    + fc[0] * log10(norm_rhyps[idx] / r2) + fc[1] * (norm_rhyps[idx] - r2)
    
    # get residuals
    logres = logamps - predamps
    
    
    if pltTrue == True:
        plt.subplot(3,4,p+1)
        plt.plot(mags, logres, '+', c='0.6', lw=0.5, ms=6)
        plt.ylabel(str(freq), fontsize=8)
    
    
    # bin data
    logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(mrng, mags, logres)
    
    # fit mag linear
    magc = polyfit(medx, logmedamp, 1)
    
    '''
    if pltTrue == True:
        plt.plot(medx, logmedamp, 'rs', ms=6.5)

    # plot mag fit
    if pltTrue == True:
        yrng = magc[0] * mrng + magc[1]
        plt.plot(mrng, yrng, 'k-', lw=2)
    '''
    ###############################################################################
    # make coeffs dict
    ###############################################################################
    
    coeffs[p]['nc1s'] = n1
    #coeffs[p]['mc0'] = mc[0]
    coeffs[p]['mc1'] = mc[1]
    coeffs[p]['fc0'] = fc[0]
    coeffs[p]['fc1'] = fc[1]
    coeffs[p]['magc0'] = magc[0]
    coeffs[p]['magc1'] = magc[1]
    coeffs[p]['freq'] = freq
    
if pltTrue == True:
    plt.savefig('norm_geom_spread.png', fmt='png', dpi=150, bbox_inches='tight')
    plt.show()
        

###############################################################################
# given coeff bugs, set f < 0.1 Hz to 0.1
###############################################################################

coeff_f = dictlist2array(coeffs, 'freq')
idx = where(coeff_f <= 0.1)[0]

fidx = idx[-1]
for i in range(0, len(coeffs)):
    if coeffs[i]['freq'] < coeff_f[fidx]:
        coeffs[i]['nc1s'] = coeffs[fidx]['nc1s']
        coeffs[i]['mc0'] = coeffs[fidx]['mc0']
        coeffs[i]['mc1'] = coeffs[fidx]['mc1']
        coeffs[i]['fc0'] = coeffs[fidx]['fc0']
        coeffs[i]['fc1'] = coeffs[fidx]['fc1']
        

###############################################################################
# write params
###############################################################################
if pltTrue == False:
    pklfile = open('atten_coeffs.pkl', 'wb')
    pickle.dump(coeffs, pklfile, protocol=-1)
    pklfile.close()

'''
txt = 'ln Y = m0 + (m1-4) * M^2 +(m2-4) * M + r1 * log10(sqrt(rhyp^2 + r2^2)) | rhyp <= '+str(r1)+'\n'
txt += 'ln Y = m0 + (m1-4) * M^2 +(m2-4) * M + r1 * log10(sqrt(rref^2 + r2^2)) + r4 * log10(rhyp / rref) + r5 * (rhyp - rref) | rhyp > '+str(r1)+'\n'
txt += 'm0' + '\t' + str(m0) + '\n'
txt += 'm1' + '\t' + str(m1) + '\n'
txt += 'm2' + '\t' + str(m2) + '\n'
txt += 'r1' + '\t' + str(n0) + '\n'
txt += 'r2' + '\t' + str(nc[2]) + '\n'
txt += 'r3' + '\t' + str(n1) + '\n'
txt += 'r4' + '\t' + str(fc[0]) + '\n'
txt += 'r5' + '\t' + str(fc[1]) + '\n'
txt += 'rref' + '\t' + str(r1) + '\n'
txt += 'f' + '\t' + str('%0.3f' % freq) + '\n'
txt += 'fidx' + '\t' + str(fidx)

#print('!!!!!!!! REMEMEBR TO UPDATE OUTPUT !!!!!!!!')
f = open('basic_atten_coeffs_tmp.txt', 'w')
f.write(txt)
f.close()
'''
