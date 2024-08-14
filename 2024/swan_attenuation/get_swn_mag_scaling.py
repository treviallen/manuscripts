import pickle
from numpy import unique, array, arange, log, log10, exp, mean, nanmean, ndarray, std, sqrt, \
                  nanmedian, nanstd, vstack, pi, nan, isnan, interp, where, zeros_like, ones_like, floor, ceil, polyfit, isfinite
from sys import argv
import matplotlib.pyplot as plt
import matplotlib as mpl
from obspy import UTCDateTime
from scipy.odr import Data, Model, ODR, models
import scipy.odr.odrpack as odrpack
from get_mag_dist_terms_swan import get_distance_term
from misc_tools import get_binned_stats, dictlist2array, get_mpl2_colourlist, savitzky_golay, dictlist2array
mpl.style.use('classic')
plt.rcParams['pdf.fonttype'] = 42
import warnings
warnings.filterwarnings("ignore")


####################################################################################
# correct attenuation
####################################################################################
def parse_kappa_data():
    from numpy import array, loadtxt
    
    kapdat = []
    # read parameter file
    lines = open('site_kappa.csv').readlines()[1:]
    for line in lines:
        dat = line.split(',')
        #kap = {'sta':dat[0], 'kappa_0': float(dat[1]), 'kappa_r': float(dat[2])}
        kap = {'sta':dat[0], 'kappa0': float(dat[1]), 'cnt': float(dat[2])}
    
        kapdat.append(kap)
    
    return kapdat
    
def parse_filtering_data():
    
    filtdat = []
    # read parameter file
    lines = open('event_filtering_lookup.csv').readlines()[1:]
    for line in lines:
        print(line)
        dat = line.split(',')
        filt = {'ev':UTCDateTime(dat[0]), 'minf': float(dat[1]), 'maxf': float(dat[2]), 'qual':float(dat[3])}
    
        filtdat.append(filt)
    
    return filtdat
    
def correct_atten(rec, coeffs, kapdat):
    channel = rec['channels'][0]
    raw_fds = rec[channel]['swave_spec']#[:-8]
    sn_ratio = rec[channel]['sn_ratio']#[:-8]
    freqs = rec[channel]['freqs']#[:-8]
    
    sn_thresh = 4. # used 4 in model regression
    max_dist = 500
    
    # loop thru freqs
    distterms = []
    for c in coeffs:
        # get distance term
        
        distterm = get_distance_term(rec['rhyp'], c)
        
        distterms.append(distterm)
    
    #distterms.append(nan) # as no coeffs for last freq
    
    #	get distance independent kappa
    #smedian_kappa = 0.072214 # should read this from file
    kappa = kapdat[-1]['kappa0'] # default kappa

    k_term = log10(exp(-1 * pi * freqs * kappa))
    
    # correct to source
    cor_fds = 10**(log10(raw_fds) - distterms - k_term)
    
    # get data exceeding SN ratio
    idx = where(sn_ratio < sn_thresh)[0]
    cor_fds_nan = cor_fds.copy()
    cor_fds_nan[idx] = nan
    
    # ensure high SN for periods affected by 2ndary microseisms
    #idx = where((sn_ratio < 20) & (freqs >= 0.09) & (freqs <= 0.35))[0]
    #cor_fds_nan[idx] = nan
    
    # if short period only use f > 0.5
    if  channel.startswith('SH') or channel.startswith('EH'):
        idx = where(freqs < 0.5)[0]
        cor_fds_nan[idx] = nan
        
    # ignore dodgy CMSA data
    if rec['sta'] == 'CMSA_FIXED' or rec['sta'] == 'PI207': 
        idx = where(freqs < 0.5)[0]
        cor_fds_nan[idx] = nan
        
    if rec['sta'] == 'NPS': 
        idx = where(freqs < 100)[0]
        cor_fds_nan[idx] = nan
        
    if rec['sta'] == 'AS17' and rec['ev'].year > 2013 and rec['ev'].year < 2018: 
        idx = where(freqs < 100)[0]
        cor_fds_nan[idx] = nan
        
    # for all remove 0.1-0.3 Hz?
            
    return cor_fds, cor_fds_nan, freqs

def parse_filtering_data():
    
    filtdat = []
    # read parameter file
    lines = open('event_filtering_lookup.csv').readlines()[1:]
    for line in lines:
        dat = line.split(',')
        filt = {'ev':UTCDateTime(dat[0]), 'minf': float(dat[1]), 'maxf': float(dat[2]), 'qual':float(dat[3])}
    
        filtdat.append(filt)
    
    return filtdat

####################################################################################
# get mag data
####################################################################################

def parse_swn_brune_data(csvfile):
    
    mwdat = []
    # read parameter file
    lines = open(csvfile).readlines()[1:]
    for line in lines:
        dat = line.split(',')
        filt = {'ev':dat[0], 'minf': float(dat[-3]), 'maxf': float(dat[-2]), 'qual':float(dat[-1]),
        	      'mw': float(dat[6]), 'sd': float(dat[7])}
    
        mwdat.append(filt)
    
    return mwdat
    
def parse_brune_data(csvfile):
    
    mwdat = []
    # read parameter file
    lines = open(csvfile).readlines()[1:]
    for line in lines:
        dat = line.split(',')
        filt = {'ev':dat[0], 'minf': float(dat[-3]), 'maxf': float(dat[-2]), 'qual':float(dat[-1]),
        	      'mw': float(dat[7]), 'sd': float(dat[8])}
    
        mwdat.append(filt)
    
    return mwdat


swn_mwdat = parse_swn_brune_data('brune_stats.csv')
sd_mwdat = parse_brune_data('../../2023/au_stress_drop/brune_stats.csv')

###############################################################################
# set data to use
###############################################################################

# remove bad recs
keep_nets = set(['AU', 'IU', 'S1', 'II', 'G', 'MEL', 'ME', '2O', 'AD', 'SR', 'UM', 'AB', 'VI', 'GM' \
                 '1P', '1P', '2P', '6F', '7K', '7G', 'G', '7B', '4N', '7D', '', 'OZ', 'OA', 'WG', 'XX'])
# get stas to ignore
ignore_stas = open('sta_ignore.test').readlines()
ignore_stas = set([x.strip() for x in ignore_stas])

recs = pickle.load(open('fft_data.pkl', 'rb' ))

# load atten coeffs
coeffs = pickle.load(open('atten_coeffs.pkl', 'rb' ))

# get kappas
kapdat = parse_kappa_data()

# get event filters
filtdat = parse_filtering_data()

####################################################################################
# loop thru data
####################################################################################
print('Use site term too???')
                            
events = unique(dictlist2array(recs, 'ev'))
stations = unique(dictlist2array(recs, 'sta'))
rhyps = dictlist2array(recs, 'rhyp')

log_stack_logfds = []
fds_mw = []

cnt = 0
max_dist = 500.
for e, event in enumerate(events): # [::-1]): #[-2:-1]:
    print(event)
    
    # get upper & lower f for filtering
    minf = 0.8
    maxf = 12
    qual = 1
    for fdat in filtdat:
        if fdat['ev'] == event:
            minf = fdat['minf']
            maxf = fdat['maxf']
            qual = fdat['qual']

    # get swn mw
    m = nan
    for md in sd_mwdat:
        if event == md['ev'][0:16].replace(':','.'):
           if md['qual'] > 0:
               m = md['mw']
               print(m)
               cnt += 1
               
    if isnan(m):
        for md in swn_mwdat:
            if event == md['ev']:
               if md['qual'] > 0:
                   m = md['mw']
    #crash      
    
    '''
    sp += 1
    plt.subplot(2,3,sp)
    log_stack_logfds = []
    handles1 = []
    labels1 = []
    '''
    ceil_log_fds = -99
    floor_log_fds = 99
    
    i = 0
    pltboxes = True
    for rec in recs:
        
        #for rec in recs:
        if rec['net'] in keep_nets:
            if not rec['sta'] in ignore_stas:
                if len(rec['channels']) > 0:
                    if rec['ev'] == event: # will need to cahnge to rec['datetime']
                        print('   '+rec['sta'])
                        
                        # do skip sta checks
                        skip_sta = False
                        if rec['sta'] == 'QIS' and rec['ev'].year <= 1999:
                            skip_sta = True
                        elif rec['sta'] == 'RMQ' and rec['ev'].year <= 1998:
                            skip_sta = True 
                        elif rec['rhyp'] > max_dist:
                            skip_sta = True 
                        
                        if skip_sta == False:
                            cor_fds, cor_fds_nan, freqs = correct_atten(rec, coeffs, kapdat)
                            
                            # set bad data to nan
                            fidx = where((freqs < minf) & (freqs > maxf))[0]
                            cor_fds_nan[fidx] = nan
                            
                            # stack
                            if log_stack_logfds == []:
                                log_stack_logfds = log10(cor_fds_nan)
                            else:
                                log_stack_logfds = vstack((log_stack_logfds, log10(cor_fds_nan)))
                                
                            fds_mw.append(m)
fds_mw = array(fds_mw)
                            
####################################################################################
# loop thru freqs and regress m-scaling
####################################################################################
# = True of plotting/testing; = False for full regression
pltTrue = argv[1]

if pltTrue == 'False':
    fidx=arange(0,len(freqs),1) # for regressing all coeffs
    #fidx = fidx[::4]
    pltTrue = False
    sg_window = 41
    sg_poly = 3

else:
    fig = plt.figure(1, figsize=(18,11))
    fidx=arange(0,90,8)+35 # for testing
    #fidx=arange(0,12,1)+0 # for testing
    pltTrue = True
    sg_window = 3
    sg_poly = 1
    

coeffs = []    
m0 = []
mrng = arange(2.0, 5.5, 0.1)
for p, fi in enumerate(fidx):
    print('Freq = '+str(freqs[fi]))
    print('Freq Idx = '+str(fi))
    # get feq data for plotting
    log_famp = log_stack_logfds[:,fi]

    # fit with quadratic
    idx = where((isnan(log_famp) == False) & (isnan(fds_mw) == False) & (fds_mw >= 2.5))[0]
    print(len(idx))
    qc = polyfit(fds_mw[idx], log_famp[idx], 1)
    #fitted_mscaling = qc[0]*mrng**2 + qc[1]*mrng + qc[2]
    fitted_mscaling = qc[0]*mrng + qc[1]
    m0.append(qc[0])
    coeffs.append({'m0':qc[0]})
    
    # bin data
    logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(mrng, fds_mw[idx], log_famp[idx])
        
    if pltTrue == True:
        plt.subplot(3,4,p+1)
        plt.ylabel(str(freqs[fi]), fontsize=7)
        plt.plot(fds_mw[idx], log_famp[idx], '+', c='0.6', lw=0.5, ms=6)
        plt.plot(medx, logmedamp, 'rs', ms=6.5)
        plt.plot(mrng, fitted_mscaling, 'k-', lw=2)
        plt.xlim([2.0, 5.5])

m0 = array(m0)
smooth_m0 = savitzky_golay(m0, sg_window, sg_poly)
for i in range(0, len(coeffs)):
    coeffs[i]['m0s'] = smooth_m0[i]

# refit intercept
m1 = []
for p, fi in enumerate(fidx):
    m0 = coeffs[p]['m0s']
    
    def refit_mag_scaling(c, x):
        from numpy import sqrt, log10
        
        ans = m0 * x + c[0]
            
        return ans 
    
    log_famp = log_stack_logfds[:,fi]
    	
    idx = where((isnan(log_famp) == False) & (isnan(fds_mw) == False) & (fds_mw >= 2.5))[0]

    data = odrpack.RealData(fds_mw[idx], log_famp[idx])

    intercept = odrpack.Model(refit_mag_scaling)
    odr = odrpack.ODR(data, intercept, beta0=[-15])

    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
    out = odr.run()
    m1.append(out.beta[0])

m1 = array(m1)
smooth_m1 = savitzky_golay(m1, sg_window, sg_poly)

for i in range(0, len(coeffs)):
    coeffs[i]['m1'] = m1[i]
    coeffs[i]['m1s'] = smooth_m1[i]

if pltTrue == False:
    pklfile = open('mag_coeffs.pkl', 'wb')
    pickle.dump(coeffs, pklfile, protocol=-1)
    pklfile.close()
    
else:
    plt.savefig('swn_mag_scaling.png', fmt='png', bbox_inches='tight')
    plt.show()
    