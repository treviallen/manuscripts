import pickle
from numpy import unique, array, arange, log, log10, exp, mean, nanmean, ndarray, std, sqrt, \
                  nanmedian, nanstd, vstack, pi, nan, isnan, interp, where, zeros_like, ones_like, floor, ceil
from scipy.stats import linregress, trim_mean
from misc_tools import get_binned_stats, dictlist2array, get_mpl2_colourlist
from mag_tools import nsha18_mb2mw, nsha18_ml2mw
from mag_tools import m02mw
from js_codes import get_average_Q_list, extract_Q_freq
from mapping_tools import distance
from scipy.odr import Data, Model, ODR, models
import scipy.odr.odrpack as odrpack
from ltsfit import lts_linefit
import matplotlib.pyplot as plt
import matplotlib as mpl
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
        dat = line.split(',')
        filt = {'ev':dat[0], 'minf': float(dat[1]), 'maxf': float(dat[2])}
    
        filtdat.append(filt)
    
    return filtdat
    
def correct_atten(rec, coeffs, kapdat):
    channel = rec['channels'][0]
    raw_fds = rec[channel]['swave_spec']
    sn_ratio = rec[channel]['sn_ratio']
    freqs = rec[channel]['freqs']
    
    sn_thresh = 5. # used 4 in model regression
    
    # loop thru freqs
    distterms = []
    for c in coeffs:
        # get geometric spreading
        D1 = sqrt(rec['rhyp']**2 + c['nref']**2)
        if rec['rhyp'] <= c['r1']:
            distterm = c['nc0s'] * log10(rec['rhyp']) #+ c['nc1s']
        
        # set mid-field
        elif rec['rhyp'] > c['r1'] and rec['rhyp'] <= c['r2']:
            D1 = sqrt(c['r1']**2 + c['nref']**2)
            distterm = c['nc0s'] * log10(D1) \
                       + c['mc0'] * log10(rec['rhyp'] / c['r1']) + c['mc1'] * (rec['rhyp'] - c['r1'])
        
        # set far-field
        elif rec['rhyp'] > c['r2']:
            D1 = sqrt(c['r1']**2 + c['nref']**2)
            distterm = c['nc0s'] * log10(D1) \
                       + c['mc0'] * log10(c['r2'] / c['r1']) + c['mc1'] * (c['r2'] - c['r1']) \
                       + c['fc0'] * log10(rec['rhyp'] / c['r2']) + c['fc1'] * (rec['rhyp'] - c['r2'])
        
        distterms.append(distterm)
    
    distterms.append(nan) # as no coeffs for last freq
    
    #	get distance independent kappa
    #smedian_kappa = 0.072214 # should read this from file
    kappa = kapdat[-1]['kappa0'] # default kappa
    
    # get site kappa
    for kap in kapdat:
        if kap['sta'] == rec['sta']:
            kappa = kap['kappa0'] # + kap['kappa_r'] * rec['rhyp'] 
    
    k_term = log10(exp(-1 * pi * freqs * kappa))
        
    cor_fds = 10**(log10(raw_fds) - distterms - k_term)
    
    #print(distterms)
    
    # get data exceeding SN ratio
    idx = where(sn_ratio < sn_thresh)[0]
    cor_fds_nan = cor_fds.copy()
    cor_fds_nan[idx] = nan
    
    # if short period only use f > 0.5
    if  channel.startswith('SH') or channel.startswith('EH'):
        idx = where(freqs < 0.8)[0]
        cor_fds_nan[idx] = nan
            
    return cor_fds, cor_fds_nan, freqs

# now fit Brune model
def fit_brune_model(c, f):
    from numpy import array, log
    '''
    c[0] = omega0
    c[1] = f0
    f    = frequency
    '''
    # set constants
    vs = 3.6 # km/s
    vsm = vs*1000.
    rho = 2800 # kg/m^3
    C = 4. * pi * rho * vsm**3 * 1000. / (0.55 * 2.0 * 0.71)
    #f = array(exp(logf))
    #print(f

    # fit curve
    FittedCurve = log(c[0] / (1 + (f / (c[1]))**2))
    #FittedCurve = C * omega / (1 + (f / c[1])**2)
    
    return FittedCurve
   
####################################################################################
# set def params
####################################################################################

vs = 3.6 # km/s
vsm = vs*1000.
rho = 2800 # kg/m^3
C = 4. * pi * rho * vsm**3 * 1000. / (0.55 * 2.0 * 0.707)

'''
kg/m**3 * (m/s)**3 = (kg . m**3) / (m**3 . s**3) = kg / s**3 * m-s = kg-m  / s**2

1 N-m = 1 (kg-m**2) / s**2
'''

RP=0.55
VHC=0.7071068
FSE=2.0
shear_vel=3.6
density=2.8
R0=1.0
C_USGS = 1/ (RP * VHC * FSE / (4 * pi * density * shear_vel**3 * R0) * 1e-20)
'''
print(C, C_USGS)

nf = 12 # 5 Hz
mindist = minr
maxdist = maxr
minmag = 3.8
minf = 12
maxf = 22
'''

###############################################################################
# load pickle 
###############################################################################

#recs = pickle.load(open('onger_fft_data.pkl', 'rb' ))
recs = pickle.load(open('fft_data.pkl', 'rb' ))

# convert mags to MW
for i, rec in enumerate(recs):
    if rec['magType'].startswith('mb'):
        recs[i]['mag'] = nsha18_mb2mw(rec['mag'])
    elif rec['magType'].startswith('ml'):
        # additional fix for use of W-A 2800 magnification pre-Antelope
        if UTCDateTime(datetimes[i]) < UTCDateTime(2008, 1, 1):
            recs[i]['mag'] -= 0.07
        
        # now fix ML
        recs[i]['mag'] = nsha18_ml2mw(rec['mag'])

# load atten coeffs
coeffs = pickle.load(open('atten_coeffs.pkl', 'rb' ))

###############################################################################
# set data to use
###############################################################################

# remove bad recs
keep_nets = set(['AU', 'IU', 'S1', 'G', 'MEL', 'ME', '2O', 'AD', 'SR', 'UM', 'AB', 'VI', 'GM' \
                 '1P', '1P', '2P', '6F', '7K', '7G', 'G', '7B', '4N', '7D', '', 'OZ'])

# get stas to ignore
ignore_stas = open('sta_ignore.txt').readlines()
ignore_stas = set([x.strip() for x in ignore_stas])

####################################################################################
# start main
####################################################################################
fig = plt.figure(1, figsize=(18,11))

# loop through & plot station  data
cs = get_mpl2_colourlist()
events = unique(dictlist2array(recs, 'ev'))
stations = unique(dictlist2array(recs, 'sta'))

'''
# make event list - needed once only
outtxt = 'ORIGIN,MINF,MAXF\n'
for event in events:
    outtxt += ','.join((event, '1.0,7.0')) + '\n'

f = open('event_filtering_lookup.hold', 'w')
f.write(outtxt)
f.close()
crash
'''

# get kappas
kapdat = parse_kappa_data()

# get event filters
filtdat = parse_filtering_data()

sp = 0
ii = 1	
for event in events[::-1]: #[-2:-1]:
    print(event)
    sp += 1
    plt.subplot(2,3,sp)
    log_stack_logfds = []
    handles1 = []
    labels1 = []
    
    ceil_log_fds = -99
    floor_log_fds = 99
    
    i = 0
    for rec in recs:
        if rec['net'] in keep_nets:
            if not rec['sta'] in ignore_stas:
                if len(rec['channels']) > 0:
                    if rec['ev'] == event:
                        print('   '+rec['sta'])
                        cor_fds, cor_fds_nan, freqs = correct_atten(rec, coeffs, kapdat)
                        
                        # get max/min
                        if ceil(max(log10(cor_fds[32:-64]))) > ceil_log_fds:
                            ceil_log_fds = ceil(max(log10(cor_fds[32:-64])))
                            
                        if floor(min(log10(cor_fds[32:-64]))) < floor_log_fds:
                            floor_log_fds = floor(min(log10(cor_fds[32:-64])))
                        
                        if i <= 9:
                            h1, = plt.loglog(freqs, cor_fds_nan,'-', c=cs[i], lw=1, label=rec['sta'])
                        elif i <= 19:
                            h1, = plt.loglog(freqs, cor_fds_nan,'--', c=cs[i-10], lw=1, label=rec['sta'])
                        elif i <= 29:
                            h1, = plt.loglog(freqs, cor_fds_nan,'-.', c=cs[i-20], lw=1, label=rec['sta'])
                        elif i <= 39:
                            linestyle = (0, (3, 5, 1, 5, 1, 5))
                            h1, = plt.loglog(freqs, cor_fds_nan, linestyle=linestyle, c=cs[i-30], lw=1, label=rec['sta'])
                        else:
                            linestyle = (0, (3, 5, 1, 5))
                            h1, = plt.loglog(freqs, cor_fds_nan, linestyle=linestyle, c=cs[i-40], lw=1) # don't write sta
                            if i == 49:
                                i = 40
                        
                        handles1.append(h1)
                        labels1.append(rec['sta'])
                        
                        # stack
                        if log_stack_logfds == []:
                            log_stack_logfds = log(cor_fds_nan)
                        else:
                            log_stack_logfds = vstack((log_stack_logfds, log(cor_fds_nan)))
                            
                        i += 1
                        
                        '''
                        if rec['sta'] == 'ONGER':
                        	plt.subplot(2,3,2)
                        	channel = rec['channels'][0]
                        	raw_fds = rec[channel]['swave_spec']
                        	freqs = rec[channel]['freqs']
                        	
                        	plt.subplot(2,3,2)
                        	plt.loglog(freqs, raw_fds,'b-')
                        	plt.subplot(2,3,sp)
                        '''	
    
    leg1 = plt.legend(handles=handles1, loc=3, fontsize=6, ncol=4)

    # get mean spectra
    if log_stack_logfds != []:
        if len(log_stack_logfds.shape) == 1:
            mean_fds = exp(log_stack_logfds)
        else:
            mean_fds = exp(nanmean(log_stack_logfds, axis=0))
            
        h2, = plt.loglog(freqs, mean_fds,'--', color='0.2', lw=1.5, label='Mean Source Spectrum')
    
        # get upper & lower f for filtering
        for fdat in filtdat:
            if fdat['ev'] == event:
                minf = fdat['minf']
                maxf = fdat['maxf']
        
        # fit mean curve
        fidx = where((freqs >= minf) & (freqs <= maxf) & (isnan(mean_fds) == False))[0]
        #sfidx = where((freqs >= minf) & (freqs <= maxsf))[0] # for labelling curves
        
        data = odrpack.RealData(freqs[fidx], log(mean_fds[fidx]))
        fitted_brune = odrpack.Model(fit_brune_model)
        odr = odrpack.ODR(data, fitted_brune, beta0=[1E-2,1.])
        odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq
        out = odr.run()
        
        omega0 = out.beta[0]
        f0 = abs(out.beta[1])
        print('f0', f0)
        
        m0 = C * omega0
        mw =  m02mw(m0)
        print('Mw', m02mw(m0))
        
        # calc stress drop
        r0 = 2.34 * vsm / (2 * pi * f0)
        
        sd = 7. * m0 / (16. * r0**3) / 10**6 # in MPa
        print('SD', sd, 'MPa')
        print('f0', f0, 'Hz\n' )
    
        # plot fitted curve
        fitted_curve = omega0 / (1 + (freqs / f0)**2)
        h3, = plt.loglog(freqs, fitted_curve, 'k-', lw=1.5, label='Fitted Brune Model')
        plt.legend(handles=[h2, h3], loc=1, fontsize=8)
        plt.title('; '.join((event, 'MW '+str('%0.2f' % mw), 'SD '+str('%0.2f' % sd)+' MPa')), fontsize=10)
        
        if sp == 1 or sp == 4:
           plt.ylabel('Fourier Displacement Spectra (m-s)')
        if sp >= 4:
           plt.xlabel('Frequency (Hz)')
        #plt.title(' - '.join((ev, 'M'+str('%0.2f' % mag), place)), fontsize=10)
    
    plt.gca().add_artist(leg1)
    if f0 < 0.7:
        plt.xlim([0.1, 20])
    else:
        plt.xlim([0.3, 20])
        
    plt.ylim([10**(ceil_log_fds-4), 10**ceil_log_fds])
    plt.grid(which='both', color='0.75')
    
    
    if sp == 6:
        plt.savefig('brune_fit/brune_fit_'+str(ii)+'.png', fmt='png', bbox_inches='tight')
        sp = 0
        ii += 1
        fig = plt.figure(ii, figsize=(18,11))
        
plt.savefig('brune_fit/brune_fit_'+str(ii)+'.png', fmt='png', dpi=150, bbox_inches='tight')
plt.show()    


    