import pickle
from numpy import unique, array, arange, log, log10, exp, mean, nanmean, ndarray, \
                  nanmedian, vstack, pi, nan, isnan, interp, where, zeros_like, sqrt
from misc_tools import get_binned_stats, dictlist2array, get_mpl2_colourlist
#from js_codes import get_average_Q_list, extract_Q_freq
from mag_tools import nsha18_mb2mw, nsha18_ml2mw
from mapping_tools import distance
from scipy.odr import Data, Model, ODR, models
import scipy.odr.odrpack as odrpack
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('classic')
plt.rcParams['pdf.fonttype'] = 42
import warnings
warnings.filterwarnings("ignore")
from sys import argv

####################################################################################
# correct Wei et al (2017) attenuation
####################################################################################

def correct_atten(rec, coeffs):
    channel = rec['channels'][0]
    raw_fds = rec[channel]['swave_spec']
    sn_ratio = rec[channel]['sn_ratio']
    freqs = rec[channel]['freqs']
    
    sn_thresh = 5.
    
    # loop thru freqs
    distterms = []
    for c in coeffs:
        # get geometric spreading
        D1 = sqrt(rec['rhyp']**2 + c['nref']**2)
        if rec['rhyp'] <= c['r1']:
            distterm = c['nc0'] * log10(D1) + c['nc1']
        
        # set mid-field
        elif rec['rhyp'] > c['r1'] and rec['rhyp'] <= c['r2']:
            D1 = sqrt(c['r1']**2 + c['nref']**2)
            distterm = c['nc0'] * log10(D1) + c['nc1'] \
                       + c['mc0'] * log10(rec['rhyp'] / c['r1']) + c['mc1'] * (rec['rhyp'] - c['r1'])
        
        # set far-field
        elif rec['rhyp'] > c['r2']:
            D1 = sqrt(c['r1']**2 + c['nref']**2)
            distterm = c['nc0'] * log10(D1) + c['nc1'] \
                       + c['mc0'] * log10(c['r2'] / c['r1']) + c['mc1'] * (c['r2'] - c['r1']) \
                       + c['fc0'] * log10(rec['rhyp'] / c['r2']) + c['fc1'] * (rec['rhyp'] - c['r2'])
        
        distterms.append(distterm)
    
    distterms.append(nan) # as no coeffs for last freq
    
    cor_fds = 10**(log10(raw_fds) - distterms) # - log10(k_term))
    
    # get data exceeding SN ratio
    idx = where(sn_ratio < sn_thresh)[0]
    cor_fds_nan = cor_fds.copy()
    cor_fds_nan[idx] = nan
            
    return cor_fds, cor_fds_nan, freqs
 
def regress_kappa(freqs, mean_fas):
    from numpy import log, isnan, nan, pi
    from scipy.stats import linregress
    
    # regress record kappa
    ridx = where((freqs >= 5.) & (freqs <= 12.) & (isnan(mean_fas) == False))[0] # assume 5Hz past corner
    if len(ridx) <= 5: # need at least 5 data points
        lamda = nan
        kappa = nan
        kappa_intercept = nan
        ridx = []
    else:
        lamda = linregress(freqs[ridx], log(mean_fas[ridx]))
        #print(lamda)
        kappa = -lamda[0] / pi
        kappa_intercept = lamda[1]
            
    return kappa, kappa_intercept, ridx
    
####################################################################################
# set def params
####################################################################################

vs = 3.6 # km/s
vsm = vs*1000.
rho = 2800 # kg/m^3
C = 4. * pi * rho * vsm**3 * 1000. / (0.55 * 2.0 * 0.707)

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

####################################################################################
# start main
####################################################################################
pltTrue = argv[1]
if pltTrue == 'True':
    pltTrue = True
else:
    pltTrue = False

if pltTrue == True:
    fig = plt.figure(1, figsize=(18,11))

# loop through & plot station  data
cs = get_mpl2_colourlist()

handles1 = []
labels1 = []

events = unique(dictlist2array(recs, 'ev'))
stations = unique(dictlist2array(recs, 'sta'))

# load station sets
'''
keep_nets = set(['AU', 'IU', 'S1', 'G', 'MEL', 'ME', '20', 'AD', 'SR', 'UM', 'AB', \
                 '1P', '1K', '1P', '2P', '6F', '7K', '7G', 'G', '7B', '4N', '7D'])
if rec['net'] in keep_nets:
'''

lines = open('station_sets.csv').readlines()
sta_sets = []
for line in lines:
    sta_sets.append(set(line.strip().split(',')))

kappa_txt = 'STA,KAPPA,CNT\n'
kappa_list = []

i = 1
ii = 1
for sta in stations:
    cnt = 0
    ax = plt.subplot(3,3,i)
    
    log_stack_logfas = []
    
    # check station sets
    sta_set = set([sta])
    
    for ss in sta_sets:
        if sta in ss:
            sta_set = ss
    
    print(sta_set)

    # get all records for each sta
    for rec in recs:
        if rec['sta'] in sta_set and len(rec['channels']) > 0:
            
            print(rec['ev'] + ' ' +rec['sta'])
            # get corrected fas
            cor_fds, cor_fds_nan, freqs = correct_atten(rec, coeffs)
            cor_fas = 2 * pi * freqs * cor_fds
            cor_fas_nan = cor_fds_nan * (2 * pi * freqs)**2
            
            # normalise at 6.3 Hz
            if isnan(cor_fas_nan[72]):
                norm_fas = zeros_like(cor_fas_nan) * nan
            else:
                norm_fas = cor_fas_nan / cor_fas_nan[72]
                
            norm_fas_show = cor_fas / cor_fas[72] # for data that does not meet S/N thresh
                
            # stack
            if log_stack_logfas == []:
                log_stack_logfas = log(norm_fas)
            else:
                log_stack_logfas = vstack((log_stack_logfas, log(norm_fas)))
            
            # plot fas
            if pltTrue == True:
                plt.semilogy(freqs, norm_fas_show, '-', c='0.6', lw=0.4)
                plt.semilogy(freqs, norm_fas, 'r-', lw=0.4)
                #plt.loglog(freqs, norm_fas, '-', c='0.6', lw=0.4)
            
            if isinstance(norm_fas, ndarray):
                cnt += 1
    
    # get mean spectra
    mean_fas = exp(nanmean(log_stack_logfas, axis=0))
    if isinstance(mean_fas, ndarray):
        if pltTrue == True:
            plt.semilogy(freqs, mean_fas, 'k-', lw=1.5)
    
    elif len(log_stack_logfas) == 0:
        ignore=0
    
    elif len(log_stack_logfas.shape) == 1:
        mean_fas = exp(log_stack_logfas)
        if pltTrue == True:
            plt.semilogy(freqs, mean_fas, 'k-', lw=1.5)
    
    nidx = where(isnan(mean_fas) == False)[0]
    # calculate kappa
    if len(nidx) >= 5:
        kappa, kappa_intercept, ridx = regress_kappa(freqs, mean_fas)
        kappa_txt += ','.join((sta, str('%0.6f' % kappa), str(cnt))) + '\n'
        kappa_list.append(kappa)
        print('kappa = ' + str('%0.4f' % kappa))
        
        # plot kappa
        if pltTrue == True:
            kap_plt = exp(-1 * pi * freqs * kappa + kappa_intercept)
            plt.semilogy(freqs[ridx], kap_plt[ridx], 'g-', lw=2)
    
    if pltTrue == True:
        plt.title(sta+' '+str(cnt))
        plt.xlim([0, 20])
        # prep for next subplot
        i += 1
        if i > 9:
            plt.savefig('kappa/site_kappa_'+str(ii)+'.png', fmt='png', bbox_inches='tight')
            i = 1
            ii += 1
            fig = plt.figure(ii, figsize=(18,11))
    
    print(sta+' '+str(cnt))
    
if pltTrue == True:
    plt.savefig('kappa/site_kappa_'+str(ii)+'.png', fmt='png', bbox_inches='tight')
    

# add mean kappa to list
mean_kappa = nanmean(array(kappa_list))
kappa_txt += ','.join(('MEAN_SITE', str('%0.6f' % mean_kappa), 'nan')) + '\n'

# write csv
f = open('site_kappa.csv', 'w')
f.write(kappa_txt)
f.close()