import pickle
from numpy import unique, array, arange, log, log10, exp, mean, nanmean, ndarray, std, sqrt, \
                  nanmedian, nanstd, vstack, pi, nan, isnan, interp, where, zeros_like, ones_like, floor, ceil, \
                  argsort
from scipy.stats import linregress, trim_mean
from misc_tools import get_binned_stats, dictlist2array, get_mpl2_colourlist
from mag_tools import nsha18_mb2mw, nsha18_ml2mw
from mag_tools import m02mw
from get_mag_dist_terms_swan import get_distance_term
#from js_codes import get_average_Q_list, extract_Q_freq
from mapping_tools import distance
from scipy.odr import Data, Model, ODR, models
import scipy.odr.odrpack as odrpack
#from ltsfit import lts_linefit
import matplotlib.pyplot as plt
import matplotlib as mpl
from obspy import UTCDateTime
mpl.style.use('classic')
plt.rcParams['pdf.fonttype'] = 42
import warnings
warnings.filterwarnings("ignore")

import shapefile
from shapely.geometry import Point, Polygon
from mapping_tools import get_field_data

shpfile = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/zones/shapefiles/NSHA13_Background/NSHA13_BACKGROUND_NSHA18_May2016.shp'

sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))

zone_code = get_field_data(sf, 'CODE', 'str')
zone_trt = get_field_data(sf, 'TRT', 'str')
used_zones = set(['EBGZ', 'CBGZ', 'NCCZ'])

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
    
    '''
    # get site kappa - only use mean kappa as geetting spectral correction
    for kap in kapdat:
        if kap['sta'] == rec['sta']:
            if not isnan(kap['kappa0']):
                kappa = kap['kappa0'] # + kap['kappa_r'] * rec['rhyp'] 
                #print(kap['kappa0'])
    '''
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


    
def fit_brune_model_fixed_omega(c, f):
    from numpy import array, log
    '''
    c[0] = omega0
    c[1] = f0
    f    = frequency
    '''
    
    fixed_omega = 2.6 # from dist corrected stacked spectra
    
    # set constants
    vs = 3.6 # km/s
    vsm = vs*1000.
    rho = 2800 # kg/m^3
    C = 4. * pi * rho * vsm**3 * 1000. / (0.55 * 2.0 * 0.71)
    #f = array(exp(logf))
    #print(f

    # fit curve
    FittedCurve = log(fixed_omega / (1 + (f / (c[0]))**2))
    #FittedCurve = C * omega / (1 + (f / c[1])**2)
    
    return FittedCurve
    
def fit_brune_model_fixed_omega_petermann(c, f):
    from numpy import array, log
    '''
    c[0] = omega0
    c[1] = f0
    f    = frequency
    '''
    
    fixed_omega = 0.32 # from dist corrected stacked spectra
    
    # set constants
    vs = 3.6 # km/s
    vsm = vs*1000.
    rho = 2800 # kg/m^3
    C = 4. * pi * rho * vsm**3 * 1000. / (0.55 * 2.0 * 0.71)
    #f = array(exp(logf))
    #print(f

    # fit curve
    FittedCurve = log(fixed_omega / (1 + (f / (c[0]))**2))
    #FittedCurve = C * omega / (1 + (f / c[1])**2)
    
    return FittedCurve
    
def fit_brune_model_fixed_omega_marblebar(c, f):
    from numpy import array, log
    '''
    c[0] = omega0
    c[1] = f0
    f    = frequency
    '''
    
    fixed_omega = 0.018 # from dist corrected stacked spectra
    
    # set constants
    vs = 3.6 # km/s
    vsm = vs*1000.
    rho = 2800 # kg/m^3
    C = 4. * pi * rho * vsm**3 * 1000. / (0.55 * 2.0 * 0.71)
    #f = array(exp(logf))
    #print(f

    # fit curve
    FittedCurve = log(fixed_omega / (1 + (f / (c[0]))**2))
    
    return FittedCurve
    
def fit_brune_model_fixed_omega_2020bowen(c, f):
    from numpy import array, log
    '''
    c[0] = omega0
    c[1] = f0
    f    = frequency
    '''
    
    fixed_omega = 0.0081 # from dist corrected stacked spectra
    
    # set constants
    vs = 3.6 # km/s
    vsm = vs*1000.
    rho = 2800 # kg/m^3
    C = 4. * pi * rho * vsm**3 * 1000. / (0.55 * 2.0 * 0.71)
    #f = array(exp(logf))
    #print(f

    # fit curve
    FittedCurve = log(fixed_omega / (1 + (f / (c[0]))**2))
    
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
            recs[i]['mag'] -= 0.05
        
        # now fix ML
        recs[i]['mag'] = nsha18_ml2mw(rec['mag'])

# load atten coeffs
coeffs = pickle.load(open('atten_coeffs.pkl', 'rb' ))

###############################################################################
# set data to use
###############################################################################

# remove bad recs
keep_nets = set(['AU', 'IU', 'S1', 'II', 'G', 'MEL', 'ME', '2O', 'AD', 'SR', 'UM', 'AB', 'VI', 'GM' \
                 '1P', '1P', '2P', '6F', '7K', '7G', 'G', '7B', '4N', '7D', '', 'OZ', 'OA', 'WG', 'XX'])
# get stas to ignore
ignore_stas = open('sta_ignore.test').readlines()
ignore_stas = set([x.strip() for x in ignore_stas])

'''
SWN20, SWNNG

'''

####################################################################################
# start main
####################################################################################
print('DO NOT EXCLUDE LOW F DATA, BUT DO NOT ALLOW FITTING FROM 0.08 - 0.4 HZ')
fig = plt.figure(1, figsize=(18,11))

# loop through & plot station  data
cs = get_mpl2_colourlist()
events = unique(dictlist2array(recs, 'ev'))
stations = unique(dictlist2array(recs, 'sta'))
rhyps = dictlist2array(recs, 'rhyp')

# sort by rhyps
rhyp_sort_idx = argsort(rhyps)

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

events_dict = []

sp = 0
ii = 1	
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

    #print(','.join((event,str(minf),str(maxf))))

    sp += 1
    plt.subplot(2,3,sp)
    log_stack_logfds = []
    handles1 = []
    labels1 = []
    
    ceil_log_fds = -99
    floor_log_fds = 99
    
    i = 0
    pltboxes = True
    for rsi in rhyp_sort_idx:
        rec = recs[rsi]
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
                            
                            # get max/min
                            if ceil(max(log10(cor_fds[30:120]))) > ceil_log_fds:
                                ceil_log_fds = ceil(max(log10(cor_fds[30:120])))
                                
                                if max(cor_fds[30:120]) < 2*10**(ceil_log_fds-1):
                                    ceil_log_fds = log10(2*10**(ceil_log_fds-1))
                                    
                                elif max(cor_fds[30:120]) < 4*10**(ceil_log_fds-1):
                                    ceil_log_fds = log10(4*10**(ceil_log_fds-1))
                                
                            if floor(min(log10(cor_fds[30:120]))) < floor_log_fds:
                                floor_log_fds = floor(min(log10(cor_fds[30:120])))
                                
                            if i <= 9:
                                h1, = plt.loglog(freqs, cor_fds_nan,'-', c=cs[i], lw=1, label='-'.join((rec['sta'],str('%0.0f' % rec['rhyp']))))
                            elif i <= 19:
                                h1, = plt.loglog(freqs, cor_fds_nan,'--', c=cs[i-10], lw=1, label='-'.join((rec['sta'],str('%0.0f' % rec['rhyp']))))
                            elif i <= 29:
                                h1, = plt.loglog(freqs, cor_fds_nan,'-.', c=cs[i-20], lw=1, label='-'.join((rec['sta'],str('%0.0f' % rec['rhyp']))))
                            elif i <= 39:
                                linestyle = (0, (3, 5, 1, 5, 1, 5))
                                h1, = plt.loglog(freqs, cor_fds_nan, linestyle=linestyle, c=cs[i-30], lw=1, label='-'.join((rec['sta'],str('%0.0f' % rec['rhyp']))))
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
                        
                        # set evmag
                        evmag = rec['mag']
                        evmagtype = rec['magType']
                        evlon = rec['eqlo']
                        evlat = rec['eqla']
                        evdep = rec['eqdp']
                        evdt = rec['ev']
                        #evmb = rec['mb']
    
    leg1 = plt.legend(handles=handles1, loc=3, fontsize=6, ncol=4)
	
    # get mean spectra
    if log_stack_logfds != []:
        if len(log_stack_logfds.shape) == 1:
            mean_fds = exp(log_stack_logfds)
        else:
            mean_fds = exp(nanmean(log_stack_logfds, axis=0))
        
        if log_stack_logfds.shape[0] <= 2:
            qual = 2
        
        h2, = plt.loglog(freqs, mean_fds,'--', color='0.2', lw=1.5, label='Mean Source Spectrum')
        
        # don't fit these freqs
        idx = where((freqs >= 0.08) & (freqs <= 0.4))[0]
        mean_fds[idx] = nan
        
        # fit mean curve
        fidx = where((freqs >= minf) & (freqs <= maxf) & (isnan(mean_fds) == False))[0]
        
        #sfidx = where((freqs >= minf) & (freqs <= maxsf))[0] # for labelling curves
        
        data = odrpack.RealData(freqs[fidx], log(mean_fds[fidx]))
        
        # do special case for Broome
        if event == UTCDateTime('2019-07-14T05:39:24.991000Z'):
            print('M6.6 Broome')
            fitted_brune = odrpack.Model(fit_brune_model_fixed_omega)
            odr = odrpack.ODR(data, fitted_brune, beta0=[0.3])
            odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq
            out = odr.run()
            
            fixed_omega = 2.6
            omega0 = fixed_omega
            f0 = abs(out.beta[0])
            print('f0', f0)
        
        elif event == UTCDateTime('2016-05-20T18:14:02.000000Z'):
            print('M6.0 Petermann')
            fitted_brune = odrpack.Model(fit_brune_model_fixed_omega_petermann)
            odr = odrpack.ODR(data, fitted_brune, beta0=[0.3])
            odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq
            out = odr.run()
            
            fixed_omega = 0.32
            omega0 = fixed_omega
            f0 = abs(out.beta[0])
            print('f0', f0)    
        elif event == UTCDateTime('2021-11-13T13:05:52.663000Z'):
            print('M5.3 Marble Bar')
            fitted_brune = odrpack.Model(fit_brune_model_fixed_omega_marblebar)
            odr = odrpack.ODR(data, fitted_brune, beta0=[0.3])
            odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq
            out = odr.run()
            
            fixed_omega = 0.018
            omega0 = fixed_omega
            f0 = abs(out.beta[0])
            print('f0', f0)    
        elif event == UTCDateTime('2020-04-15T07:11:04.955000Z'):
            print('2020 Bowen')
            fitted_brune = odrpack.Model(fit_brune_model_fixed_omega_marblebar)
            odr = odrpack.ODR(data, fitted_brune, beta0=[0.3])
            odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq
            out = odr.run()
            
            fixed_omega = 0.0081
            omega0 = fixed_omega
            f0 = abs(out.beta[0])
            print('f0', f0) 
        else:
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
        
        # add to events
        edict = {}
        edict['evstr'] = event
        edict['evdt'] = evdt
        edict['lon'] = evlon
        edict['lat'] = evlat
        edict['dep'] = evdep
        edict['omag'] = evmag
        edict['omag_type'] = evmagtype
        #edict['mb'] = evmb
        edict['brune_mw'] = mw
        edict['brune_sd'] = sd
        edict['brune_f0'] = f0
        edict['minf'] = minf
        edict['maxf'] = maxf
        edict['qual'] = qual
        edict['stas'] = labels1
        edict['sta_spectra'] = log_stack_logfds
        edict['freqs'] = freqs
        
        
        if log_stack_logfds.shape[0] == 150:
            nrecs = 1
        else:
            nrecs = log_stack_logfds.shape[0]
        edict['nrecs'] = nrecs
        
        events_dict.append(edict)
    
        # plot fitted curve
        fitted_curve = omega0 / (1 + (freqs / f0)**2)
        h3, = plt.loglog(freqs, fitted_curve, 'k-', lw=1.5, label='Fitted Brune Model')
        plt.legend(handles=[h2, h3], loc=1, fontsize=8)
        
        edict['fitted_spectra'] = fitted_curve
        
        if qual == 0:
            plt.title('; '.join((evmagtype+str('%0.1f' % evmag), str(event)[0:16], 'MW '+str('%0.2f' % mw), 'SD '+str('%0.2f' % sd)+' MPa')), fontsize=10, color='red')
        else:
            plt.title('; '.join((evmagtype+str('%0.1f' % evmag), str(event)[0:16], 'MW '+str('%0.2f' % mw), 'SD '+str('%0.2f' % sd)+' MPa')), fontsize=10, color='k')
        
        if sp == 1 or sp == 4:
           plt.ylabel('Fourier Displacement Spectra (m-s)')
        if sp >= 4:
           plt.xlabel('Frequency (Hz)')
        #plt.title(' - '.join((ev, 'M'+str('%0.2f' % mag), place)), fontsize=10)
    
    plt.gca().add_artist(leg1)
    if f0 < 1. or omega0 >= 0.005:
        plt.xlim([0.03, 20])
    else:
        plt.xlim([0.1, 20])
    
    # set ylims based on omega0
    ceil_log_fds = log10(omega0) + 1.1
    ymin = 10**(ceil_log_fds-4)  
    ymax = 10**ceil_log_fds
    plt.ylim([ymin, ymax])
    
    # plot ignored freqs
    plt.fill([0.01, 0.01, minf, minf, 0.01], [ymin, ymax, ymax, ymin, ymin], '0.9', ec='0.9', zorder=0)
    plt.fill([maxf, maxf, 20, 20, maxf], [ymin, ymax, ymax, ymin, ymin], '0.9', ec='0.9', zorder=1)
                            
    plt.grid(which='both', color='0.7')
        
    if sp == 6:
        plt.savefig('brune_fit/brune_fit_'+str(ii)+'.png', fmt='png', bbox_inches='tight')
        sp = 0
        ii += 1
        fig = plt.figure(ii, figsize=(18,11))

plt.savefig('brune_fit/brune_fit_'+str(ii)+'.png', fmt='png', dpi=150, bbox_inches='tight')


##########################################################################################

# export Brune data
pklfile = open('brune_data.pkl', 'wb')
pickle.dump(events_dict, pklfile, protocol=-1)
pklfile.close()

# write csv
txt = 'EVENT,LON,LAT,DEP,OMAG,OMAG_TYPE,BRUNE_MAG,STRESS_DROP,CORN_FREQ,NRECS,FMIN,FMAX,QUALITY\n'

for ev in events_dict:
    txt += ','.join((str(ev['evdt']),str(ev['lon']),str(ev['lat']),str(ev['dep']),str(ev['omag']),ev['omag_type'], \
                     str(ev['brune_mw']),str(ev['brune_sd']),str(ev['brune_f0']),str(ev['nrecs']), str(ev['minf']),str(ev['maxf']),str(ev['qual']))) + '\n'

# write to file
f = open('brune_stats.csv', 'w')
f.write(txt)
f.close()

# now show figs        
plt.show()    

    