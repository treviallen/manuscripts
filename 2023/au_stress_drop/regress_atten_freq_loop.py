import pickle
from numpy import unique, array, arange, log, log10, logspace, exp, mean, nanmean, ndarray, \
                  nanmedian, hstack, pi, nan, isnan, interp, polyfit, where, zeros_like, polyfit, sqrt
from misc_tools import get_binned_stats, dictlist2array, get_mpl2_colourlist
from mag_tools import nsha18_mb2mw, nsha18_ml2mw
#from js_codes import get_average_Q_list, extract_Q_freq
#from mapping_tools import distance
from scipy.odr import Data, Model, ODR, models
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import linregress
import scipy.odr.odrpack as odrpack
from obspy import UTCDateTime
mpl.style.use('classic')
plt.rcParams['pdf.fonttype'] = 42
import warnings
warnings.filterwarnings("ignore")

def highside(x, hx):
    from numpy import zeros_like
    xmod = zeros_like(x)

    idx = x >= hx
    xmod[idx] = 1
    return xmod

def lowside(x, hx):
    from numpy import zeros_like
    xmod = zeros_like(x)

    idx = x <= hx
    xmod[idx] = 1
    return xmod

def bilinear_reg_free(c, x):
    """ sets up bilinear equation with a free hinge position"""
    from numpy import zeros_like
    hx = c[3] # hinge magnitude
    ans2 = zeros_like(x)
    ans1 = zeros_like(x)

    #idx1 = x <= hx
    #idx2 = x >= hx

    modx_lo = lowside(x, hx)
    modx_hi = highside(x, hx)

    ans1 = modx_lo * (c[0] * x + c[1])
    yhinge = c[0] * hx + c[1]
    ans2 = modx_hi * (c[2] * (x-hx) + yhinge)

    return ans1 + ans2

def bilinear_Q_regression(freqs_log, Q_list_log):
    freqs_log = array(freqs_log)
    
    #data = odrpack.RealData(x[12:], y[12:])
    data = odrpack.RealData(freqs_log, Q_list_log)

    #x_range = arange(-0.5, 1.0, step=0.01) # x range 

    bilin_reg = odrpack.Model(bilinear_reg_free)
    odr = odrpack.ODR(data, bilin_reg, beta0=[0.4, 3.0, 0.3, -0.5])

    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as least squares
    out = odr.run()
    #print('\nbilinear auto\n')
    #out.pprint()

    a = out.beta[0]
    b = out.beta[1]
    c = out.beta[2]
    hx = out.beta[3] # x hinge point

    log_q_fit = b + a * freqs_log # get y values from bilinear
    yhinge = b + a * hx
    idx = freqs_log > hx
    log_q_fit[idx] = c * (freqs_log[idx]-hx) + yhinge

    return log_q_fit
    
# get near-source fit
nref = 10.
def fit_near_source_saturation(c, x):
    from numpy import sqrt
    
    D = sqrt((10**x)**2 + nref**2) # were c1 ~= 5-10?
    #D = sqrt((10**x)**2 + nref**2)
    
    ans = c[0] * log10(D) + c[1]
    
    return ans


            
###############################################################################
# grunt defs
###############################################################################

def normalise_data(p, recs):
    log_norm_amps = []

    chan = recs[0]['channels'][-1]
    
    print("Reg Freq = " +str('%0.3f' % freq))
    
    i = 1
    for m in mrng:
        cnt = 0
        
        mrhyps = []
        mamps  = []
        mmags = []
    
        # get all records for each sta
        for rec in recs:
            if len(rec['channels']) > 0 and rec['mag'] >= m-0.05 and rec['mag'] < m+0.05:
                if rec['net'] in keep_nets:
                    channel = rec['channels'][0]
                    
                    if rec[channel]['sn_ratio'][fidx[p]] >= 5.:
                        mrhyps.append(rec['rhyp'])
                        mmags.append(rec['mag'])
                        mamps.append(rec[channel]['swave_spec'][fidx[p]])
        
        mrhyps = array(mrhyps)
        mamps = array(mamps)
        mmags = array(mmags)
        
        if len(mrhyps) > 0:
            
            i += 1
            
            # get binned data
            logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log10(mrhyps), log10(mamps))
            #print (m, 10**binstrp)
            # normalise data @ 500 km
            nidx = where((binstrp > 2.69) & (binstrp < 2.71))[0]
            nidx = where((binstrp > 2.49) & (binstrp < 2.51))[0]
            
            if len(nidx) > 0:
                #print (m, nidx)
                namps = log10(mamps) - logmedamp[nidx]
                
                if len(log_norm_amps) == 0:
                    log_norm_amps = namps
                    logamps = log10(mamps)
                    norm_rhyps = mrhyps
                    mags = mmags
                else:
                    log_norm_amps = hstack((log_norm_amps, namps))
                    logamps = hstack((logamps, log10(mamps)))
                    norm_rhyps = hstack((norm_rhyps, mrhyps))
                    mags = hstack((mags, mmags))
    
    return log_norm_amps, norm_rhyps, mags, logamps
    
def fit_near_source_atten(plt, norm_rhyps, log_norm_amps, r1):
    if pltTrue == True:
        plt.loglog(norm_rhyps, 10**log_norm_amps, '+', c='0.6', lw=0.5, ms=6)
        plt.xlim([4, 2250])
        #plt.ylim([5E-3, 500])
    
    # get binned data
    log_norm_rhyps = log10(norm_rhyps)
    logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log_norm_rhyps, log_norm_amps)
    if pltTrue == True:
        plt.loglog(10**medx, 10**logmedamp, 'rs', ms=7)

    # fit all data - keep me!
    didx = where((medx >= log10(minr)) & (medx < log10(r1)) & (nperbin > 2))[0] 
    data = odrpack.RealData(medx[didx], logmedamp[didx])
    
    #didx = where((log_norm_rhyps >= log10(minr)) & (log_norm_rhyps <= log10(r1)))[0]
    #data = odrpack.RealData(log10(norm_rhyps[didx]), log_norm_amps[didx])
    
    afit = odrpack.Model(fit_near_source_saturation)
    #odr = odrpack.ODR(data, afit, beta0=[-1., 2, 1.])
    odr = odrpack.ODR(data, afit, beta0=[-1., 2])
    
    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
    out = odr.run()
    nc = out.beta
    
    return nc

#nref2 = 180.    
def fit_mid_field_atten(plt, nc, norm_rhyps, log_norm_amps):
    # fit all data
    log_norm_rhyps = log10(norm_rhyps)
    
    logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log_norm_rhyps, log_norm_amps)
    
    # fit all data - keep me!
    didx = where((medx >= log10(r1)) & (medx < log10(r2)) & (nperbin > 2))[0]
    data = odrpack.RealData(medx[didx], logmedamp[didx])
    
    #didx = where((log_norm_rhyps >= log10(minr)) & (log_norm_rhyps < log10(minr)))[0] # & (nperbin > 2))[0] #log10(r1))[0]
    #data = odrpack.RealData(log_norm_rhyps[didx], log_norm_amps[didx])

    afit = odrpack.Model(fit_mid_field)
    odr = odrpack.ODR(data, afit, beta0=[0., -0.001])
    
    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
    out = odr.run()
    mc = out.beta
    
    # plot
    xrng_ff = arange(log10(r1), log10(r2), 0.02)
    #D = sqrt(r1**2 + nc[2]**2)
    D1 = sqrt(r1**2 + nref**2)
    yrng_ff = nc[0] * log10(D1) + nc[1] + mc[0] * log10(10**xrng_ff / r1) + mc[1] * (10**xrng_ff - r1)
        
    if pltTrue == True:
       plt.loglog(10**xrng_ff, 10**yrng_ff, 'g-', lw=2)
    
    return mc
    
def fit_far_field_atten(plt, nc, mc, norm_rhyps, log_norm_amps):
    # fit all data
    log_norm_rhyps = log10(norm_rhyps)
    
    logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log_norm_rhyps, log_norm_amps)
    
    # fit all data - keep me!
    didx = where((medx >= log10(r2)) & (medx < log10(r3)) & (nperbin > 2))[0] 
    data = odrpack.RealData(medx[didx], logmedamp[didx])
    
    #didx = where((log_norm_rhyps >= log10(minr)) & (log_norm_rhyps < log10(minr)))[0] # & (nperbin > 2))[0] #log10(r1))[0]
    #data = odrpack.RealData(log_norm_rhyps[didx], log_norm_amps[didx])

    afit = odrpack.Model(fit_far_field)
    odr = odrpack.ODR(data, afit, beta0=[-1., -0.001])
    
    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
    out = odr.run()
    fc = out.beta
    
    # plot
    xrng_ff = arange(log10(r2), log10(r3), 0.02)
    D1 = sqrt(r1**2 + nref**2)
    yrng_ff = nc[0] * log10(D1) + nc[1] \
              + mc[0] * log10(r2 / r1) + mc[1] * (r2 - r1) \
              + fc[0] * log10(10**xrng_ff / r2) + fc[1] * (10**xrng_ff - r2)
        
    if pltTrue == True:
        plt.loglog(10**xrng_ff, 10**yrng_ff, 'b-', lw=2)
    
    return fc

def get_distance_residuals(nc, mc, fc):
     D1 = sqrt((10**x)**2 + nref**2)
     ans = nc[0] * log10(D1) + nc[1]
     
     # set mid-field
     idx = where((10**x > r1) & (10**x <= r2))[0]
     D1 = sqrt(r1**2 + nref**2)
     ans[idx] = nc[0] * log10(D1) + nc[1] + mc[0] * log10(10**x[idx] / r1) + mc[1] * (10**x[idx] - r1)
     
     idx = where(10**x > r2)[0]
     ans[idx] = nc[0] * log10(D1) + nc[1] \
                + mc[0] * log10(r2 / r1) + mc[1] * (r2 - r1) \
                + fc[0] * log10(10**x[idx] / r2) + fc[1] * (10**x[idx] - r2)
   
###############################################################################
# load pickle 
###############################################################################

recs = pickle.load(open('fft_data.pkl', 'rb' ))

# remove bad recs
keep_nets = set(['AU', 'IU', 'S1', 'G', 'MEL', 'ME', '20', 'AD', 'SR', 'UM', 'OA', \
                 '1P', 'OA', '1K', '1P', '2P', '6F', '7K', '7G', 'G', '7B', '4N'])

####################################################################################
# start main
####################################################################################

events = unique(dictlist2array(recs, 'ev'))
mags = dictlist2array(recs, 'mag')
magTypes = dictlist2array(recs, 'magType')
rhyp = dictlist2array(recs, 'rhyp')
datetimes = dictlist2array(recs, 'ev')

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
maxDist = 2200
minRegDist = 100
maxRegDist = 1000

bins = arange(log10(minDist), log10(maxDist), 0.1)
    
# loop through freqs
rec = recs[0]
chan = rec['channels'][-1]
freqs = rec[chan]['freqs']
fidx=arange(0,100,10)+6 # for testing
fidx=arange(14,95,2) # for regressing all coeffs

pltTrue = True
if len(fidx) > 12:
    pltTrue = False

if pltTrue == True:
    fig = plt.figure(1, figsize=(18,11))
    
coeffs = []

for p, freq in enumerate(freqs[fidx]):
    if pltTrue == True:
        fig = plt.figure(1, figsize=(18,11))
        plt.subplot(3,4,p+1)
        plt.ylabel(str(freq), fontsize=7)
    
    minr = 5.
    r1 = 70 # max dist for near source
    r2 = 600
    r3 = 2100
    
    # get normalised amplitudes
    log_norm_amps, norm_rhyps, mags, logamps = normalise_data(p, recs)
    
    ###############################################################################
    # plt near-field GR
    ###############################################################################
    
    #plt.ylabel('Normalised Spectral Amplitude')
    #plt.xlabel('Hypocentral Distance (km)')
    
    nc = fit_near_source_atten(plt, norm_rhyps, log_norm_amps, r1)
    print(nc)
    
    # plot
    xrng_nf = log10(arange(1, r1+1))
    
    #D = sqrt((10**xrng_nf)**2 + nc[2]**2)
    D1 = sqrt((10**xrng_nf)**2 + nref**2)
    yrng_nf = nc[0] * log10(D1) + nc[1]
    
    if pltTrue == True:
        plt.loglog(10**xrng_nf, 10**yrng_nf, 'k-', lw=2)
    
    ###############################################################################
    # plt mid-field GR
    ###############################################################################
    
    # get far-field fit
    #nref2 = 100.
    def fit_mid_field(c, x):
        from numpy import sqrt, log10
        
        #D = sqrt((10**x)**2 + nc[2]**2) # were c1 ~= 5-10?
        D1 = sqrt((10**x)**2 + nref**2)
        ans = nc[0] * log10(D1) + nc[1]
        
        idx = where((10**x > r1) & (10**x <= r2))[0]
        D1 = sqrt(r1**2 + nref**2)
        
        ans[idx] = nc[0] * log10(D1) + nc[1] + c[0] * log10(10**x[idx] / r1) + c[1] * (10**x[idx] - r1)
            
        return ans
    
    mc = fit_mid_field_atten(plt, nc, norm_rhyps, log_norm_amps)
    print(mc)
    
    ###############################################################################
    # plt far-field GR
    ###############################################################################
    
    # get far-field fit
    #nref2 = 100.
    def fit_far_field(c, x):
        from numpy import sqrt, log10
        
        # set near-field
        D1 = sqrt((10**x)**2 + nref**2)
        ans = nc[0] * log10(D1) + nc[1]
        
        # set mid-field
        idx = where((10**x > r1) & (10**x <= r2))[0]
        D1 = sqrt(r1**2 + nref**2)
        ans[idx] = nc[0] * log10(D1) + nc[1] + mc[0] * log10(10**x[idx] / r1) + mc[1] * (10**x[idx] - r1)
        
        idx = where(10**x > r2)[0]
        ans[idx] = nc[0] * log10(D1) + nc[1] \
                   + mc[0] * log10(r2 / r1) + mc[1] * (r2 - r1) \
                   + c[0] * log10(10**x[idx] / r2) + c[1] * (10**x[idx] - r2)
            
        return ans
    
    fc = fit_far_field_atten(plt, nc, mc, norm_rhyps, log_norm_amps)
    print(fc)
    
    #plt.savefig('norm_geom_spread.png', fmt='png', bbox_inches='tight')
    #plt.show()
    
    ###############################################################################
    # fit mag intercept - mag scaling only needs to be rough as is for testing 
    # instrument response
    ###############################################################################
    
    fig = plt.figure(2, figsize=(18,11))
    def fit_mag_intercept(c, x):
        from numpy import sqrt, log10
        
        D1 = sqrt(x**2 + nref**2) # were c1 ~= 5-10?
        
        ans = nc[0] * log10(D1) + nc[0]
        
        idx = x >= r1
        D1 = sqrt(r1**2 + nc[2]**2) # were c1 ~= 5-10?
        #D = sqrt(r1**2 + nref**2)
        
        ans[idx] = nc[0] * log10(D1) + c[0] + fc[0] * log10(x[idx] / r1) + fc[1] * (x[idx] - r1)
        
        return ans
    
    # get residual
    D1 = sqrt(norm_rhyps**2 + nref**2)
    predamps = nc[0] * log10(D1) + nc[1]
    
    # set mid-field
    idx = where((norm_rhyps > r1) & (norm_rhyps <= r2))[0]
    D1 = sqrt(r1**2 + nref**2)
    predamps[idx] = nc[0] * log10(D1) + nc[1] + mc[0] * log10(norm_rhyps[idx] / r1) + mc[1] * (norm_rhyps[idx] - r1)
    
    # set far-field
    idx = where(norm_rhyps > r2)[0]
    predamps[idx] = nc[0] * log10(D1) + nc[1] \
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
    if pltTrue == True:
        plt.plot(medx, logmedamp, 'rs', ms=6.5)
    
    # fit mag linear
    magc = polyfit(medx, logmedamp, 1)
    
    # plot mag fit
    if pltTrue == True:
        yrng = magc[0] * mrng + magc[1]
        plt.plot(mrng, yrng, 'k-', lw=2)
    
    ###############################################################################
    # make coeffs dict
    ###############################################################################
    
    coeffs.append({'nc0': nc[0], 'nc1': nc[1],
                   'mc0': mc[0], 'mc1': mc[1],
                   'fc0': fc[0], 'fc1': fc[1],
                   'magc0': magc[0], 'magc1': magc[1],
                   'r1': r1, 'r2': r2, 'nref':nref, 'freq':freq})
    
    '''
    # fit mag-dependent offset
     didx = where((medx >= 0.1) & (nperbin > 2))[0] #log10(r1))[0]
    
     data = odrpack.RealData(10**medx[didx], logmedamp[didx])
     
     afit = odrpack.Model(fit_mag_intercept)
     odr = odrpack.ODR(data, afit, beta0=[1.])
     
     odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
     out = odr.run()
     mc = out.beta
     '''
    ###############################################################################
    # build mag regression data
    ###############################################################################
    
    """            
    m_intercept = []
    m_reg = []
    #fig = plt.figure(3, figsize=(18,11))
    i = 1
    for m in mrng:
        cnt = 0
        #ax = plt.subplot(4,5,i)
        
        mrhyps = []
        mamps  = []
    
        # get all records for each sta
        for rec in recs:
            if rec['net'] in keep_nets:
                if len(rec['channels']) > 0 and rec['mag'] >= m-0.05 and rec['mag'] < m+0.05:
                    
                    channel = rec['channels'][0]
                    
                    if rec[channel]['sn_ratio'][fidx[p]] >= 4.:
                        mrhyps.append(rec['rhyp'])
                        mamps.append(rec[channel]['swave_spec'][fidx[p]])
        
        mrhyps = array(mrhyps)
        mamps = array(mamps)
        
        if len(mrhyps) > 0:
            '''
            plt.loglog(mrhyps, mamps, '+', c='0.6', lw=0.5, ms=5)
            plt.xlim([100, 2250])
            plt.ylim([1E-8, 1E-3])
            plt.ylabel(str(round(m,1)))
            i += 1
            '''
            # get binned data
            logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log10(mrhyps), log10(mamps))
            #plt.loglog(10**medx, 10**logmedamp, 'rs', ms=6.5)
            
            # normalise to 500 km
            #namps = log10(mamps) - logmedamp[nidx]
            '''
            if len(log_norm_amps) == 0:
                log_norm_amps = namps
                norm_rhyps = mrhyps
            else:
                log_norm_amps = hstack((log_norm_amps, namps))
                norm_rhyps = hstack((norm_rhyps, mrhyps))
            '''
            
            # fit mag-dependent offset
            didx = where((medx >= 0.1) & (nperbin > 2))[0] #log10(r1))[0]
    
            data = odrpack.RealData(10**medx[didx], logmedamp[didx])
            
            afit = odrpack.Model(fit_mag_intercept)
            odr = odrpack.ODR(data, afit, beta0=[1.])
            
            odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
            out = odr.run()
            mc = out.beta
            
            if round(m,1) != 5.5 and round(m,1) != 6.0 and round(m,1) != 6.2:
                m_intercept.append(mc[0])
                m_reg.append(m)
            
            #plt.loglog(10**xplt, 10**yplt, 'k-', lw=2)
                
    #plt.savefig('fixed_geom_spread.png', fmt='png', bbox_inches='tight')
    #plt.show()
    """
    ###############################################################################
    # get mag scaling
    ###############################################################################
    """
    '''
    fig = plt.figure(4, figsize=(8, 8))
    
    plt.plot(m_reg, m_intercept, 'rs')
    '''
    # linear
    #mag_fit = linregress(array(m_reg), array(m_intercept))
    #yplt = mag_fit[0] * mrng + mag_fit[1]
    
    # quadratic
    m1, m2, m0 = polyfit((array(m_reg)-4.), array(m_intercept), 2)
    yplt = m0 + m1 * (mrng-4.)**2. + m2 * (mrng-4.)
    
    # also fit linear
    m2_lin, m0_lin = polyfit((array(m_reg)-4.), array(m_intercept), 1)
    yplt_lin = m0_lin + m2_lin * (mrng-4.)
    
    # check concavity of slope - use quadratic
    if yplt[0] < yplt_lin[0] and yplt[-1] < yplt_lin[-1]:                    
        '''
        m0_array.append(m0)
        m1_array.append(m1)
        m2_array.append(m2)
        '''
        plt.plot(mrng, yplt, 'k-', lw=2)
    # use linear
    else:
        '''
        m0_array.append(m0_lin)
        m1_array.append(0.0)
        m2_array.append(m2_lin)
        '''
        plt.plot(mrng, yplt_lin, 'k-', lw=2)
    '''
    plt.savefig('mag_scaling.png', fmt='png', bbox_inches='tight')
    plt.show()
    '''
    """
    ###############################################################################
    # get residuals and far-field correction
    ###############################################################################
    #y = m0 + m1 * (mags-4.)**2. + m2 * (mags-4.) + r1 * log10(rhyp)
    """
    yres = [] 
    rhyps = []
    stas = []
    evdt = []
    for rec in recs:
        if rec['net'] in keep_nets:
            try:
                channel = rec['channels'][0]
                    
                if rec[channel]['sn_ratio'][fidx[p]] >= 4.:
                    rhyps.append(rec['rhyp'])
                    
                    #!!! use new atten coeffs !!!
                    if rec['rhyp'] <= r1:
                        D = sqrt(rec['rhyp']**2 + nc[2]**2)
                        y = m0 + m1 * (rec['mag']-4.)**2. + m2 * (rec['mag']-4.) \
                            + nc[0] * log10(D)
                    else:
                        D = sqrt(r1**2 + nc[2]**2)
                        y = m0 + m1 * (rec['mag']-4.)**2. + m2 * (rec['mag']-4.) \
                            + nc[0] * log10(D) + fc[0] * log10(rec['rhyp'] / r1) + fc[1] * (rec['rhyp'] - r1)
                    
                    yobs = log10(rec[channel]['swave_spec'][fidx[p]])
                    yres.append(yobs - y)
                    stas.append(rec['sta'])
                    evdt.append(rec['ev'])
                    
                else:
                    yres.append(nan)
                    rhyps.append(rec['rhyp'])
                    stas.append(rec['sta'])
                    evdt.append(rec['ev'])
            
            except:
                print('No data')
    
    # get binned data   
    bins = arange(log10(minDist), log10(maxDist), 0.1)
    medbin, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log10(rhyps), yres)
    #plt.plot(medx, medbin, 'rs', ms=6.5)
    
    y_at_xhinge = 0
    xhinge = 3.
    xmax = log10(10**3.35)
    def fit_fixed_intercept(c, x):
        '''
        x = array fo x values
        y_at_xmax = y value at xmax
        '''
        #print(xhinge)
        #xmax = 10**logxmax # hardwired in distance atten
        
        ans = c * (x - xhinge) + y_at_xhinge
        
        return ans
    
    ridx = where((binstrp >= xhinge) & (binstrp < xmax))[0]
    data = odrpack.RealData(binstrp[ridx], medbin[ridx])
    
    truncdist = odrpack.Model(fit_fixed_intercept)
    odr = odrpack.ODR(data, truncdist, beta0=[0.])
    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, ODR = 0
    out = odr.run()
    
    c = out.beta[0] # slope of dist atten
    '''          
    fig = plt.figure(1, figsize=(18,5))
    
    plt.plot(log10(rhyps), yres, '+', c='0.5')
    
    # identify outliers here!!!
    oidx = where(array(yres) > 1.25)[0]
    for oi in oidx:
        plt.text(log10(rhyps[oi]), yres[oi], ' '.join((evdt[oi], stas[oi])), fontsize=5)
    oidx = where(array(yres) < -1.25)[0]
    for oi in oidx:
        plt.text(log10(rhyps[oi]), yres[oi], ' '.join((evdt[oi], stas[oi])), fontsize=5)
    
    plt.plot(log10(array([5, 2000])), [0, 0], 'k--')
    plt.xlim(log10(array([5, 2000])))
    plt.ylim([-3, 3])
    
    plt.savefig('annotated_residuals.png', fmt='png', dpi=300, bbox_inches='tight')
    plt.show()
    '''
    """
plt.show()
###############################################################################
# write params
###############################################################################

pklfile = open('atten_coeffs.pkl', 'wb')
pickle.dump(coeffs, pklfile, protocol=-1)
pklfile.close()

'''
txt = 'ln Y = m0 + (m1-4) * M^2 +(m2-4) * M + r1 * log10(sqrt(rhyp^2 + r2^2)) | rhyp <= '+str(r1)+'\n'
txt += 'ln Y = m0 + (m1-4) * M^2 +(m2-4) * M + r1 * log10(sqrt(rref^2 + r2^2)) + r4 * log10(rhyp / rref) + r5 * (rhyp - rref) | rhyp > '+str(r1)+'\n'
txt += 'm0' + '\t' + str(m0) + '\n'
txt += 'm1' + '\t' + str(m1) + '\n'
txt += 'm2' + '\t' + str(m2) + '\n'
txt += 'r1' + '\t' + str(nc[0]) + '\n'
txt += 'r2' + '\t' + str(nc[2]) + '\n'
txt += 'r3' + '\t' + str(nc[1]) + '\n'
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
