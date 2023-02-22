import pickle
from numpy import unique, array, arange, log, log10, logspace, exp, mean, nanmean, ndarray, \
                  nanmedian, hstack, pi, nan, isnan, interp, polyfit, where, zeros_like, polyfit, sqrt
from misc_tools import get_binned_stats, dictlist2array, get_mpl2_colourlist, savitzky_golay
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

"""
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
    #sets up bilinear equation with a free hinge position
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
"""    
# get near-source fit
nref = 1.
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
    
    #print('!!!!!!! UNCOMMENT FILTER BY SAMPLE RATE !!!!!!!')
    
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
                            
                        if rec[channel]['sn_ratio'][fidx[p]] >= 5. and addData == True:
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
    odr = odrpack.ODR(data, afit, beta0=[-1.5, 2])
    
    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
    out = odr.run()
    nc = out.beta
    
    return nc
    
def refit_near_source_atten(plt, norm_rhyps, log_norm_amps, r1, n0):
    
    # get binned data
    log_norm_rhyps = log10(norm_rhyps)
    logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log_norm_rhyps, log_norm_amps)
    if pltTrue == True:
        plt.loglog(10**medx, 10**logmedamp, 'rs', ms=7)
    
    def refit_near_source_saturation(c, x):
        from numpy import sqrt
        
        D = sqrt((10**x)**2 + nref**2) # were c1 ~= 5-10?
        
        ans = n0 * log10(D) + c[0]
        
        return ans
        
    # fit all data - keep me!
    didx = where((medx >= log10(minr)) & (medx < log10(r1)) & (nperbin > 2))[0] 
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

#nref2 = 180.    
def fit_mid_field_atten(plt, norm_rhyps, log_norm_amps, n0, n1):
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
    yrng_ff = n0 * log10(D1) + n1 + mc[0] * log10(10**xrng_ff / r1) + mc[1] * (10**xrng_ff - r1)
        
    if pltTrue == True:
       plt.loglog(10**xrng_ff, 10**yrng_ff, 'g-', lw=2)
    
    return mc
    
def fit_far_field_atten(plt, mc, norm_rhyps, log_norm_amps, n0, n1):
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
fidx=arange(0,150,17)+3 # for testing
fidx=arange(0,len(freqs)-1,1) # for regressing all coeffs

pltTrue = True
if len(fidx) > 12:
    pltTrue = False

if pltTrue == True:
    fig = plt.figure(1, figsize=(18,11))
    
coeffs = []
nc0_array = []

for p, freq in enumerate(freqs[fidx]):
    if pltTrue == True:
        fig = plt.figure(1, figsize=(18,11))
        plt.subplot(3,4,p+1)
        plt.ylabel(str(freq), fontsize=7)
    
    minr = 5.
    r1 = 80 # max dist for near source
    r2 = 500
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

    coeffs.append({'nc0': nc[0], 'r1': r1, 'r2': r2, 'nref':nref, 'freq':freq})
    
    # set slope array for smoothing
    nc0_array.append(nc[0])

###############################################################################
# smooth nc coeff
###############################################################################
if pltTrue == True:
    sg_window = 3
else:
    sg_window = 41
    
sg_poly = 2
smooth_nc0 = savitzky_golay(array(nc0_array[3:]), sg_window, sg_poly) # slope

for i, c in enumerate(coeffs):
    if i < 3:
        coeffs[i]['nc0s'] = coeffs[i]['nc0']
    else:
        coeffs[i]['nc0s'] = smooth_nc0[i-3]

# using smoothed n0 values, refit
for p, freq in enumerate(freqs[fidx]):
    
    # set n0 - use smoothed vals
    n0 = coeffs[p]['nc0s']
    
    # get normalised amplitudes
    log_norm_amps, norm_rhyps, mags, logamps = normalise_data(p, recs)
    
    ###############################################################################
    # plt near-field GR
    ###############################################################################
    
    #plt.ylabel('Normalised Spectral Amplitude')
    #plt.xlabel('Hypocentral Distance (km)')
    
    n1 = refit_near_source_atten(plt, norm_rhyps, log_norm_amps, r1, n0)[0]
    print(n1)
 
    # plot
    xrng_nf = log10(arange(1, r1+1))
    
    #D = sqrt((10**xrng_nf)**2 + nc[2]**2)
    D1 = sqrt((10**xrng_nf)**2 + nref**2)
    yrng_nf = n0 * log10(D1) + n1
    
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
        ans = n0 * log10(D1) + n1
        
        idx = where((10**x > r1) & (10**x <= r2))[0]
        D1 = sqrt(r1**2 + nref**2)
        
        ans[idx] = n0 * log10(D1) + n1 + c[0] * log10(10**x[idx] / r1) + c[1] * (10**x[idx] - r1)
            
        return ans
    
    mc = fit_mid_field_atten(plt, norm_rhyps, log_norm_amps, n0, n1)
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
        ans = n0 * log10(D1) + n1
        
        # set mid-field
        idx = where((10**x > r1) & (10**x <= r2))[0]
        D1 = sqrt(r1**2 + nref**2)
        ans[idx] = n0 * log10(D1) + n1 + mc[0] * log10(10**x[idx] / r1) + mc[1] * (10**x[idx] - r1)
        
        idx = where(10**x > r2)[0]
        ans[idx] = n0 * log10(D1) + n1 \
                   + mc[0] * log10(r2 / r1) + mc[1] * (r2 - r1) \
                   + c[0] * log10(10**x[idx] / r2) + c[1] * (10**x[idx] - r2)
            
        return ans
    
    fc = fit_far_field_atten(plt, mc, norm_rhyps, log_norm_amps, n0, n1)
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
    
    coeffs[p]['nc1s'] = n1
    coeffs[p]['mc0'] = mc[0]
    coeffs[p]['mc1'] = mc[1]
    coeffs[p]['fc0'] = fc[0]
    coeffs[p]['fc1'] = fc[1]
    coeffs[p]['magc0'] = magc[0]
    coeffs[p]['magc1'] = magc[1]
    


if pltTrue == True:
    plt.show()
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
