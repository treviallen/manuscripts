import pickle
from numpy import unique, array, arange, log, log10, logspace, exp, mean, nanmean, ndarray, \
                  nanmedian, hstack, pi, nan, isnan, interp, where, zeros_like, polyfit, sqrt
from misc_tools import get_binned_stats, dictlist2array, get_mpl2_colourlist
#from js_codes import get_average_Q_list, extract_Q_freq
#from mapping_tools import distance
from scipy.odr import Data, Model, ODR, models
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.stats import linregress
import scipy.odr.odrpack as odrpack
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
#fig = plt.figure(1, figsize=(18,11))

events = unique(dictlist2array(recs, 'ev'))
mags = dictlist2array(recs, 'mag')
rhyp = dictlist2array(recs, 'rhyp')
stations = unique(dictlist2array(recs, 'sta'))

mrng = arange(3.8, 6.9, 0.1)
minDist = 10**0.5
maxDist = 2200
minRegDist = 100
maxRegDist = 1000

bins = arange(log10(minDist), log10(maxDist), 0.1)
    
log_norm_amps = []

# set freq index
# idx 35 = 0.74989421 Hz
fidx = 35
#hardwire nf corner distance
#fidx = 40 # 1 Hz
#fidx = 56 # 2.0 Hz
#fidx = 67 # 2.0 Hz
#fidx = 20

minr = 5.
maxr = 120 # max dist for near source
nref = 5.

chan = recs[0]['channels'][0]
freq = recs[0][chan]['freqs'][fidx]
print("Reg Freq = " +str('%0.3f' % freq))

reg_freqs = logspace(-1,1,31)

i = 1
for m in mrng:
    cnt = 0
    #ax = plt.subplot(4,5,i)
    
    mrhyps = []
    mamps  = []

    # get all records for each sta
    for rec in recs:
        if len(rec['channels']) > 0 and rec['mag'] >= m-0.05 and rec['mag'] < m+0.05:
            if rec['net'] in keep_nets:
                channel = rec['channels'][0]
                
                if rec[channel]['sn_ratio'][fidx] >= 4.:
                    mrhyps.append(rec['rhyp'])
                    mamps.append(rec[channel]['swave_spec'][fidx])
    
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
        
        # normalise data @ 630 km
        #nidx = where((binstrp > 2.79) & (binstrp < 2.81))[0]
        # normalise data @ 500 km
        nidx = where((binstrp > 2.69) & (binstrp < 2.71))[0]
        
        if len(nidx) > 0:
            print (m, nidx)
            namps = log10(mamps) - logmedamp[nidx]
            
            if len(log_norm_amps) == 0:
                log_norm_amps = namps
                norm_rhyps = mrhyps
            else:
                log_norm_amps = hstack((log_norm_amps, namps))
                norm_rhyps = hstack((norm_rhyps, mrhyps))
        '''    
        print(m)
        print(10**binstrp)
        '''
        # regress mag-dependent data and get intercept
        didx = where((medx >= 2) & (medx <= 3.35))[0]
        if len(didx) > 1:
            gr_fit = linregress(medx[didx], logmedamp[didx])
            
            xplt = array([2, 3])
            yplt = gr_fit[0] * xplt + gr_fit[1]
            
            #plt.loglog(10**xplt, 10**yplt, 'k-', lw=2)
            
#plt.savefig('mag_specific_geom_spread.png', fmt='png', bbox_inches='tight')
#plt.show()

###############################################################################
# plt near-field GR
###############################################################################

fig = plt.figure(2, figsize=(9,7))
#print('2016-05-20T18.10.AU.WRKA?') 

plt.loglog(norm_rhyps, 10**log_norm_amps, '+', c='0.6', lw=0.5, ms=6)
plt.xlim([5, 2250])
plt.ylim([5E-3, 500])
plt.ylabel('Normalised Spectral Amplitude (0.75 Hz)')
plt.xlabel('Hypocentral Distance (km)')

# get binned data
logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log10(norm_rhyps), log_norm_amps)
plt.loglog(10**medx, 10**logmedamp, 'rs', ms=7)

didx = where((medx >= log10(minr)) & (medx <= log10(maxr)) & (nperbin > 2))[0]
gr_fit = linregress(medx[didx], logmedamp[didx])

# get near-source fit
def fit_near_source_saturation(c, x):
    from numpy import sqrt
    
    D = sqrt((10**x)**2 + c[2]**2) # were c1 ~= 5-10?
    #D = sqrt((10**x)**2 + nref**2)
    
    ans = c[0] * log10(D) + c[1]
    
    return ans

# fit all data
data = odrpack.RealData(medx[didx], logmedamp[didx])

afit = odrpack.Model(fit_near_source_saturation)
odr = odrpack.ODR(data, afit, beta0=[-1., 2, 1.])
#odr = odrpack.ODR(data, afit, beta0=[-1., 2])

odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
out = odr.run()
nc = out.beta

# plot
#xrng = arange(0, log10(maxr)+0.02, 0.02)
xrng_nf = log10(arange(1, maxr+1))

D = sqrt((10**xrng_nf)**2 + nc[2]**2)
#D = sqrt((10**xrng_nf)**2 + nref**2)
yrng_nf = nc[0] * log10(D) + nc[1]

plt.loglog(10**xrng_nf, 10**yrng_nf, 'k-', lw=2)

###############################################################################
# plt far-field GR
###############################################################################

# get near-source fit
def fit_far_field(c, x):
    from numpy import sqrt, log10
    
    D = sqrt((10**x)**2 + nc[2]**2) # were c1 ~= 5-10?
    #D = sqrt((10**x)**2 + nref**2)
    ans = nc[0] * log10(D) + nc[1]
    
    idx = 10**x >= maxr
    D = sqrt(maxr**2 + nc[2]**2) # were c1 ~= 5-10?
    #D = sqrt(maxr**2 + nref**2)
    ans[idx] = nc[0] * log10(D) + nc[1] + c[0] * log10(10**x[idx] / maxr) + c[1] * (10**x[idx] - maxr)
    
    return ans
    
# fit all data
didx = where((medx >= log10(minr)) & (nperbin > 2))[0] #log10(maxr))[0]
print(medx, didx)
data = odrpack.RealData(medx[didx], logmedamp[didx])

afit = odrpack.Model(fit_far_field)
odr = odrpack.ODR(data, afit, beta0=[-1., -0.001])

odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
out = odr.run()
fc = out.beta

# plot
xrng_ff = arange(log10(maxr), log10(2100), 0.02)
D = sqrt(maxr**2 + nc[2]**2)
#D = sqrt(maxr**2 + nref**2)
yrng_ff = nc[0] * log10(D) + nc[1] + fc[0] * log10(10**xrng_ff / maxr) + fc[1] * (10**xrng_ff - maxr)

plt.loglog(10**xrng_ff, 10**yrng_ff, 'g-', lw=2)
    
plt.savefig('norm_geom_spread.png', fmt='png', bbox_inches='tight')
plt.show()

###############################################################################
# fit mag intercept
###############################################################################

def fit_mag_intercept(c, x):
    from numpy import sqrt, log10
    
    D = sqrt(x**2 + nc[2]**2) # were c1 ~= 5-10?
    #D = sqrt((10**x)**2 + nref**2)
    
    ans = nc[0] * log10(D) + c[0]
    
    idx = x >= maxr
    D = sqrt(maxr**2 + nc[2]**2) # were c1 ~= 5-10?
    #D = sqrt(maxr**2 + nref**2)
    
    ans[idx] = nc[0] * log10(D) + c[0] + fc[0] * log10(x[idx] / maxr) + fc[1] * (x[idx] - maxr)
    
    return ans

###############################################################################
# build mag regression data
###############################################################################

            
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
                
                if rec[channel]['sn_ratio'][fidx] >= 4.:
                    mrhyps.append(rec['rhyp'])
                    mamps.append(rec[channel]['swave_spec'][fidx])
    
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
        didx = where((medx >= 0.1) & (nperbin > 2))[0] #log10(maxr))[0]

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

###############################################################################
# get mag scaling
###############################################################################

fig = plt.figure(4, figsize=(8, 8))

plt.plot(m_reg, m_intercept, 'rs')

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

plt.savefig('mag_scaling.png', fmt='png', bbox_inches='tight')
plt.show()

###############################################################################
# get residuals and far-field correction
###############################################################################
#y = m0 + m1 * (mags-4.)**2. + m2 * (mags-4.) + r1 * log10(rhyp)

yres = [] 
rhyps = []
stas = []
evdt = []
for rec in recs:
    if rec['net'] in keep_nets:
        try:
            channel = rec['channels'][0]
                
            if rec[channel]['sn_ratio'][fidx] >= 4.:
                rhyps.append(rec['rhyp'])
                
                #!!! use new atten coeffs !!!
                if rec['rhyp'] <= maxr:
                    D = sqrt(rec['rhyp']**2 + nc[2]**2)
                    y = m0 + m1 * (rec['mag']-4.)**2. + m2 * (rec['mag']-4.) \
                        + nc[0] * log10(D)
                else:
                    D = sqrt(maxr**2 + nc[2]**2)
                    y = m0 + m1 * (rec['mag']-4.)**2. + m2 * (rec['mag']-4.) \
                        + nc[0] * log10(D) + fc[0] * log10(rec['rhyp'] / maxr) + fc[1] * (rec['rhyp'] - maxr)
                
                yobs = log10(rec[channel]['swave_spec'][fidx])
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
plt.plot(medx, medbin, 'rs', ms=6.5)

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

###############################################################################
# write scaling params
###############################################################################
txt = 'ln Y = m0 + (m1-4) * M^2 +(m2-4) * M + r1 * log10(sqrt(rhyp^2 + r2^2)) | rhyp <= '+str(maxr)+'\n'
txt += 'ln Y = m0 + (m1-4) * M^2 +(m2-4) * M + r1 * log10(sqrt(rref^2 + r2^2)) + r4 * log10(rhyp / rref) + r5 * (rhyp - rref) | rhyp > '+str(maxr)+'\n'
txt += 'm0' + '\t' + str(m0) + '\n'
txt += 'm1' + '\t' + str(m1) + '\n'
txt += 'm2' + '\t' + str(m2) + '\n'
txt += 'r1' + '\t' + str(nc[0]) + '\n'
txt += 'r2' + '\t' + str(nc[2]) + '\n'
txt += 'r3' + '\t' + str(nc[1]) + '\n'
txt += 'r4' + '\t' + str(fc[0]) + '\n'
txt += 'r5' + '\t' + str(fc[1]) + '\n'
txt += 'rref' + '\t' + str(maxr) + '\n'
txt += 'f' + '\t' + str('%0.3f' % freq) + '\n'
txt += 'fidx' + '\t' + str(fidx)

#print('!!!!!!!! REMEMEBR TO UPDATE OUTPUT !!!!!!!!')
f = open('basic_atten_coeffs_tmp.txt', 'w')
f.write(txt)
f.close()

