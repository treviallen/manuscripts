import pickle
from numpy import unique, array, arange, log, log10, logspace, exp, mean, nanmean, ndarray, \
                  nanmedian, hstack, pi, nan, isnan, interp, where, zeros_like, polyfit
from misc_tools import get_binned_stats, dictlist2array, get_mpl2_colourlist
#from js_codes import get_average_Q_list, extract_Q_freq
#from mapping_tools import distance
#from scipy.odr import Data, Model, ODR, models
#import scipy.odr.odrpack as odrpack
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

####################################################################################
# start main
####################################################################################
#fig = plt.figure(1, figsize=(18,11))

events = unique(dictlist2array(recs, 'ev'))
mags = dictlist2array(recs, 'mag')
rhyp = dictlist2array(recs, 'rhyp')
stations = unique(dictlist2array(recs, 'sta'))

mrng = arange(4.4, 6.9, 0.1)
minDist = 10
maxDist = 2200
minRegDist = 100
maxRegDist = 1000

bins = arange(log10(minDist), log10(maxDist), 0.1)                                                                              
    
log_norm_amps = []

# set freq index
# idx 35 = 0.74989421 Hz
fidx = 35
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
        nidx = where((binstrp > 2.79) & (binstrp < 2.81))[0]
        
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
# plt normalised GR
###############################################################################

fig = plt.figure(2, figsize=(9,7))
#print('2016-05-20T18.10.AU.WRKA?') 

plt.loglog(norm_rhyps, 10**log_norm_amps, '+', c='0.6', lw=0.5, ms=6)
plt.xlim([10, 2250])
plt.ylim([5E-3, 500])
plt.ylabel('Normalised Spectral Amplitude (0.75 Hz)')
plt.xlabel('Hypocentral Distance (km)')

# get binned data
logmedamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, log10(norm_rhyps), log_norm_amps)
plt.loglog(10**medx, 10**logmedamp, 'rs', ms=7)

didx = where((medx >= 2) & (medx <= 3.35))[0]
gr_fit = linregress(medx[didx], logmedamp[didx])

plt.savefig('norm_geom_spread.png', fmt='png', bbox_inches='tight')
plt.show()

###############################################################################
# fix GR and refit data
###############################################################################
meanslope = gr_fit[0]
def linear_fixed_slope(c, x):
        
    return meanslope * x + c[0]

            
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
        
        # normalise data @ 630 km
        nidx = where((binstrp > 2.79) & (binstrp < 2.81))[0]
        print(nidx)
        
        if len(nidx) > 0:
            
            namps = log10(mamps) - logmedamp[nidx]
            
            if len(log_norm_amps) == 0:
                log_norm_amps = namps
                norm_rhyps = mrhyps
            else:
                log_norm_amps = hstack((log_norm_amps, namps))
                norm_rhyps = hstack((norm_rhyps, mrhyps))
        
        # regress mag-dependent data and get intercept
        didx = where((medx >= 2) & (medx <= 3.35))[0]
       
        data = odrpack.RealData(medx[didx], logmedamp[didx])
        intfit = odrpack.Model(linear_fixed_slope)
        odr = odrpack.ODR(data, intfit, beta0=[5.0])
        
        odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
        out = odr.run()
        intercept = out.beta
        
        if round(m,1) != 4.9 and round(m,1) != 5.5: #or round(m,1) > 6.6
            m_intercept.append(intercept[0])
            m_reg.append(m)
        
        xplt = array([2, 3])
        yplt = meanslope * xplt + intercept[0]
        
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
        
plt.plot(mrng, yplt, 'k-', lw=2)

plt.savefig('mag_scaling.png', fmt='png', bbox_inches='tight')
plt.show()

###############################################################################
# get residuals and far-field correction
###############################################################################
r1 = meanslope
#y = m0 + m1 * (mags-4.)**2. + m2 * (mags-4.) + r1 * log10(rhyp)

yres = [] 
rhyps = []
for rec in recs:
    try:
        channel = rec['channels'][0]
            
        if rec[channel]['sn_ratio'][fidx] >= 4.:
            rhyps.append(rec['rhyp'])
            
            y = m0 + m1 * (rec['mag']-4.)**2. + m2 * (rec['mag']-4.) + r1 * log10(rec['rhyp'])
            
            yobs = log10(rec[channel]['swave_spec'][fidx])
            yres.append(yobs - y)
            
        else:
            yres.append(nan)
            rhyps.append(rec['rhyp'])
    
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
    print(xhinge)
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
plt.plot(log10(array([10, 2000])), [0, 0], 'k--')
plt.xlim(log10(array([10, 2000])))
plt.ylim([-3, 3])

# plot correction
xplt = array([xhinge, xmax])
yplt = c * (xplt - xhinge)
plt.plot(xplt, yplt, 'k-', lw=2)
plt.show()

###############################################################################
# write scaling params
###############################################################################
txt = 'ln Y = m0 + (m1-4) * M^2 +(m2-4) * M + r1*log10(rhyp) + r2*(log10(rhyp) - 3)\n'
txt += 'm0' + '\t' + str(m0) + '\n'
txt += 'm1' + '\t' + str(m1) + '\n'
txt += 'm2' + '\t' + str(m2) + '\n'
txt += 'r1' + '\t' + str(meanslope) + '\n'
txt += 'r2' + '\t' + str(c) + '\n'

print('!!!!!!!! REMEMEBR TO UPDATE OUTPUT !!!!!!!!')
f = open('basic_atten_coeffs_tmp.txt', 'w')
f.write(txt)
f.close()

