import pickle
from numpy import unique, array, arange, log, log10, exp, mean, nanmean, ndarray, std, \
                  nanmedian, nanstd, vstack, pi, nan, isnan, interp, where, zeros_like, ones_like
from scipy.stats import linregress, trim_mean
from misc_tools import get_binned_stats, dictlist2array, get_mpl2_colourlist
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
# correct Wei et al (2017) attenuation
####################################################################################
def build_Q_data():
    from numpy import array, loadtxt
    
    # read parameter file
    wei_freqs = []
    filenums = []
    lines = open('wei_etal_2017_lgq/parameters.txt').readlines()[1:]
    for line in lines:
        dat = line.split()
        if float(dat[1]) >= 0.1:
            filenums.append(dat[0].zfill(2))
            wei_freqs.append(float(dat[1]))
    
    # now load all data
    qDatList = []
    for freq, filenum in zip(wei_freqs, filenums):
        file_name = "wei_etal_2017_lgq/Q." + filenum +  "..txt"
        data = loadtxt(file_name)
        qDatList.append(data)
        
    return qDatList, wei_freqs

def get_centre_ray(elon, elat, slon, slat):
    """ locates the centre point of a ray.  An alternitive to this
    would be to locate the deepest point of the ray, but this is
    unessecery
    """
    from mapping_tools import distance, reckon
    
    rngkm, az, baz = distance(elat, elon, slat, slon)
    clon, clat = reckon(elat, elon, rngkm/2., az)
    
    return clon, clat

def get_all_Qs_centre(clat, clon, qDatList, half_bin_size=2.5):
    """
    OPTION 2:
    function draws a default 6 degree polygon around the point at the centre 
    of a ray and randomly samples points within the polygon.  
    Sized based on the spacing in the Q model.  Also returns a standard deviation
    Q value calculated in log space.  

    This option chosen over the ray as there is less variation in the
    values, and the code is much quicker to run (smaller nested loop).  
    """
    from numpy import array, mean, std, where

    minlon = clon - half_bin_size
    maxlon = clon + half_bin_size
    minlat = clat - half_bin_size
    maxlat = clat + half_bin_size
    
    Q_list = []
    stdevQ_list = []
    
    #print('looping freqs')
    for qdat in qDatList:
        lons = qdat[:,0]
        lats = qdat[:,1]
    
        # get q vals within region
        idx = where((lons >= minlon) & (lons <= maxlon) \
                    & (lats >= minlat) & (lats <= maxlat))[0] 
        
        Q_vals = qdat[:,2][idx]
        	
        # get mean & stdev
        Q_list.append(mean(Q_vals))
        stdevQ_list.append(std(Q_vals))
        
    return array(Q_list), array(stdevQ_list)
    
    
def correct_wei_atten(rec, qDatList, wei_freqs):
    channel = rec['channels'][0]
    raw_fds = rec[channel]['swave_spec']
    sn_ratio = rec[channel]['sn_ratio']
    freqs = rec[channel]['freqs']
    
    sn_thresh = 4
    
    # get Q(f)
    clon, clat = get_centre_ray(rec['eqlo'], rec['eqla'], rec['stlo'], rec['stla'])
    Q_list, stdevQ_list = get_all_Qs_centre(clat, clon, qDatList)
    wei_freqs_log = log10(wei_freqs)
    Q_list_log = log10(Q_list)
    
    # regress values and fit bilinear
    log_Qfact = bilinear_Q_regression(wei_freqs_log, Q_list_log)
    #print(log_Qfact)
    
    # interpolate to freqs
    Qinterp = 10**(interp(log10(freqs), wei_freqs_log, log_Qfact))
    
    #print(Qinterp)
    Qterm = exp(-1 * pi * freqs * rec['rhyp'] / (vs * Qinterp))
    #print(log10(Qterm))
    
    cor_fds = raw_fds / Qterm
    
    # convert to acceleration
    cor_fas = cor_fds * (2.*pi*freqs)**2
    
    # now set spectra < sn_thresh to nan
    idx = where(sn_ratio < sn_thresh)[0]
    cor_fas_nan = cor_fas.copy()
    cor_fas_nan[idx] = nan
            
    return cor_fas, cor_fas_nan, freqs, log_Qfact

def regress_rec_kappa(rec):
    from numpy import log, isnan, nan, pi
    from scipy.stats import linregress
    
    channel = rec['channels'][0]
    cor_fds = rec[channel]['swave_spec']
    sn_ratio = rec[channel]['sn_ratio']
    sample_rate = rec[channel]['sample_rate']
    freqs = rec[channel]['freqs']
    
    sn_thresh = 5
    
    cor_fas = cor_fds * (2.*pi*freqs)**2
    
    # now set spectra < sn_thresh to nan
    maxfreq = min([0.425*sample_rate, 10.])
    #print(maxfreq)
    idx = where((sn_ratio < sn_thresh) | (freqs >= maxfreq))[0]
    cor_fas_nan = cor_fas.copy()
    cor_fas_nan[idx] = nan
    
    # regress record kappa
    ridx = where((freqs >= 4.) & (isnan(cor_fas_nan) == False))[0] # assume 5Hz past corner
    if len(ridx) <= 5: # need at least 5 data points
        lamda = nan
        kappa = nan
    else:
        lamda = linregress(freqs[ridx], log(cor_fas_nan[ridx]))
        #print(lamda)
        kappa = -lamda[0] / pi
            
    return kappa, rec['rhyp']

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

'''
def get_freq_Q_log(stlo, stla, eqlo, eqla, evdp):
    """ run functions to get Q and freq as lists"""

    freqs = extract_Q_freq.get_freq_list()
    freqs_log = get_average_Q_list.convert_list_to_log(freqs)
    #print(freqs_log)
    Q_list, stdevQ_list = get_average_Q_list.get_all_Qs_centre(eqla, eqlo, stla, stlo, evdp)
    print(Q_list)
    Q_list_log = get_average_Q_list.convert_list_to_log(Q_list)
    print(Q_list_log)
    
    return freqs_log, Q_list_log, stdevQ_list
'''

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

####################################################################################
# start main
####################################################################################
fig = plt.figure(1, figsize=(18,11))

# loop through & plot station  data
cs = get_mpl2_colourlist()

handles1 = []
labels1 = []

events = unique(dictlist2array(recs, 'ev'))
stations = unique(dictlist2array(recs, 'sta'))

# build Wei et al Q data to avoid reading for every record
qDatList, wei_freqs = build_Q_data()

stas = []
kappa_zeros = []
kappa_slopes = []

i = 1
ii = 1
kap_txt = 'STA,LTS_K0, LTS_KR\n'
for sta in stations:
    cnt = 0
    ax = plt.subplot(3,3,i)
    kappas = []
    rhyps = []

    # get all records for each sta
    for rec in recs:
        if rec['sta'] == sta and len(rec['channels']) > 0:
            
            # get corrected fas
            cor_fas, cor_fas_nan, freqs, log_Qfact = correct_wei_atten(rec, qDatList, wei_freqs)	
            
            # get kappa & dist
            rec_kappa, rhyp = regress_rec_kappa(rec)
            
            kappas.append(rec_kappa)
            rhyps.append(rhyp)
            
            print(rec['ev'] + ' ' + rec['sta'] + ' ' + str('%0.5f' % rec_kappa))
            
    # plot fas
    plt.plot(rhyps, kappas, 'o', c='0.6', lw=0.4)
    
    # now regress for kappa zero
    kappas = array(kappas)
    rhyps = array(rhyps)
    nidx = where(isnan(kappas) == False)[0]
    
    lts_slope = nan
    lts_inter = nan
    lts_rval  = nan

    if len(nidx) > 0:
        kappa_regress = linregress(rhyps[nidx], kappas[nidx])
        kappa_zero = kappa_regress[1]
        kappa_slope = kappa_regress[0]
        
        # do least trimmed squares
        '''
        sigx = ones_like(rhyps[nidx]) * 10.
        sigy = ones_like(rhyps[nidx]) * 0.01
        lts_fit = lts_linefit.lts_linefit(rhyps[nidx], kappas[nidx], sigx, sigy)
        lts_inter, lts_slope = lts_fit.ab
        '''
        
        # do trimmed lin reg
        mod_pred = kappa_regress[1] + kappa_regress[0] * rhyps
        mod_res = kappas - mod_pred
        mod_std = nanstd(mod_res)
        
        # remove vals outside 2 std
        idx = where((abs(mod_res) <= 2*mod_std) & (isnan(kappas)== False))[0]
        plt.plot(rhyps[idx], kappas[idx], 'ro', lw=0.4)
        
        # re-regress
        if len(idx) > 0:
           lts_fit = linregress(rhyps[idx], kappas[idx])
           lts_slope = lts_fit[0]
           lts_inter = lts_fit[1]
           lts_rval = lts_fit[2]
           
        else:
           lts_slope = nan
           lts_inter = nan
           lts_rval  = nan
        
        # do some checks 
        if kappa_zero < -0.15 or kappa_slope < 0 or len(kappas[nidx]) < 3:
           # take mean
           kappa_zero = trim_mean(kappas[nidx], 0.125)
           kappa_slope = 0.
           
        # do some more checks
        if lts_inter < -0.15 or lts_slope < 0 or len(kappas[idx]) < 3:
           # take mean
           lts_inter = trim_mean(kappas[nidx], 0.125)
           lts_slope = 0.
           
    else:
        kappa_zero = nan
        kappa_slope = nan
    
    kap_txt += ','.join((sta, str(lts_inter), str(lts_slope))) + '\n'
    # plt kappa_r
    xplt = array([0, 2500])
    yplt = kappa_zero + kappa_slope*xplt
    plt.plot(xplt, yplt, 'r-', lw=2)
    
    yplt = lts_inter + lts_slope*xplt
    plt.plot(xplt, yplt, 'g-', lw=2)
        
    plt.title(sta)
    plt.xlim([0, 2500])
    plt.ylim([-0.15, 0.3])
    # prep for next subplot
    i += 1
    if i > 9:
        plt.savefig('kappa/site_kappa_dist_'+str(ii)+'.png', fmt='png', bbox_inches='tight')
        i = 1
        ii += 1
        fig = plt.figure(ii, figsize=(18,11))
    
    print(sta+' '+str(cnt))
    
plt.savefig('kappa/site_kappa_dist_'+str(ii)+'.png', fmt='png', bbox_inches='tight')

# write txt
f = open('site_kappas.csv', 'w')
f.write(kap_txt)
f.close()