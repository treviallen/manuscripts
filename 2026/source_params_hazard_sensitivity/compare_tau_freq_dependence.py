import pickle
from numpy import unique, array, arange, log, log10, exp, mean, nanmean, ndarray, sqrt, \
                  nanmedian, hstack, pi, nan, isnan, interp, where, zeros_like, polyfit, \
                  hstack, nanstd, loadtxt, nanmedian, nanvar
from misc_tools import get_binned_stats, dictlist2array, get_mpl2_colourlist, savitzky_golay, get_log_xy_locs
from mag_tools import nsha18_mb2mw, nsha18_ml2mw
from get_mag_dist_terms import get_distance_term, get_magnitude_term, get_kappa_term, get_regional_term
from scipy.stats import linregress, median_abs_deviation, f
import scipy.odr.odrpack as odrpack
from obspy import UTCDateTime
from datetime import timedelta
from response import stationlist2dict
import warnings
warnings.filterwarnings("ignore")
import shapefile
from sys import argv
from shapely.geometry import Point, Polygon
from mapping_tools import get_field_data

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import colorbar
import matplotlib as mpl
from gmt_tools import cpt2colormap
mpl.style.use('classic')
plt.rcParams['pdf.fonttype'] = 42

cor_dist = 100
def fit_regional_correction(medx, logmedamp):
    idx = where(10**medx >= cor_dist)[0]

    data = odrpack.RealData(medx[idx], logmedamp[idx])

    # fit all as free params
    afit = odrpack.Model(correct_far_field)
    odr = odrpack.ODR(data, afit, beta0=[0.0])

    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
    out = odr.run()
    ffc = out.beta

    return ffc[0]
    
keep_nets = set(['AU', 'IU', 'S1', 'II', 'G', 'MEL', 'ME', '2O', 'AD', 'SR', 'UM', 'AB', 'VI', 'GM', 'M8', 'DU', 'WG', '4N', \
                 '1P', '1P', '2P', '6F', '7K', '7G', 'G', '7B', '4N', '7D', '', 'OZ', 'OA', 'WG', 'XX', 'AM', 'YW', '3B', '1K', \
                 '1Q', '3O', '7F', '6K', '5G', '5C'])

ignore_stas = open('../../2023/au_stress_drop/sta_ignore.txt').readlines()
#ignore_stas = open('sta_ignore.test').readlines()
ignore_stas = set([x.strip() for x in ignore_stas])

###############################################################################
# load datasets 
###############################################################################

recs = pickle.load(open('../../2023/au_stress_drop/fft_data.pkl', 'rb' ))
#recs = pickle.load(open('fft_data_mag_match.pkl', 'rb' ))
chan = recs[0]['channels'][0]
freqs = recs[0][chan]['freqs']

# load final Mw estimates
lines = open('brune_stats_cluster.csv').readlines()[1:]
brunedat = []
for line in lines:
    dat = line.strip().split(',')
    tmp = {'datetime':UTCDateTime(dat[0]), 'mw':float(dat[8]), 'qual':int(float(dat[-5])), \
           'lon':float(dat[2]), 'lat':float(dat[3]), 'sd':float(dat[10]), 'cluster':int(float(dat[-1]))}
    
    brunedat.append(tmp)
    
def get_brune_deets(fft_datetime):
    bruneStats = {'mw':nan, 'qual':0, 'cluster':nan, 'datetime':UTCDateTime(2599,12,31), 'lon':nan, 'sd':nan}
    for evnum, ev in enumerate(brunedat): 
        #ev['datetime'] = UTCDateTime(2024,2,27,16,4,9)
        
        if fft_datetime > UTCDateTime(ev['datetime']-timedelta(seconds=601)) \
           and fft_datetime < UTCDateTime(ev['datetime']+timedelta(seconds=300)):
            
               bruneStats = ev
    if isnan(bruneStats['mw']):
        print('Event not found: ',fft_datetime)
                       
    return bruneStats

# load atten coeffs
print('run compare_tau_freq_dependence.py ../../2023/au_stress_drop/atten_coeffs_1.3_5km.pkl')
coeff_pkl = argv[1]
coeffs = pickle.load(open(coeff_pkl, 'rb' ))
 
###############################################################################
# do ml-mw conversions
###############################################################################
 
s0, s1, s2 = loadtxt('../../2023/au_stress_drop/mw-ml_coeffs_2800.csv', delimiter=',', skiprows=1)   

lines = open('ml_mw_bias_cluster.csv').readlines()
cmean = []
for line in lines[1:]:
    dat = line.strip().split(' +- ')[0].split(',')
    cmean.append(float(dat[1]))
    
print(cmean)

'''
cdata = loadtxt('ml_mw_bias_cluster.csv', delimiter=',', skiprows=1) 
cmean = cdata[:,1]
'''

def get_ml2mw_cluster(ml, cluster=0):
    #print(ml)
    return s0 * ml**2 + s1 * ml + s2 - cmean[cluster], cmean[cluster] 


# load ml-mw data from "get_cluster_ml_mw_bias.py"
csvfile = '../../2023/au_stress_drop/ml_mw_stats.csv'

print(csvfile)
lines = open(csvfile).readlines()
ml = []
mldate = []

for line in lines[1:]:
    dat = line.strip().split(',')
    mldate.append(dat[0])
    ml.append(float(dat[2]))
    
###############################################################################
# set datasets
###############################################################################
'''
events = unique(dictlist2array(recs, 'ev'))
mags = dictlist2array(recs, 'mag')
stations = unique(dictlist2array(recs, 'sta'))
eqlo = dictlist2array(recs, 'eqlo')
eqla = dictlist2array(recs, 'eqla')
stlo = dictlist2array(recs, 'stlo')
stla = dictlist2array(recs, 'stla')
'''

### !!! COMMENT OUT AFTER FIRST RUN !!!!
'''
# set brune mags
print('Setting brune mags and MLs ...')
for i, rec in enumerate(recs):
    bruneStats = get_brune_deets(UTCDateTime(rec['evdt']))
    
    recs[i]['mwb'] = bruneStats['mw']
    recs[i]['qual'] = bruneStats['qual']
    recs[i]['sd'] = bruneStats['sd']
    recs[i]['mag_cluster'] = bruneStats['cluster']
    
    # lmatch events
    ml_match = nan
    for k, mld in enumerate(mldate):
        if UTCDateTime(rec['evdt']) == mld:
            ml_match = ml[k]
            #print(ml_match)
    recs[i]['ml'] = ml_match
'''
###############################################################################
# loop thru freqs
###############################################################################
fidxs = arange(48, 135,3)
#fidxs = [69,99]

def get_inter_event_terms(magType, recs, maxdist):
    '''
    magType: 1 = Brune mags
             2 = ML-MW cluster correction (MW - dMW) where dMW from mean MW - MW(conv)
             3 = Use converted magnitude
    '''
    resDict = []
    for f, c in enumerate(coeffs):
    
        if f in set(fidxs):
            print("Coeffs Freq = " +str('%0.3f' % c['freq']))
                
            chan = recs[0]['channels'][0]
            freq = recs[0][chan]['freqs'][f]
            print("Reg Freq = " +str('%0.3f' % freq))
            
            if not freq == c['freq']:
               print('\n!!!! Frequency of coefficients inconsistent with data index !!!!\n')
               crash
            
            ###############################################################################
            # parse coefs and get model prediction
            ###############################################################################
            rhyps = []
            yres = []
            rmags = []
            revent = []
            rsd = []
            revdt = []
            for i, rec in enumerate(recs):
                
                # if qual == 1, "mag" is Mwb
                if rec['net'] in keep_nets and rec['qual'] == 1 \
                    and rec['rhyp'] < maxdist \
                    and freq >= rec['ev_fmin'] and freq <= rec['ev_fmax']:
                    
                    if not rec['sta'] in ignore_stas:
                        if len(rec['channels']) > 0:     
                        	  #try:
                            channel = rec['channels'][0]
                    
                            # filter by instrument type
                            addData = True
                            if channel.startswith('SH') or channel.startswith('EH'):
                                if rec[channel]['freqs'][f] < 0.9 and rec['pazfile'].endswith('s6000-2hz.paz'):
                                    addData = False
                                elif rec[channel]['freqs'][f] < 0.4:
                                    addData = False
                            
                            # filer by sample-rate
                            if rec[channel]['freqs'][f] > (0.4 * rec[channel]['sample_rate']):
                                addData = False
                                
                            if rec[channel]['sn_ratio'][f] >= 4. and addData == True:
                                rhyps.append(rec['rhyp'])
                    
                                #magterm = get_magnitude_term(bruneStats['mw'], c)
                                if magType == 1:
                                    prefmw = rec['mag']
                                elif magType == 2:
                                    mwconv, cor = get_ml2mw_cluster(rec['ml2800'], cluster=rec['cluster']) # only use to get correction
                                    #print(rec['cluster'], cor)
                                    prefmw = rec['mag'] - cor 
                                    #print(rec['mag'], cor, prefmw, mwconv, rec['cluster'],rec['evdt'])
                                                                        
                                elif magType == 3:
                                    prefmw, cor = get_ml2mw_cluster(rec['ml2800'], 0) # subtract 1 as we want index 0
                                    
                                    
                                # get mag term
                                magterm = get_magnitude_term(prefmw, c)
                                #print(magterm)
                                
                                # get distance term
                                distterm = get_distance_term(rec['rhyp'], c)
                                
                                #	get distance independent kappa
                                kapterm = get_kappa_term(rec['sta'], c['freq'])
                                
                                #	get regional term
                                regterm = get_regional_term(rec['rhyp'], c, rec['eqdom'])
                                
                                # get total correction
                                ypred = magterm + distterm + kapterm + regterm
                                #print(ypred)
                    
                                yobs = log10(rec[channel]['p-swave_spec'][f])
                                #yobs = log10(rec[channel]['swave_spec'][f])
                                yres.append(yobs - ypred)
                                rmags.append(rec['mag'])
                                revent.append(rec['ev'])
                                revdt.append(rec['evdt'])
                                rsd.append(rec['sd'])
                                #except:
                                #    # do nothing
                                #    dummy = 0
                        
            resData = {'yres':array(yres), 'mags':array(rmags), 'ev':array(revent), \
                       'rhyps':array(rhyps), 'rsd':array(rsd), 'evdt':array(revdt)}
            
            resDict.append(resData)
                    
    return resDict

###############################################################################
# functions for trilinear fit
###############################################################################

def highside(x, hx):
    from numpy import zeros_like
    xmod = zeros_like(x)
    
    idx = x >= hx
    xmod[idx] = 1
    return xmod
    
def midside(x, hx1, hx2):
    from numpy import zeros_like, where
    xmod = zeros_like(x)
    
    idx = where((x >= hx1) & (x < hx2))[0]
    xmod[idx] = 1
    return xmod

def lowside(x, hx):
    from numpy import zeros_like
    xmod = zeros_like(x)
    
    idx = x < hx
    xmod[idx] = 1
    return xmod

def trilinear_reg_fix_intercept_slope(c, x):
    from numpy import zeros_like
    hx1 = c[1] # hinge magnitude
    hx2 = c[2] # hinge magnitude
    #print(hx1, hx2)
    ans3 = zeros_like(x)
    ans2 = zeros_like(x)
    ans1 = zeros_like(x)
    
    modx_lo = lowside(x, hx1)
    modx_md = midside(x, hx1, hx2)
    modx_hi = highside(x, hx2)
    
    ans1 = modx_lo * 0.0
    hy1 = 0.0
    ans2 = modx_md * (c[0] * (x-hx1) + hy1)
    hy2 = c[0] * (hx2-hx1)
    ans3 = modx_hi * hy2
    
    return ans1 + ans2 + ans3

mx1 = 4.75
mx2 = 5.5
    
def trilinear_reg_fix_intercept_slope_mags(c, x):
    from numpy import zeros_like
    hx1 = mx1 # hinge magnitude - average from 0.75-10 Hz
    hx2 = mx2 # hinge magnitude - average from 0.75-10 Hz
    #print(hx1, hx2)
    ans3 = zeros_like(x)
    ans2 = zeros_like(x)
    ans1 = zeros_like(x)
    
    modx_lo = lowside(x, hx1)
    modx_md = midside(x, hx1, hx2)
    modx_hi = highside(x, hx2)
    
    ans1 = modx_lo * 0.0
    hy1 = 0.0
    ans2 = modx_md * (c[0] * (x-hx1) + hy1)
    hy2 = c[0] * (hx2-hx1)
    ans3 = modx_hi * hy2
    
    return ans1 + ans2 + ans3
        
trifitx = arange(3., 7.01, 0.01)

def regress_mbias(plt_event_mags, plt_event_terms_sd_corrected, trifitx):
    bins = arange(3.8, 6.8, 0.2)
    medamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, plt_event_mags, plt_event_terms_sd_corrected)
    realData = odrpack.RealData(medx, medamp)
    
    afit = odrpack.Model(trilinear_reg_fix_intercept_slope_mags)
    odr = odrpack.ODR(realData, afit, beta0=[-0.3])
    
    odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
    out = odr.run()
    c = out.beta
    
    trifity = zeros_like(trifitx)
    idx = where((trifitx >= mx1) & (trifitx < mx2))[0]
    trifity[idx] = c[0] * (trifitx[idx]-mx1)
    hy = c[0] * (mx2-mx1)
    idx = trifitx >= mx2
    trifity[idx] = hy

    return trifity, c
    

        
###############################################################################
# calculate different scenarios
###############################################################################
maxdist = 2000
'''
# ergodic
print('Getting raw residual data ...')
resDict = get_inter_event_terms(1, recs, maxdist)
pklfile = open('residual_data_fdep_'+str(maxdist)+'.pkl', 'wb')
pickle.dump(resDict, pklfile, protocol=-1)
pklfile.close()

#  mw using regional correction
print('Getting cluster mag residual data ...')
resDict = get_inter_event_terms(2, recs, maxdist)
pklfile = open('residual_ml_cluster_data_fdep_'+str(maxdist)+'.pkl', 'wb')
pickle.dump(resDict, pklfile, protocol=-1)
pklfile.close()

# ML2MW no corrections
print('Getting ml2mw residual data ...')
resDict = get_inter_event_terms(3, recs, maxdist)
pklfile = open('residual_ml_nocluster_data_fdep_'+str(maxdist)+'.pkl', 'wb')
pickle.dump(resDict, pklfile, protocol=-1)
pklfile.close()
'''

resDictRaw = pickle.load(open('residual_data_fdep_'+str(maxdist)+'.pkl', 'rb' ))
resDictCluster = pickle.load(open('residual_ml_cluster_data_fdep_'+str(maxdist)+'.pkl', 'rb' ))
resDictNoCluster = pickle.load(open('residual_ml_nocluster_data_fdep_'+str(maxdist)+'.pkl', 'rb' ))

###############################################################################
# get within event sigmas
###############################################################################
def get_inter_event_residuals(resDict):
    
    sigmaDict = []
    sigma_be = []
    sigma_we = []
    mags_list = []
    sds_list = []
    ev_list = []
    ev_term_list = []

    for i, resData in enumerate(resDict):
        
        # get unique events
        events = array(resData['ev'])
        uevents = unique(array(resData['ev']))
        evdt = array(resData['evdt'])
        uevdt = unique(array(resData['evdt']))
        yres = array(resData['yres'])
        mags = array(resData['mags'])
        rhyps = array(resData['rhyps'])
        sds = array(resData['rsd'])
        #print(sds)
        #print(i, len(uevents))
        
        yres_weterm = zeros_like(yres)
        
        # loop thru events
        event_terms = []
        event_mags = []
        event_sds = []
        uevdt = []
        
        for ue in uevents:
            
            #print(ue)
            # get records for given event
            ridx = where(events == ue)[0]
            event_mag = mags[ridx][0]
            
            #print(event_mag)
            event_mags.append(event_mag)
            event_sds.append(sds[ridx][0])
            uevdt.append(evdt[ridx][0])
            
            # get event terms
            event_term = nanmean(log(10**yres[ridx]))
            event_terms.append(event_term)
            
            # remove event terms
            yres_weterm[ridx] = log(10**yres[ridx]) - event_term
            
            '''
            print(event_mag, sds[ridx][0], event_term)
            print(10**yres[ridx])
            print(rhyps[ridx])
            '''
        sigma_be.append(nanstd(array(event_terms)))
        sigma_we.append(nanstd(array(yres_weterm)))   
        
        #sigma_be.append(median_abs_deviation(array(event_terms), nan_policy='omit'))
        #sigma_we.append(median_abs_deviation(array(yres_weterm), nan_policy='omit'))   
        
        mags_list.append(array(event_mags))
        sds_list.append(array(event_sds)) 
        ev_list.append(array(uevdt)) 
        ev_term_list.append(array(event_terms))
                
    return sigma_be, sigma_we, mags_list, sds_list, ev_list, ev_term_list

###############################################################################
# function for f-test
###############################################################################

def f_test(group1, group2):
    import scipy.stats as stats
    
    '''
    in this case, group1 is the ergodic case and group2 is the alternative
    '''
    
    # Calculate the sample variances
    variance1 = nanvar(group1, ddof=1)
    variance2 = nanvar(group2, ddof=1)
    
    # Calculate the F-statistic
    f_value = variance1 / variance2
    
    # Calculate the degrees of freedom
    df1 = len(group1) - 1
    df2 = len(group2) - 1
    
    # two-tailed Calculate the p-value
    #p_value = stats.f.cdf(f_value, df1, df2)
    
    # Upper one-tailed p-value
    p_value = stats.f.sf(f_value, df1, df2)

    
    # Print the results
    '''
    print('Degree of freedom 1:',df1)
    print('Degree of freedom 2:',df2)
    print("F-statistic:", f_value)
    print("p-value:", p_value)
    '''
    return p_value
    
###############################################################################
# get stress drop correction
###############################################################################

plt_freqs = freqs[fidxs]

sigma_be_list, sigma_we_list, mags_list, sds_list, ev_list, ergodic_ev_term_list \
    = get_inter_event_residuals(resDictRaw)

etl =''
for et, sd in zip(ergodic_ev_term_list[1], sds_list[1]):
	etl += ','.join((str(sd), str(et)+'\n'))
	
f = open('etl.csv','w')
f.write(etl)
f.close()

#crash

sdreg_coeffs = []  
sd_corrected_event_terms = []  
sd_corrected_tau = []
clust_sd_corrected_tau = []
ergodic_tau = []
ergodic_ftest_data = []
sd_corrected_ftest = []
clust_sd_corrected_ftest = []
sd_mag_corrected_tau = []
sd_mag_corrected_ftest = []

# get f-dependent dBe vs SD
dBe_sd_r2 = []
mmax = []
mmin = []

for i, f in enumerate(plt_freqs):

    # regress
    idx = where(isnan(ergodic_ev_term_list[i]) == False)[0]
    
    ergodic_tau.append(nanstd(ergodic_ev_term_list[i][idx]))
    ergodic_ev_terms = ergodic_ev_term_list[i][idx]
    #ergodic_tau.append(median_abs_deviation(ev_term_list[i][idx], nan_policy='omit'))
    
    sdreg = linregress(log10(sds_list[i][idx]), ergodic_ev_terms)
    #print(sdreg)
    
    sdreg_coeffs.append(sdreg)
    dBe_sd_r2.append(sdreg[2]**2)
    #print(f, sdreg[2]**2)
    
    # correct event terms for SD
    sd_corrected_event_terms = ergodic_ev_term_list[i][idx] - (sdreg[0] * log10(sds_list[i][idx]) + sdreg[1])
    
    sd_corrected_tau.append(nanstd(ergodic_ev_term_list[i][idx] - (sdreg[0] * log10(sds_list[i][idx]) + sdreg[1])))
    #sdct = median_abs_deviation((ev_term_list[i][idx] - (sdreg[0] * log10(sds_list[i][idx]) + sdreg[1])), \
    #                            nan_policy='omit')
    #sd_corrected_tau.append(sdct)
    
    sd_corrected_ftest.append(f_test(ergodic_ev_terms, sd_corrected_event_terms))
    
    ###############################################################################
    # get mag bias correction
    ###############################################################################
    sd_mag_corrected_event_terms, tc = regress_mbias(mags_list[i][idx], sd_corrected_event_terms, trifitx)
    '''
    if f >= 0.75 and f <= 10:
        mmin.append(tc[1])
        mmax.append(tc[2])
        print(f, tc)
    
    # do some checks
    if tc[0] < -1.0:
        tc[0] = 0
        tc[2] = tc[1]
    elif tc[2] < tc[1]:
        tc[0] = 0
        tc[2] = tc[1]
    '''
    # correct data for mag bias
    fmags_list = mags_list[i][idx]
    sd_mag_corrected_event_terms = sd_corrected_event_terms
    midx = where((fmags_list >= mx1) & (fmags_list < mx2))[0]
    sd_mag_corrected_event_terms[midx] -= tc[0] * (fmags_list[midx]-mx1)
    hy = tc[0] * (mx2-mx1)
    midx = fmags_list >= mx2
    sd_mag_corrected_event_terms[midx] -= hy

    sd_mag_corrected_tau.append(nanstd(sd_mag_corrected_event_terms))
    sd_mag_corrected_ftest.append(f_test(ergodic_ev_terms, sd_mag_corrected_event_terms))
    
    ###############################################################################
    # get cluster SD correction
    ###############################################################################
    
    # get clusters
    clusters = []
    for ev in ev_list[i]:
        bruneStats = get_brune_deets(UTCDateTime(ev))
        clusters.append(bruneStats['cluster'])

    uclusters = unique(clusters)
    #clust_logsd = array(clust_logsd)

    mean_cluster_log_sd = []
    for cluster in uclusters:
        idx = where(clusters == cluster)[0]
        
        mean_cluster_log_sd.append(nanmean(log10(sds_list[i][idx])))

    # now get regional correction    
    ev_term = ergodic_ev_term_list[i].copy()
    sd_clust_cor = zeros_like(ev_term)
    for j, cluster in enumerate(uclusters):
        idx = where(clusters == cluster)[0]
        sd_clust_cor[idx] = ev_term[idx] - (sdreg[0] * mean_cluster_log_sd[j] + sdreg[1])
        
    clust_sd_corrected_tau.append(nanstd(sd_clust_cor))  
    #clust_sd_corrected_tau.append(median_abs_deviation(sd_clust_cor, nan_policy='omit'))
    
    clust_sd_corrected_ftest.append(f_test(ergodic_ev_terms, sd_clust_cor))
    
###############################################################################
# now get alt mag clusterd taus
###############################################################################

sigma_be_list, sigma_we_list, mags_list, sds_list, ev_list, mconv_ev_term_list \
    = get_inter_event_residuals(resDictNoCluster)

mconv_tau = []
mconv_ftest = []
for i, f in enumerate(plt_freqs):

    # regress
    idx = where(isnan(mconv_ev_term_list[i]) == False)[0]
    
    mconv_tau.append(nanstd(mconv_ev_term_list[i][idx])) 
    mconv_ev_terms = mconv_ev_term_list[i][idx]
    #mconv_tau.append(median_abs_deviation(ev_term_list[i][idx], nan_policy='omit'))
    
    mconv_ftest.append(f_test(ergodic_ev_term_list[i], mconv_ev_terms))


sigma_be_list, sigma_we_list, mags_list, sds_list, ev_list, m_clust_ev_term_list \
    = get_inter_event_residuals(resDictCluster)

delta_m_tau = []
delta_m_ftest = []
for i, f in enumerate(plt_freqs):

    # regress
    idx = where(isnan(m_clust_ev_term_list[i]) == False)[0]
    
    delta_m_tau.append(nanstd(m_clust_ev_term_list[i][idx])) 
    m_clust_ev_terms = m_clust_ev_term_list[i][idx]
    #delta_m_tau.append(median_abs_deviation(ev_term_list[i][idx], nan_policy='omit'))
    
    delta_m_ftest.append(f_test(ergodic_ev_term_list[i], m_clust_ev_terms))
    
###############################################################################
# now plot
###############################################################################
ncols = 7
cptfile = '/Users/trev/Documents/DATA/GMT/cpt/keshet.cpt'
cmap, zvals = cpt2colormap(cptfile, ncols)
cs = (cmap(arange(ncols)))

syms = ['o', '^', 's', 'd', 'v', 'h', '<', '>', 'p']

plt.clf()
fig = plt.figure(1, figsize=(8,9))
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)
fig = plt.gcf()
plt.clf()
fig.set_size_inches(8,9)

plt.subplot(211)
plt.cla()
plt.semilogx(plt_freqs, ergodic_tau, syms[0], ls='-', c=cs[0], lw=2, \
             ms=6, mec=cs[0], mfc='w', mew=2, markevery=2, label='Ergodic')
plt.semilogx(plt_freqs, mconv_tau, syms[1], ls='-', c=cs[1], lw=2, \
             ms=6, mec=cs[1], mfc='w', mew=2, markevery=2, label='$\mathregular{M_{conv}}$')
plt.semilogx(plt_freqs, sd_corrected_tau, syms[2], ls='-', c=cs[2], lw=2, \
             ms=6, mec=cs[2], mfc='w', mew=2, markevery=2, label=r"$\Delta\sigma_i$ adjustment")
plt.semilogx(plt_freqs, sd_mag_corrected_tau, syms[3], ls='-', c=cs[3], lw=2, \
             ms=6, mec=cs[3], mfc='w', mew=2, markevery=2, label=r"$\Delta\sigma_{i,M}$ adjustment")
plt.semilogx(plt_freqs, delta_m_tau, syms[4], ls='-', c=cs[4], lw=2, \
             ms=6, mec=cs[4], mfc='w', mew=2, markevery=2, label=r"M($\Delta$" + '$\mathregular{M_k}$' + ') adjustment')
plt.semilogx(plt_freqs, clust_sd_corrected_tau, syms[5], ls='-', c=cs[5], lw=2, \
             ms=6, mec=cs[5], mfc='w', mew=2, markevery=2, label=r"$\Delta\sigma_k$ adjustment")

#plt.xlabel('Frequency (Hz)', fontsize=15)
plt.ylabel(r"$\tau$", weight='bold', fontsize=16)
plt.legend(loc=2, fontsize=9.5, ncol=3, numpoints=1)
plt.grid(which='both')
plt.xlim([0.3, 10])
plt.ylim([0.2, 0.7])


# plot f-test p-value
plt.subplot(212)
plt.cla()
plt.semilogx([plt_freqs[0], plt_freqs[-1]], [0.05, 0.05], 'k--', lw=2)
plt.semilogx(plt_freqs, mconv_ftest, syms[1], ls='-', c=cs[1], lw=2, \
             ms=6, mec=cs[1], mfc='w', mew=2, markevery=2)
plt.semilogx(plt_freqs, sd_corrected_ftest, syms[2], ls='-', c=cs[2], lw=2, \
             ms=6, mec=cs[2], mfc='w', mew=2, markevery=2)
plt.semilogx(plt_freqs, sd_mag_corrected_ftest, syms[3], ls='-', c=cs[3], lw=2, \
             ms=6, mec=cs[3], mfc='w', mew=2, markevery=2)
plt.semilogx(plt_freqs, delta_m_ftest, syms[4], ls='-', c=cs[4], lw=2, \
             ms=6, mec=cs[4], mfc='w', mew=2, markevery=2)
plt.semilogx(plt_freqs, clust_sd_corrected_ftest, syms[5], ls='-', c=cs[5], lw=2, \
             ms=6, mec=cs[5], mfc='w', mew=2, markevery=2)

plt.xlabel('Frequency (Hz)', fontsize=14)
plt.ylabel('f-test p-value', fontsize=14)
#plt.legend(loc=2, fontsize=10, ncol=2, numpoints=1)
plt.grid(which='both')
plt.xlim([0.3, 10])
plt.ylim([-0.05, 1.05])

plt.tight_layout()
plt.savefig('figures/freq_vs_tau_f-test_'+str(maxdist)+'.png', format='png', dpi=300, bbox_inches='tight')       
plt.show()


# plot dBe vs SD regression
plt.clf()
plt.cla()
#plt.rcdefaults()
plt.figure(2, figsize=(4,2))
fig = plt.gcf()
fig.set_size_inches(4, 2.5)
plt.rc('xtick',labelsize=11)
plt.rc('ytick',labelsize=11)
plt.semilogx(plt_freqs, dBe_sd_r2, 'ko')
plt.grid(which='both')
plt.ylabel('$\mathregular{r^2}$', fontsize=14)
plt.xlabel('Frequency (Hz)', fontsize=14)
plt.xlim([0.2, 10])
plt.ylim([0, 0.9])

plt.savefig('figures/r_squared_vs_dBe_SD_reg_'+str(maxdist)+'.png', format='png', dpi=300, bbox_inches='tight')       
plt.show()

print('Mmin', mean(array(mmin)))
print('Mmax', mean(array(mmax)))
