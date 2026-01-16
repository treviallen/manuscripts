import pickle
from numpy import unique, array, arange, log, log10, exp, mean, nanmean, ndarray, sqrt, \
                  nanmedian, hstack, pi, nan, isnan, interp, where, zeros_like, polyfit, \
                  hstack, nanstd, loadtxt, nanmedian
from misc_tools import get_binned_stats, dictlist2array, get_mpl2_colourlist, savitzky_golay, get_log_xy_locs
from mag_tools import nsha18_mb2mw, nsha18_ml2mw
from get_mag_dist_terms import get_distance_term, get_magnitude_term, get_kappa_term, get_regional_term
from scipy.stats import linregress
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
    tmp = {'datetime':UTCDateTime(dat[0]), 'mw':float(dat[8]), 'qual':int(float(dat[-2])), \
           'lon':float(dat[2]), 'lat':float(dat[3]), 'sd':float(dat[10]), 'cluster':int(float(dat[-1]))}
    
    brunedat.append(tmp)
    
def get_brune_deets(fft_datetime):
    bruneStats = {'mw':nan, 'qual':0}
    for evnum, ev in enumerate(brunedat): 
        #ev['datetime'] = UTCDateTime(2024,2,27,16,4,9)
        
        if fft_datetime > UTCDateTime(ev['datetime']-timedelta(seconds=601)) \
           and fft_datetime < UTCDateTime(ev['datetime']+timedelta(seconds=300)):
            
               bruneStats = ev
    if isnan(bruneStats['mw']):
        print('Event not found: ',fft_datetime)
                       
    return bruneStats

# load atten coeffs
print('../../2023/au_stress_drop/atten_coeffs_1.3_5km.pkl')
coeff_pkl = argv[1]
coeffs = pickle.load(open(coeff_pkl, 'rb' ))
 
###############################################################################
# do ml-mw conversions
###############################################################################
 
s0, s1, s2 = loadtxt('../../2023/au_stress_drop/mw-ml_coeffs_2800.csv', delimiter=',', skiprows=1)   

cdata = loadtxt('ml_mw_bias_cluster.csv', delimiter=',', skiprows=1) 
cmean = cdata[:,1]

def get_ml2mw_cluster(ml, cluster=0):
    #print(ml)
    return s0 * ml**2 + s1 * ml + s2 - cmean[cluster+1], cmean[cluster+1] # add 1 because first row is all data


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
def get_inter_event_terms(magType, recs):
    '''
    magType: 1 = Brune mags
             2 = ML-MW cluster correction (MW - dMW) where dMW from mean MW - MW(conv)
             3 = Use converted magnitude
    '''
    resDict = []
    for f, c in enumerate(coeffs):
    
        if f == 69 or f == 99:
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
                if rec['net'] in keep_nets and rec['qual'] == 1:
                    
                    if not rec['sta'] in ignore_stas:
                        try:
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
                                    prefmw = rec['mag'] - cor 
                                    #print(rec['mag'], cor, prefmw, rec['cluster'],rec['evdt'])
                                    
                                elif magType == 3:
                                    prefmw, cor = get_ml2mw_cluster(rec['ml2800'], -1) # subtract 1 as we want index 0
                                    
                                    
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
                    
                                yobs = log10(rec[channel]['swave_spec'][f])
                                yres.append(yobs - ypred)
                                rmags.append(rec['mag'])
                                revent.append(rec['ev'])
                                revdt.append(rec['evdt'])
                                rsd.append(rec['sd'])
                        except:
                            # do nothing
                            dummy = 0
                            
                        
            resData = {'yres':array(yres), 'mags':array(rmags), 'ev':array(revent), \
                       'rhyps':array(rhyps), 'rsd':array(rsd), 'evdt':array(revdt)}
            
            resDict.append(resData)
                    
    return resDict
    
###############################################################################
# calculate different scenarios
###############################################################################
'''
# raw
print('Getting raw residual data ...')
resDict = get_inter_event_terms(1, recs)
pklfile = open('residual_data.pkl', 'wb')
pickle.dump(resDict, pklfile, protocol=-1)
pklfile.close()

#  mw using regional correction
print('Getting cluster mag residual data ...')
resDict = get_inter_event_terms(2, recs)
pklfile = open('residual_ml_cluster_data.pkl', 'wb')
pickle.dump(resDict, pklfile, protocol=-1)
pklfile.close()

# ML2MW no corrections
print('Getting ml2mw residual data ...')
resDict = get_inter_event_terms(3, recs)
pklfile = open('residual_ml_nocluster_data.pkl', 'wb')
pickle.dump(resDict, pklfile, protocol=-1)
pklfile.close()
'''

resDictRaw = pickle.load(open('residual_data.pkl', 'rb' ))
resDictCluster = pickle.load(open('residual_ml_cluster_data.pkl', 'rb' ))
resDictNoCluster = pickle.load(open('residual_ml_nocluster_data.pkl', 'rb' ))

###############################################################################
# get within event sigmas
###############################################################################
def get_inter_event_residuals(resDict):
    
    sigmaDict = []
    sigma_be = []
    sigma_we = []
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
        
        sigma_be.append(nanstd(array(event_terms)))
        sigma_we.append(nanstd(array(yres_weterm)))    
        
        # keep some data for plotting
        if i == 0:
            plt_event_terms1 = event_terms
            plt_event_mags1 = event_mags
            plt_event_sds1 = event_sds
            plt_yres_weterm1 = yres_weterm
            plt_mags1 = mags
            plt_rhyps1 = rhyps
            plt_events1 = uevdt
        elif i == 1:
            plt_event_terms2 = event_terms
            plt_event_mags2 = event_mags
            plt_event_sds2 = event_sds
            plt_yres_weterm2 = yres_weterm
            plt_mags2 = mags
            plt_rhyps2 = rhyps
            plt_events2 = uevdt
    
    return plt_events1, plt_events2, array(plt_event_terms1), array(plt_event_mags1), array(plt_event_sds1), \
           array(plt_event_terms2), array(plt_event_mags2), array(plt_event_sds2), array(sigma_be), \
           array(plt_rhyps1), array(plt_rhyps2), array(plt_yres_weterm1), array(plt_yres_weterm2), \
           array(plt_mags1), array(plt_mags2)

###############################################################################
# get stress drop correction
###############################################################################
fig = plt.figure(1, figsize=(14,4))
plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)
#plt.tight_layout() 
props = dict(boxstyle='round', facecolor='w', alpha=1)

pltlett = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)']

plt_events1, plt_events2, plt_event_terms1, plt_event_mags1, plt_event_sds1, plt_event_terms2, \
plt_event_mags2, plt_event_sds2, sigma_be, plt_rhyps1, plt_rhyps2, plt_yres_weterm1, plt_yres_weterm2, plt_mags1, plt_mags2 \
    = get_inter_event_residuals(resDictRaw)
    
ax = plt.subplot(1,2,1)
norm2 = mpl.colors.Normalize(vmin=3, vmax=7)
plt.scatter(plt_event_sds1, plt_event_terms1, c=plt_event_mags1, marker='o', s=30, edgecolor='none', cmap='viridis_r', norm=norm2, alpha=1)
#plt.scatter(plt_event_sds1, plt_event_terms1, c=clusters, marker='o', s=30, edgecolor='none', cmap='viridis_r', norm=norm2, alpha=1)
ax.set_xscale('log')
plt.xlim([0.05, 100])
plt.ylim([-2, 2])
plt.grid(which='both')
plt.ylabel('Between-Event\n(ln Residual)', fontsize=14)
plt.xlabel('Stress Drop (MPa)', fontsize=14)

# regress
idx = where(isnan(plt_event_terms1) == False)[0]
sdreg = linregress(log10(plt_event_sds1[idx]), plt_event_terms1[idx])
print(sdreg)
xfit = array([log10(0.05), 2])
yfit = sdreg[0] * xfit + sdreg[1]
plt.semilogx(10**xfit, yfit, 'k--', lw = 3)

insettxt = '$\mathregular{r^2}$' + ' = '+  str('%0.2f' % sdreg[2])

xpos = get_log_xy_locs([0.05, 100], 0.04)
ypos = (4*0.94) - 2
plt.text(xpos, ypos, insettxt, ha='left', va='top', fontsize=14, bbox=props)

# correct event terms for SD
plt_event_terms_sd_corrected1 = plt_event_terms1 - (sdreg[0] * log10(plt_event_sds1) + sdreg[1])

ax = plt.subplot(1,2,2)
norm2 = mpl.colors.Normalize(vmin=3, vmax=7)
plt.scatter(plt_event_sds2, plt_event_terms2, c=plt_event_mags2, marker='o', s=30, edgecolor='none', cmap='viridis_r', norm=norm2, alpha=1)
ax.set_xscale('log')
plt.xlim([0.05, 100])
plt.ylim([-2, 2])
plt.grid(which='both')
plt.xlabel('Stress Drop (MPa)', fontsize=14)

# regress
idx = where(isnan(plt_event_terms2) == False)[0]
sdreg = linregress(log10(plt_event_sds2[idx]), plt_event_terms2[idx])
print(sdreg)
xfit = array([log10(0.05), 2])
yfit = sdreg[0] * xfit + sdreg[1]
plt.semilogx(10**xfit, yfit, 'k--', lw = 3)

insettxt = '$\mathregular{r^2}$' + ' = '+  str('%0.2f' % sdreg[2])

xpos = get_log_xy_locs([0.05, 100], 0.04)
ypos = (4*0.94) - 2
plt.text(xpos, ypos, insettxt, ha='left', va='top', fontsize=14, bbox=props)

# set cbars
cax = fig.add_axes([0.91,0.15,0.02,0.7]) # setup colorbar axes.
cb = colorbar.ColorbarBase(cax, cmap=plt.cm.viridis_r, orientation='vertical', alpha=1, norm=norm2)
ticks = [3.0, 4.0, 5.0, 6.0, 7.0]
cb.set_ticks(ticks)
cb.set_ticklabels([str(x) for x in ticks])
cb.set_label('Moment Magnitude', rotation=270, labelpad=20, fontsize=15)

plt.savefig('figures/event_term_vs_stressdrop.png', fmt='png', dpi=300, bbox_inches='tight')       
plt.show()

# correct event terms for SD
plt_event_terms_sd_corrected2 = plt_event_terms2 - (sdreg[0] * log10(plt_event_sds2) + sdreg[1])

###############################################################################
# get cluster correction
###############################################################################
# get clusters
clusters1 = []
clust_ev = []
clust_lon = []
clust_logsd = []
for plt_event in plt_events1:
    bruneStats = get_brune_deets(UTCDateTime(plt_event))
    clusters1.append(bruneStats['cluster'])
    
    
clusters2 = []
for plt_event in plt_events2:
    bruneStats = get_brune_deets(UTCDateTime(plt_event))
    clusters2.append(bruneStats['cluster'])
    clust_ev.append(bruneStats['datetime'])
    clust_lon.append(bruneStats['lon'])
    clust_logsd.append(log10(bruneStats['sd']))

clusters1 = array(clusters1)
clusters2 = array(clusters2)    
uclusters = unique(clusters2)
clust_logsd = array(clust_logsd)

cluster_correction1 = []
cluster_correction2 = []
for cluster in uclusters:
    idx = where(clusters1 == cluster)[0]
    meancluster = nanmean(plt_event_terms1[idx])
    cluster_correction1.append(meancluster)
    
    idx = where(clusters2 == cluster)[0]
    meancluster = nanmean(plt_event_terms2[idx])
    cluster_correction2.append(meancluster)
    meanlogstress = nanmean(clust_logsd[idx])
    print(cluster+1, meancluster, meanlogstress, len(clust_logsd[idx]))

# now make regional correction    
plt_event_terms_cluster_cor1 = plt_event_terms1.copy()
plt_event_terms_cluster_cor2 = plt_event_terms2.copy()

for i, cluster in enumerate(uclusters):
    idx = where(clusters1 == cluster)[0]
    plt_event_terms_cluster_cor1[idx] = plt_event_terms1[idx] - cluster_correction1[i]
    
    idx = where(clusters2 == cluster)[0]
    plt_event_terms_cluster_cor2[idx] = plt_event_terms2[idx] - cluster_correction2[i]
    print(plt_event_terms2[idx])
    print(cluster, cluster_correction2[i])
    print(plt_event_terms_cluster_cor2[idx])
    #print(clust_lon[idx])

###############################################################################
# plot ergodic inter-event residuals
###############################################################################
fig = plt.figure(2, figsize=(14,8))
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
#plt.tight_layout() 
props = dict(boxstyle='round', facecolor='w', alpha=1)

pltlett = ['a)', 'b)', 'c)', 'd)', 'e)', 'f)']

ax = plt.subplot(2,2,1)
norm1 = mpl.colors.Normalize(vmin=-1, vmax=1.8)
plt.scatter(plt_event_mags1, plt_event_terms1, c=log10(plt_event_sds1), marker='o', s=30, edgecolor='none', cmap='plasma_r', norm=norm1, alpha=1)
plt.plot([3, 7],[0,0], 'k--', lw=1.)

plt.xlabel('Moment Magnitude', fontsize=14)
plt.ylabel('Between-Event\n(ln Residual)', fontsize=14)

xticks = [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
ax.set_xticks(xticks)
ax.set_xticklabels([str(x) for x in xticks])

plt.xlim([3, 7])
plt.ylim([-2.5, 2.5])
insettxt = '$\mathregular{M_{Brune}}$\n' \
            + 'f = 0.75 Hz\n' \
            + r"$\tau$" + ' = ' + str('%0.2f' % sigma_be[0])

xpos = 3 + 4 * 0.97
ypos = (5*0.95) - 2.5
plt.text(xpos, ypos, insettxt, ha='right', va='top', fontsize=12, bbox=props)

xpos = 3 - (7.-3.)*0.04
ypos = 5*1.03 - 2.5
plt.text(xpos, ypos, pltlett[0], fontsize=20, va='bottom', ha='left')

###########

ax = plt.subplot(2,2,2)
plt.scatter(plt_event_mags2, plt_event_terms2, c=log10(plt_event_sds2), marker='o', s=30, edgecolor='none', cmap='plasma_r', norm=norm1, alpha=1)
plt.plot([3, 7],[0,0], 'k--', lw=1.)

plt.xlabel('Moment Magnitude', fontsize=14)
#plt.ylabel('Between-Event\n(ln Residual)', fontsize=14)

ax.set_xticks(xticks)
ax.set_xticklabels([str(x) for x in xticks])

plt.xlim([3, 7])
plt.ylim([-2.5, 2.5])
insettxt = '$\mathregular{M_{Brune}}$\n' \
            + 'f = 3.0 Hz\n' \
            + r"$\tau$" + ' = ' + str('%0.2f' % sigma_be[1])

xpos = 3 + 4 * 0.97
ypos = (5*0.95) - 2.5
plt.text(xpos, ypos, insettxt, ha='right', va='top', fontsize=12, bbox=props)

xpos = 3 - (7.-3.)*0.04
ypos = 5*1.03 - 2.5
plt.text(xpos, ypos, pltlett[1], fontsize=18, va='bottom', ha='right')

###############################################################################
# plot ergodic intra-event residuals
###############################################################################

ax = plt.subplot(2,2,3)
norm2 = mpl.colors.Normalize(vmin=3, vmax=7)
plt.scatter(plt_rhyps1, plt_yres_weterm1, c=plt_mags1, marker='o', s=30, edgecolor='none', cmap='viridis_r', norm=norm2, alpha=1)
plt.plot([0, 2000],[0,0], 'k--', lw=1.)

plt.xlabel('Hypocentral Distance (km)', fontsize=14)
plt.ylabel('Within-Event\n(ln Residual)', fontsize=14)

plt.xlim([0, 2000])
plt.ylim([-3, 3])

sigma_we = nanstd(plt_yres_weterm1)
insettxt = '$\mathregular{M_{Brune}}$\n' \
            + 'f = 0.75 Hz\n' \
            + r"$\phi$" + ' = ' + str('%0.2f' % sigma_we)

xpos = 2000 * 0.97
ypos = (6*0.95) - 3
plt.text(xpos, ypos, insettxt, ha='right', va='top', fontsize=12, bbox=props)

xpos = -2000*0.04
ypos = 6*1.05 - 3.
plt.text(xpos, ypos, pltlett[2], fontsize=18, va='bottom', ha='right')

'''
xticks = [0, 500, 700, 1000, 1500, 2000]
ax.set_xticks(xticks)
ax.set_xticklabels([str(x) for x in xticks])
'''
############

ax = plt.subplot(2,2,4)
plt.scatter(plt_rhyps2, plt_yres_weterm2, c=plt_mags2, marker='o', s=30, edgecolor='none', cmap='viridis_r', norm=norm2, alpha=1)
plt.plot([0, 2000],[0,0], 'k--', lw=1.)

plt.xlabel('Hypocentral Distance (km)', fontsize=14)
#plt.ylabel('Within-Event\n(ln Residual)', fontsize=14)

plt.xlim([0, 2000])
plt.ylim([-3, 3])

sigma_we = nanstd(plt_yres_weterm2)
insettxt = '$\mathregular{M_{Brune}}$\n' \
            + 'f = 3.0 Hz\n' \
            + r"$\phi$" + ' = ' + str('%0.2f' % sigma_we)

xpos = 2000 * 0.97
ypos = (6*0.95) - 3
plt.text(xpos, ypos, insettxt, ha='right', va='top', fontsize=12, bbox=props)

xpos = -2000*0.04
ypos = 6*1.05 - 3.
plt.text(xpos, ypos, pltlett[3], fontsize=18, va='bottom', ha='right')

'''
xticks = [0, 500, 700, 1000, 1500, 2000]
ax.set_xticks(xticks)
ax.set_xticklabels([str(x) for x in xticks])
'''
plt.subplots_adjust(hspace=0.3)

#############
# make colorbars

plt.rc('xtick',labelsize=12)
plt.rc('ytick',labelsize=12)
cax = fig.add_axes([0.91,0.58,0.02,0.3]) # setup colorbar axes.
cb = colorbar.ColorbarBase(cax, cmap=plt.cm.plasma_r, orientation='vertical', alpha=1, norm=norm1)
ticks = [-1. , -0.5,  0. ,  0.5,  1. ,  1.5]
cb.set_ticks(ticks)
cb.set_ticklabels([str('%0.1f' % 10**x) for x in ticks])
cb.set_label('Stress Drop (MPa)', rotation=270, labelpad=20, fontsize=14)

# set cbars
cax = fig.add_axes([0.91,0.12,0.02,0.3]) # setup colorbar axes.
cb = colorbar.ColorbarBase(cax, cmap=plt.cm.viridis_r, orientation='vertical', alpha=1, norm=norm2)
ticks = [3.0, 4.0, 5.0, 6.0, 7.0]
cb.set_ticks(ticks)
cb.set_ticklabels([str(x) for x in ticks])
cb.set_label('Moment Magnitude', rotation=270, labelpad=20, fontsize=14)

plt.savefig('figures/ergodic_within-between_residuals.png', fmt='png', dpi=300, bbox_inches='tight')       
plt.show()

###############################################################################
# plot event-specific inter-event residuals 
###############################################################################
# gget data using MW converted from ML
plt_events1, plt_events2, plt_event_terms1, plt_event_mags1, plt_event_sds1, plt_event_terms2, \
plt_event_mags2, plt_event_sds2, sigma_be, plt_rhyps1, plt_rhyps2, plt_yres_weterm1, plt_yres_weterm2, plt_mags1, plt_mags2 \
    = get_inter_event_residuals(resDictNoCluster)

fig = plt.figure(3, figsize=(14,8))
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)

ax = plt.subplot(2,2,1)
norm1 = mpl.colors.Normalize(vmin=-1, vmax=1.8)
plt.scatter(plt_event_mags1, plt_event_terms1, c=log10(plt_event_sds1), marker='o', s=30, edgecolor='none', cmap='plasma_r', norm=norm1, alpha=1)
plt.plot([3, 7],[0,0], 'k--', lw=1.)

plt.ylabel('Between-Event\n(ln Residual)', fontsize=14)

xticks = [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
ax.set_xticks(xticks)
ax.set_xticklabels([str(x) for x in xticks])

plt.xlim([3, 7])
plt.ylim([-2.5, 2.5])

# = nanstd(plt_event_terms_sd_corrected1)

plt.xlim([3, 7])
plt.ylim([-2.5, 2.5])
insettxt = '$\mathregular{M_{conv}}$\n' \
            + 'f = 0.75 Hz\n' \
            + r"$\tau$" + ' = ' + str('%0.2f' % sigma_be[0])

xpos = 3 + 4 * 0.97
ypos = (5*0.95) - 2.5
plt.text(xpos, ypos, insettxt, ha='right', va='top', fontsize=12, bbox=props)

xpos = 3 - (7.-3.)*0.04
ypos = 5*1.03 - 2.5
plt.text(xpos, ypos, pltlett[0], fontsize=18, va='bottom', ha='right')


###########

ax = plt.subplot(2,2,2)
plt.scatter(plt_event_mags2, plt_event_terms2, c=log10(plt_event_sds2), marker='o', s=30, edgecolor='none', cmap='plasma_r', norm=norm1, alpha=1)
plt.plot([3, 7],[0,0], 'k--', lw=1.)

#plt.ylabel('Between-Event\n(ln Residual)', fontsize=14)

ax.set_xticks(xticks)
ax.set_xticklabels([str(x) for x in xticks])

plt.xlim([3, 7])
plt.ylim([-2.5, 2.5])

insettxt = '$\mathregular{M_{conv}}$\n' \
            + 'f = 3.0 Hz\n' \
            + r"$\tau$" + ' = ' + str('%0.2f' % sigma_be[1])

xpos = 3 + 4 * 0.97
ypos = (5*0.95) - 2.5
plt.text(xpos, ypos, insettxt, ha='right', va='top', fontsize=12, bbox=props)

xpos = 3 - (7.-3.)*0.04
ypos = 5*1.03 - 2.5
plt.text(xpos, ypos, pltlett[1], fontsize=18, va='bottom', ha='right')


####################
# do SD-inter-event

# re-read ergodic data to correct
plt_events1, plt_events2, plt_event_terms1, plt_event_mags1, plt_event_sds1, plt_event_terms2, \
plt_event_mags2, plt_event_sds2, sigma_be, plt_rhyps1, plt_rhyps2, plt_yres_weterm1, plt_yres_weterm2, plt_mags1, plt_mags2 \
    = get_inter_event_residuals(resDictRaw)


ax = plt.subplot(2,2,3)
norm1 = mpl.colors.Normalize(vmin=-1, vmax=1.8)
plt.scatter(plt_event_mags1, plt_event_terms_sd_corrected1, c=log10(plt_event_sds1), marker='o', s=30, edgecolor='none', cmap='plasma_r', norm=norm1, alpha=1)
plt.plot([3, 7],[0,0], 'k--', lw=1.)

plt.xlabel('Moment Magnitude', fontsize=14)
plt.ylabel('Between-Event\n(ln Residual)', fontsize=14)

xticks = [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
ax.set_xticks(xticks)
ax.set_xticklabels([str(x) for x in xticks])

plt.xlim([3, 7])
plt.ylim([-2.5, 2.5])
sigma = nanstd(plt_event_terms_sd_corrected1)
# r"raw string Mathtext: $\alpha > \beta$"
insettxt = r"$\Delta\sigma_i$ adjustment" + '\n' \
            + 'f = 0.75 Hz\n' \
            + r"$\tau$" + ' = ' + str('%0.2f' % sigma)

xpos = get_log_xy_locs([3, 7], 0.98)
ypos = (5*0.95) - 2.5
plt.text(xpos, ypos, insettxt, ha='right', va='top', fontsize=12, bbox=props)

xpos = 3 - (7.-3.)*0.04
ypos = 5*1.03 - 2.5
plt.text(xpos, ypos, pltlett[2], fontsize=18, va='bottom', ha='left')

###########

ax = plt.subplot(2,2,4)
plt.scatter(plt_event_mags2, plt_event_terms_sd_corrected2, c=log10(plt_event_sds2), marker='o', s=30, edgecolor='none', cmap='plasma_r', norm=norm1, alpha=1)
plt.plot([3, 7],[0,0], 'k--', lw=1.)

plt.xlabel('Moment Magnitude', fontsize=14)
#plt.ylabel('Between-Event\n(ln Residual)', fontsize=14)

ax.set_xticks(xticks)
ax.set_xticklabels([str(x) for x in xticks])

plt.xlim([3, 7])
plt.ylim([-2.5, 2.5])
sigma = nanstd(plt_event_terms_sd_corrected2)
insettxt = r"$\Delta\sigma_i$ adjustment" + '\n' \
            + 'f = 3.0 Hz\n' \
            + r"$\tau$" + ' = ' + str('%0.2f' % sigma)

xpos = get_log_xy_locs([3, 7], 0.98)
ypos = (5*0.95) - 2.5
plt.text(xpos, ypos, insettxt, ha='right', va='top', fontsize=12, bbox=props)

xpos = 3 - (7.-3.)*0.04
ypos = 5*1.03 - 2.5
plt.text(xpos, ypos, pltlett[3], fontsize=18, va='bottom', ha='right')

plt.subplots_adjust(hspace=0.3)

# set cbars
plt.rc('xtick',labelsize=13)
plt.rc('ytick',labelsize=13)
cax = fig.add_axes([0.93,0.25,0.02,0.5]) # setup colorbar axes.
cb = colorbar.ColorbarBase(cax, cmap=plt.cm.plasma_r, orientation='vertical', alpha=1, norm=norm1)
ticks = [-1. , -0.5,  0. ,  0.5,  1. ,  1.5]
cb.set_ticks(ticks)
cb.set_ticklabels([str('%0.1f' % 10**x) for x in ticks])
cb.set_label('Stress Drop (MPa)', rotation=270, labelpad=20, fontsize=15)


plt.savefig('figures/interevent_event-specific.png', fmt='png', dpi=300, bbox_inches='tight')       
plt.show()

###############################################################################
# plot inter-event residuals with clustering
###############################################################################

plt_events1, plt_events2, plt_event_terms1, plt_event_mags1, plt_event_sds1, plt_event_terms2, \
plt_event_mags2, plt_event_sds2, sigma_be, plt_rhyps1, plt_rhyps2, plt_yres_weterm1, plt_yres_weterm2, plt_mags1, plt_mags2 \
    = get_inter_event_residuals(resDictCluster)

fig = plt.figure(4, figsize=(14,8))
plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)

ax = plt.subplot(2,2,1)
norm1 = mpl.colors.Normalize(vmin=-1, vmax=1.8)
plt.scatter(plt_event_mags1, plt_event_terms1, c=log10(plt_event_sds1), marker='o', s=30, edgecolor='none', cmap='plasma_r', norm=norm1, alpha=1)
plt.plot([3, 7],[0,0], 'k--', lw=1.)

plt.ylabel('Between-Event\n(ln Residual)', fontsize=14)

xticks = [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
ax.set_xticks(xticks)
ax.set_xticklabels([str(x) for x in xticks])

plt.xlim([3, 7])
plt.ylim([-2.5, 2.5])

# = nanstd(plt_event_terms_sd_corrected1)

plt.xlim([3, 7])
plt.ylim([-2.5, 2.5])
insettxt = r"M($\Delta$" + '$\mathregular{M_k}$' + ') adjustment\n' \
            + 'f = 0.75 Hz\n' \
            + r"$\tau$" + ' = ' + str('%0.2f' % sigma_be[0])

xpos = 3 + 4 * 0.97
ypos = (5*0.95) - 2.5
plt.text(xpos, ypos, insettxt, ha='right', va='top', fontsize=12, bbox=props)

xpos = 3 - (7.-3.)*0.04
ypos = 5*1.03 - 2.5
plt.text(xpos, ypos, pltlett[0], fontsize=18, va='bottom', ha='right')


###########

ax = plt.subplot(2,2,2)
plt.scatter(plt_event_mags2, plt_event_terms2, c=log10(plt_event_sds2), marker='o', s=30, edgecolor='none', cmap='plasma_r', norm=norm1, alpha=1)
plt.plot([3, 7],[0,0], 'k--', lw=1.)

#plt.ylabel('Between-Event\n(ln Residual)', fontsize=14)

ax.set_xticks(xticks)
ax.set_xticklabels([str(x) for x in xticks])

plt.xlim([3, 7])
plt.ylim([-2.5, 2.5])

insettxt = r"M($\Delta$" + '$\mathregular{M_k}$' + ') adjustment\n' \
            + 'f = 3.0 Hz\n' \
            + r"$\tau$" + ' = ' + str('%0.2f' % sigma_be[1])

xpos = 3 + 4 * 0.97
ypos = (5*0.95) - 2.5
plt.text(xpos, ypos, insettxt, ha='right', va='top', fontsize=12, bbox=props)

xpos = 3 - (7.-3.)*0.04
ypos = 5*1.03 - 2.5
plt.text(xpos, ypos, pltlett[1], fontsize=18, va='bottom', ha='right')


####################
# do SD-inter-event
# re-read ergodic data to correct
plt_events1, plt_events2, plt_event_terms1, plt_event_mags1, plt_event_sds1, plt_event_terms2, \
plt_event_mags2, plt_event_sds2, sigma_be, plt_rhyps1, plt_rhyps2, plt_yres_weterm1, plt_yres_weterm2, plt_mags1, plt_mags2 \
    = get_inter_event_residuals(resDictRaw)

ax = plt.subplot(2,2,3)
norm1 = mpl.colors.Normalize(vmin=-1, vmax=1.8)
plt.scatter(plt_event_mags1, plt_event_terms_cluster_cor1, c=log10(plt_event_sds1), marker='o', s=30, edgecolor='none', cmap='plasma_r', norm=norm1, alpha=1)
plt.plot([3, 7],[0,0], 'k--', lw=1.)

plt.xlabel('Moment Magnitude', fontsize=14)
plt.ylabel('Between-Event\n(ln Residual)', fontsize=14)

xticks = [3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
ax.set_xticks(xticks)
ax.set_xticklabels([str(x) for x in xticks])

plt.xlim([3, 7])
plt.ylim([-2.5, 2.5])
sigma = nanstd(plt_event_terms_cluster_cor1)
# r"raw string Mathtext: $\alpha > \beta$"
insettxt = r"$\Delta\sigma_k$ adjustment" + '\n' \
            + 'f = 0.75 Hz\n' \
            + r"$\tau$" + ' = ' + str('%0.2f' % sigma)

xpos = get_log_xy_locs([3, 7], 0.98)
ypos = (5*0.95) - 2.5
plt.text(xpos, ypos, insettxt, ha='right', va='top', fontsize=12, bbox=props)

xpos = 3 - (7.-3.)*0.04
ypos = 5*1.03 - 2.5
plt.text(xpos, ypos, pltlett[2], fontsize=18, va='bottom', ha='left')

###########

ax = plt.subplot(2,2,4)
plt.scatter(plt_event_mags2, plt_event_terms_cluster_cor2, c=log10(plt_event_sds2), marker='o', s=30, edgecolor='none', cmap='plasma_r', norm=norm1, alpha=1)
plt.plot([3, 7],[0,0], 'k--', lw=1.)

plt.xlabel('Moment Magnitude', fontsize=14)
#plt.ylabel('Between-Event\n(ln Residual)', fontsize=14)

ax.set_xticks(xticks)
ax.set_xticklabels([str(x) for x in xticks])

plt.xlim([3, 7])
plt.ylim([-2.5, 2.5])
sigma = nanstd(plt_event_terms_cluster_cor2)

a = 'M'
b = "$\Delta$"
#fstr = f"M_{{s}{a}}"

insettxt =  r"$\Delta\sigma_k$ adjustment" + '\n' \
            + 'f = 3.0 Hz\n' \
            + r"$\tau$" + ' = ' + str('%0.2f' % sigma)

xpos = get_log_xy_locs([3, 7], 0.98)
ypos = (5*0.95) - 2.5
plt.text(xpos, ypos, insettxt, ha='right', va='top', fontsize=12, bbox=props)

xpos = 3 - (7.-3.)*0.04
ypos = 5*1.03 - 2.5
plt.text(xpos, ypos, pltlett[3], fontsize=18, va='bottom', ha='right')

plt.subplots_adjust(hspace=0.3)

# set cbars
plt.rc('xtick',labelsize=13)
plt.rc('ytick',labelsize=13)
cax = fig.add_axes([0.93,0.25,0.02,0.5]) # setup colorbar axes.
cb = colorbar.ColorbarBase(cax, cmap=plt.cm.plasma_r, orientation='vertical', alpha=1, norm=norm1)
ticks = [-1. , -0.5,  0. ,  0.5,  1. ,  1.5]
cb.set_ticks(ticks)
cb.set_ticklabels([str('%0.1f' % 10**x) for x in ticks])
cb.set_label('Stress Drop (MPa)', rotation=270, labelpad=20, fontsize=15)


plt.savefig('figures/interevent_clustering.png', fmt='png', dpi=300, bbox_inches='tight')       
plt.show()
