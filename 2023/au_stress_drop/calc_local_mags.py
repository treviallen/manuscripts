import pickle
from numpy import unique, array, arange, log, log10, exp, mean, nanmean, ndarray, std, sqrt, \
                  nanmedian, nanstd, vstack, pi, nan, isnan, interp, where, zeros_like, ones_like, floor, ceil, \
                  argsort, loadtxt, percentile
from scipy.stats import linregress, trim_mean
from misc_tools import get_binned_stats, dictlist2array, get_mpl2_colourlist, get_log_xy_locs
from mag_tools import nsha18_mb2mw, nsha18_ml2mw, get_au_ml_zone
from mag_tools import m02mw
from calculate_magnitudes import calc_R35, calc_HB87, calc_MLM92, calc_BJ84, calc_GG91, calc_GS86
from mapping_tools import distance
from scipy.odr import Data, Model, ODR, models
import scipy.odr.odrpack as odrpack
from ltsfit import lts_linefit
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import colorbar, colors #, cm
import cmcrameri.cm as cmc
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
# load surface ruptures
sface_csv = 'au_surface_ruptures.csv'
lines = open(sface_csv).readlines()[1:]
sface_ruptures = []
for line in lines:
    dat = line.strip().split(',')
    sface_ruptures.append(dat[0].strip('\n'))

def parse_filtering_data():
    
    mwdat = []
    # read parameter file
    lines = open('brune_stats.csv').readlines()[1:]
    for line in lines:
        dat = line.split(',')
        filt = {'ev':dat[0], 'minf': float(dat[-3]), 'maxf': float(dat[-2]), 'qual':float(dat[-1]),
        	      'mw': float(dat[8]), 'sd': float(dat[10]), 'sd_std': float(dat[11]), 'surf_rup':0}
        
        # check if sface rupture
        for sr in sface_ruptures:
            if sr == dat[0]:
                filt['surf_rup'] = 1
                print('surf_rup')
    
        mwdat.append(filt)
    
    return mwdat
    
# load surface ruptures
sface_csv = 'au_surface_ruptures.csv'
lines = open(sface_csv).readlines()[1:]
sface_ruptures = []
for line in lines:
    dat = line.split(',')
    sface_ruptures.append(dat[0])
   
####################################################################################
# set def params
####################################################################################


###############################################################################
# load pickle 
###############################################################################

#recs = pickle.load(open('onger_fft_data.pkl', 'rb' ))
recs = pickle.load(open('wa_data.pkl', 'rb' ))

###############################################################################
# set data to use
###############################################################################

# remove bad recs
keep_nets = set(['AU', 'IU', 'S1', 'II', 'G', 'MEL', 'ME', '2O', 'AD', 'SR', 'UM', 'AB', 'VI', 'GM' \
                 '1P', '1P', '2P', '6F', '7K', '7G', 'G', '7B', '4N', '7D', '', 'OZ', 'OA', 'WG', 'XX'])
# get stas to ignore
ignore_stas = open('sta_ignore.txt').readlines()
#ignore_stas = open('sta_ignore.test').readlines()
ignore_stas = set([x.strip() for x in ignore_stas])

####################################################################################
# start main
####################################################################################

# loop through & plot station  data
cs = get_mpl2_colourlist()
events = unique(dictlist2array(recs, 'evdt'))
stations = unique(dictlist2array(recs, 'sta'))
rhyps = dictlist2array(recs, 'rhyp')

# sort by rhyps
rhyp_sort_idx = argsort(rhyps)

# get event filters
mwdat = parse_filtering_data()

events_dict = []

sp = 0
ii = 1	
magcsv = 'EVENT,ML_2800,ML_2080,MW,SD,REG\n'
stacsv = 'EVENT,ML_REGION,STA,RHYP,ML_2800,ML_2080,MW\n'

print('need to get mag region')

ml_array = []
mw_array = []
sd_array = []
sd_std_array = []
mzone_array = []

for e, event in enumerate(events): # [::-1]): #[-2:-1]:
    print(event)
    
    # get upper & lower f for filtering
    qual = 0
    mw = nan
    sd = nan
    sd_std = nan
    for fdat in mwdat:
        if fdat['ev'] == event:
            minf = fdat['minf']
            maxf = fdat['maxf']
            qual = fdat['qual']
            if qual == 1.0:
                mw = fdat['mw']
                sd = fdat['sd']
                sd_std = fdat['sd_std']
            else:
                mw = nan
                sd = nan

    i = 0
    pltboxes = True
    
    bj84_2800 = []
    bj84_2080 = []
    mlm92_2800 = []
    mlm92_2080 = []
    gg91_2800 = []
    gg91_2080 = []
    gs86_2800 = []
    gs86_2080 = []
    evstas = []
    
    for rsi in rhyp_sort_idx:
        rec = recs[rsi]
        #for rec in recs:
        if rec['net'] in keep_nets:
            if not rec['sta'] in ignore_stas:
                if len(rec['channels']) > 0:
                    if rec['evdt'] == event and rec['rhyp'] <= 800:
                        print('   '+rec['sta'])
                        
                        # get NSHA23 preferred amp
                        if rec['mag'] < 4.0:
                            wa2800 = rec['wa_data']['wa_amp_2800_0.75']
                            wa2080 = rec['wa_data']['wa_amp_2080_0.75']
                        
                        elif rec['mag'] < 6.0:
                            wa2800 = rec['wa_data']['wa_amp_2800_0.5']
                            wa2080 = rec['wa_data']['wa_amp_2080_0.5']
                            
                        else:
                            try:
                                wa2800 = rec['wa_data']['wa_amp_2800_0.02']
                                wa2080 = rec['wa_data']['wa_amp_2080_0.02']
                            except:
                                wa2800 = rec['wa_data']['wa_amp_2800_0.2']
                                wa2080 = rec['wa_data']['wa_amp_2080_0.2']
                            
                        #bj84_2800.append(calc_BJ84(1, log10(wa2800), rec['rhyp']) + 0.18) # added 0.18 as mean H-V correction - see hv_ratio.png in dropbox
                        m92_2800 = calc_MLM92(0, log10(wa2800), rec['rhyp'])
                        mlm92_2800.append(m92_2800)
                        
                        #bj84_2080.append(calc_BJ84(1, log10(wa2080), rec['rhyp']) + 0.18) # added 0.18 as mean H-V correction - see hv_ratio.png in dropbox
                        m92_2080 = calc_MLM92(0, log10(wa2080), rec['rhyp'])
                        mlm92_2080.append(m92_2080)
                        
                        # get gg91
                        tmp = calc_GG91(0, log10(wa2800), rec['rhyp'])
                        gg91_2800.append(tmp)
                        
                        tmp = calc_GG91(0, log10(wa2080), rec['rhyp'])
                        gg91_2080.append(tmp)
                        
                        # get gs86
                        tmp = calc_GS86(0, log10(wa2800), rec['repi'])
                        gs86_2800.append(tmp)
                        
                        tmp = calc_GS86(0, log10(wa2080), rec['repi'])
                        gs86_2080.append(tmp)
                        
                        evstas.append(rec['sta'])
                        
                        # get mag zone 
                        mzone = get_au_ml_zone([rec['eqlo']], [rec['eqla']])
                        #print(mzone)
                        
                        if mzone[0] == 'EA':
                            staml_2800 = mlm92_2800[-1]
                            staml_2080 = mlm92_2080[-1]
                        elif mzone[0] == 'WCA':
                            staml_2800 = gg91_2800[-1] + 0.13 # correct for V2H
                            staml_2080 = gg91_2080[-1] + 0.13 # correct for V2H
                        elif mzone[0] == 'SA':
                            staml_2800 = gs86_2800[-1]
                            staml_2080 = gs86_2080[-1]
                        else:
                            staml_2800 = mlm92_2800[-1]
                            staml_2080 = mlm92_2080[-1]
                        
                        stacsv += ','.join((str(event), mzone[0], rec['sta'], str('%0.1f' % rec['rhyp']), \
                                            str('%0.2f' % staml_2800), str('%0.2f' % staml_2080), \
                                            str('%0.2f' % mw))) + '\n'
                        
    # get pref mag
    if mzone[0] == 'EA':
        prefml_2800 = mlm92_2800
        prefml_2080 = mlm92_2080
        magzone = 'EA'
    elif mzone[0] == 'WCA':
        prefml_2800 = array(gg91_2800) + 0.13 # correct for V2H
        prefml_2080 = array(gg91_2080) + 0.13 # correct for V2H
        magzone = 'WCA'
    elif mzone[0] == 'SA':
        prefml_2800 = gs86_2800
        prefml_2080 = gs86_2080
        magzone = 'SA'
    else:
        prefml_2800 = mlm92_2800
        prefml_2080 = mlm92_2080
        magzone = 'EA'
    
    if len(mlm92_2080) >= 3:
        ml_2800 = trim_mean(array(prefml_2800), 0.1)
        ml_2080 = trim_mean(array(prefml_2080), 0.1)
    else:
        ml_2800 = nan
        ml_2080 = nan

    magcsv += ','.join((str(event), str('%0.2f' % ml_2800), str('%0.2f' % ml_2080), \
                        str('%0.2f' % mw), str('%0.5f' % sd), magzone[0])) + '\n'
    
    ml_array.append(ml_2800)
    mw_array.append(mw)
    sd_array.append(sd)
    sd_std_array.append(sd_std)
    mzone_array.append(mzone)

ml_array = array(ml_array)
mw_array = array(mw_array)
sd_array = array(sd_array)
mzone_array = array(mzone_array)
##########################################################################################
# plot
fig = plt.figure(1, figsize=(10.5,9))

plt.plot([2,7], [2,7], 'k--', lw=0.5, label='1:1')
#plt.plot(ml_array, mw_array, 'o', c='darkorange', label='Data')

# plt 2023 simulated (2800) relationship
xplt = arange(3, 6.81, 0.01)
s0, s1, s2 = loadtxt('mw-ml_coeffs_2800.csv', delimiter=',', skiprows=1)
yplt_nsha23 = s0 * xplt**2 + s1 * xplt + s2

cm = plt.cm.get_cmap('RdBu_r')
sc = plt.scatter(ml_array, mw_array, c=log10(sd_array), vmin=-0.7, vmax=1.7, s=36, cmap=cm, label='Data')

# fit new data
def ortho_quad_reg(c, x):
    return c[0] * x**2 + c[1] * x + c[2]

pc84 = percentile(sd_std_array, 84)
idx = where(sd_std_array <= pc84)[0]
data = odrpack.RealData(ml_array[idx], mw_array[idx])
'''
sx = interp(mw, [min(mw), max(mw)],[0.3, 0.1])
sy = sx
data = odrpack.RealData(mb[mb>=3.5], mw[mb>=3.5], sx=sx[mb>=3.5], sy=sy[mb>=3.5])
'''
# do straight quadratic
quad_reg = odrpack.Model(ortho_quad_reg)
odr = odrpack.ODR(data, quad_reg, beta0=[0.09, 0.1, 2.0])

odr.set_job(fit_type=2) #if set fit_type=2, returns the same as least squares
out = odr.run()
c0 = out.beta[0]
c1 = out.beta[1]
c2 = out.beta[2]
print('\nQuad')
out.pprint()
yplt = c0 * xplt**2 + c1 * xplt + c2

# plt NSHA23
plt.plot(xplt, yplt_nsha23, '-', lw=2, c='k', label='NSHA23 (W-A 2800)')

# plt new regression
plt.plot(xplt, yplt, '--', lw=2, c='k', label='Present Study')

plt.xlim([2.5,7])
plt.ylim([2.5,7])
plt.xlabel('$\mathregular{M_{L(2800)}}$', fontsize=18)
plt.ylabel('$\mathregular{M_{Brune}}$', fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.legend(loc=2, fontsize=14)
plt.grid(which='both')

# make colourbar
cb = plt.colorbar(sc)

logticks = cb.get_ticks()[0:]
labels = [str('%0.1f' % 10**x) for x in logticks]
labels[0] = str('%0.2f' % 10**logticks[0])
labels[1] = str('%0.2f' % 10**logticks[1])
labels[2] = str('%0.2f' % 10**logticks[2])
labels[3] = str('%0.2f' % 10**logticks[3])
labels[4] = str('%0.2f' % 10**logticks[4])

#ticks = arange(0,len(logstressbins))
cb.set_ticks(logticks)
cb.set_ticklabels(labels)

cb.set_label('Brune Stress Drop (MPa)', rotation=270, fontsize=18, labelpad=22) # 


"""                        
##########################################################################################

# export Brune data
pklfile = open('brune_data.pkl', 'wb')
pickle.dump(events, pklfile, protocol=-1)
pklfile.close()

# write csv
txt = 'EVENT,LON,LAT,DEP,OMAG,OMAG_TYPE,BRUNE_MAG,STRESS_DROP,CORN_FREQ,FMIN,FMAX,QUALITY\n'

for ev in events_dict:
    txt += ','.join((str(ev['evdt']),str(ev['lon']),str(ev['lat']),str(ev['dep']),str(ev['omag']),ev['omag_type'], \
                     str(ev['brune_mw']),str(ev['brune_sd']),str(ev['brune_f0']),str(ev['minf']),str(ev['maxf']),str(ev['qual']))) + '\n'
"""
# write to file
f = open('ml_mw_stats.csv', 'w')
f.write(magcsv)
f.close()

# write to file
f = open('ml_sta_stats.csv', 'w')
f.write(stacsv)
f.close()


# now show figs 
plt.savefig('ml_vs_mw_brune.png', fmt='png', dpi=300, bbox_inches='tight')       
plt.show()    

##########################################################################################
# plot regional

fig = plt.figure(2, figsize=(6, 18))

#plt all data
'''
plt.subplot(221)
plt.plot([2,7], [2,7], 'k--', lw=0.5, label='1:1')
plt.plot(xplt, yplt, '-', lw=2, c='k', label='2023 Simulated (W-A 2800)')
sc = plt.scatter(ml_array, mw_array, c=log10(sd_array), vmin=-1.0, vmax=2.0, s=36, cmap=cm, label='Data')
plt.title('All Data')
plt.grid(which='both')
plt.xlabel('$\mathregular{M_L}$', fontsize=18)
plt.ylabel('Brune $\mathregular{M_W}$', fontsize=18)
plt.xlim([2.5,7])
plt.ylim([2.5,7])
'''
props = dict(boxstyle='round', facecolor='w', alpha=1)
xloc = 2.68
yloc = 6.82        
#plt EA
plt.subplot(311)
plt.plot([2,7], [2,7], 'k--', lw=0.5, label='1:1')
plt.plot(xplt, yplt, '-', lw=2, c='k', label='2023 Simulated (W-A 2800)')
idx = where(mzone_array == 'EA')[0]
sc = plt.scatter(ml_array[idx], mw_array[idx], c=log10(sd_array[idx]), vmin=-0.7, vmax=1.7, s=36, cmap=cm, label='Data')
plt.text(xloc, yloc, 'Eastern Australia', va='top', ha ='left', fontsize=12, bbox=props)
plt.xlabel('$\mathregular{M_{L(MLM92)}}$', fontsize=18)
#plt.ylabel('Brune $\mathregular{M_W}$', fontsize=18)
plt.grid(which='both')
plt.xlim([2.5,7])
plt.ylim([2.5,7])
plt.text(2.0,7.4, '(a)', va='top', ha ='left', fontsize=18)

#plt WCA
plt.subplot(312)
plt.plot([2,7], [2,7], 'k--', lw=0.5, label='1:1')
plt.plot(xplt, yplt, '-', lw=2, c='k', label='2023 Simulated (W-A 2800)')
idx = where(mzone_array == 'WCA')[0]
sc = plt.scatter(ml_array[idx], mw_array[idx], c=log10(sd_array[idx]), vmin=-0.7, vmax=1.7, s=36, cmap=cm, label='Data')
plt.text(xloc, yloc, 'Western and Central Australia', va='top', ha ='left', fontsize=12, bbox=props)
plt.xlabel('$\mathregular{M_{L(GG91)}}$', fontsize=18)
plt.ylabel('$\mathregular{M_{Brune}}$', fontsize=18)
plt.grid(which='both')
plt.xlim([2.5,7])
plt.ylim([2.5,7])
plt.text(2.0,7.4, '(b)', va='top', ha ='left', fontsize=18)

#plt SA
plt.subplot(313)
plt.plot([2,7], [2,7], 'k--', lw=0.5, label='1:1')
plt.plot(xplt, yplt, '-', lw=2, c='k', label='2023 Simulated (W-A 2800)')
idx = where(mzone_array == 'SA')[0]
sc = plt.scatter(ml_array[idx], mw_array[idx], c=log10(sd_array[idx]), vmin=-0.7, vmax=1.7, s=36, cmap=cm, label='Data')
plt.text(xloc, yloc, 'South Australia', va='top', ha ='left', fontsize=12, bbox=props)
plt.xlabel('$\mathregular{M_{L(GS86)}}$', fontsize=18)
#plt.ylabel('Brune $\mathregular{M_W}$', fontsize=18)
plt.grid(which='both')
plt.xlim([2.5,7])
plt.ylim([2.5,7])
plt.text(2.0,7.4, '(c)', va='top', ha ='left', fontsize=18)

# make colourbar
#cb = plt.colorbar(sc)

#ticks = arange(0,len(logstressbins))
plt.gcf().subplots_adjust(right=0.85)
cax = fig.add_axes([0.88,0.3,0.05,0.4]) # setup colorbar axes.
norm = colors.Normalize(vmin=-0.7, vmax=1.7)
cb = colorbar.ColorbarBase(cax, cmap=cm, norm=norm, orientation='vertical')

logticks = cb.get_ticks()[0:]
labels = [str('%0.1f' % 10**x) for x in logticks]
labels[0] = str('%0.2f' % 10**logticks[0])
labels[1] = str('%0.2f' % 10**logticks[1])
labels[2] = str('%0.2f' % 10**logticks[2])
labels[3] = str('%0.2f' % 10**logticks[3])
labels[4] = str('%0.2f' % 10**logticks[4])

cb.set_ticks(logticks)
cb.set_ticklabels(labels)

cb.set_label('Brune Stress Drop (MPa)', rotation=270, fontsize=18, labelpad=18)

# now show figs 
plt.savefig('mw_vs_ml_zone.png', fmt='png', dpi=300, bbox_inches='tight')       
plt.show()    

##########################################################################################
# plot regional residual
mw_pred = c0 * ml_array**2 + c1 * ml_array + c2 # NSHA23
mw_res = mw_array - mw_pred
xloc = 2.87
yloc = 0.736        
xloc2 = 6.63
yloc2 = -0.736

fig = plt.figure(3, figsize=(14, 8.))

#plt all data
plt.subplot(221)
plt.plot([2,7], [0,0], 'k--', lw=0.5)
sc = plt.scatter(mw_array, mw_res, c=log10(sd_array), vmin=-0.7, vmax=1.7, s=36, cmap=cm, label='Data')
#sc = plt.scatter(ml_array, mw_res, c=log10(sd_array), vmin=-0.7, vmax=1.7, s=36, cmap=cm, label='Data')
plt.text(xloc, yloc, 'All Data', va='top', ha ='left', fontsize=12, bbox=props)
plt.grid(which='both')
plt.ylabel('$\mathregular{M_{Brune}}$ - $\mathregular{M_{W(Conv)}}$', fontsize=16)
plt.xlim([2.75,6.75])
plt.ylim([-0.8,0.8])
plt.text(2.3,0.98, '(a)', va='top', ha ='left', fontsize=16)

# get mean & std
mean_res = nanmean(mw_res)
med_res = nanmedian(mw_res[idx])
std_res = nanstd(mw_res)
stat_txt = r'$\mu$ = '+str('%0.2f' % mean_res) + '\n' \
           + r'$\tilde{x}$ = '+str('%0.2f' % med_res) + '\n' \
           + r'$\tau$ = '+str('%0.2f' % std_res)
plt.text(xloc2, yloc2, stat_txt, va='bottom', ha ='right', fontsize=12, bbox=props)


#plt EA
plt.subplot(222)
plt.plot([2,7], [0,0], 'k--', lw=0.5)
idx = where(mzone_array == 'EA')[0]
sc = plt.scatter(mw_array[idx], mw_res[idx], c=log10(sd_array[idx]), vmin=-0.7, vmax=1.7, s=36, cmap=cm, label='Data')
#sc = plt.scatter(ml_array[idx], mw_res[idx], c=log10(sd_array[idx]), vmin=-0.7, vmax=1.7, s=36, cmap=cm, label='Data')
plt.text(xloc, yloc, 'Eastern Australia', va='top', ha ='left', fontsize=12, bbox=props)
#plt.xlabel('$\mathregular{M_{L(MLM92)}}$', fontsize=18)
#plt.ylabel('Brune $\mathregular{M_W}$', fontsize=18)

# get mean & std
mean_res = nanmean(mw_res[idx])
med_res = nanmedian(mw_res[idx])
std_res = nanstd(mw_res[idx])
stat_txt = r'$\mu$ = '+str('%0.2f' % mean_res) + '\n' \
           + r'$\tilde{x}$ = '+str('%0.2f' % med_res) + '\n' \
           + r'$\tau$ = '+str('%0.2f' % std_res)
plt.text(xloc2, yloc2, stat_txt, va='bottom', ha ='right', fontsize=12, bbox=props)

plt.grid(which='both')
plt.xlim([2.75,6.75])
plt.ylim([-0.8,0.8])
plt.text(2.3,0.98, '(b)', va='top', ha ='left', fontsize=16)

#plt WCA
plt.subplot(223)
plt.plot([2,7], [0,0], 'k--', lw=0.5)
idx = where(mzone_array == 'WCA')[0]
sc = plt.scatter(mw_array[idx], mw_res[idx], c=log10(sd_array[idx]), vmin=-0.7, vmax=1.7, s=36, cmap=cm, label='Data')
#sc = plt.scatter(ml_array[idx], mw_res[idx], c=log10(sd_array[idx]), vmin=-0.7, vmax=1.7, s=36, cmap=cm, label='Data')
plt.text(xloc, yloc, 'Western and Central Australia', va='top', ha ='left', fontsize=12, bbox=props)

# get mean & std
mean_res = nanmean(mw_res[idx])
med_res = nanmedian(mw_res[idx])
std_res = nanstd(mw_res[idx])
stat_txt = r'$\mu$ = '+str('%0.2f' % mean_res) + '\n' \
           + r'$\tilde{x}$ = '+str('%0.2f' % med_res) + '\n' \
           + r'$\tau$ = '+str('%0.2f' % std_res)
plt.text(xloc2, yloc2, stat_txt, va='bottom', ha ='right', fontsize=12, bbox=props)

plt.xlabel('$\mathregular{M_{Brune}}$', fontsize=16)
plt.ylabel('$\mathregular{M_{Brune}}$ - $\mathregular{M_{W(Conv)}}$', fontsize=16)
plt.grid(which='both')
plt.xlim([2.75,6.75])
plt.ylim([-0.8,0.8])
plt.text(2.3,0.98, '(c)', va='top', ha ='left', fontsize=16)

#plt SA
plt.subplot(224)
plt.plot([2,7], [0,0], 'k--', lw=0.5)
idx = where(mzone_array == 'SA')[0]
sc = plt.scatter(mw_array[idx], mw_res[idx], c=log10(sd_array[idx]), vmin=-0.7, vmax=1.7, s=36, cmap=cm, label='Data')
plt.text(xloc, yloc, 'South Australia', va='top', ha ='left', fontsize=12, bbox=props)

# get mean & std
mean_res = nanmean(mw_res[idx])
med_res = nanmedian(mw_res[idx])
std_res = nanstd(mw_res[idx])
stat_txt = r'$\mu$ = '+str('%0.2f' % mean_res) + '\n' \
           + r'$\tilde{x}$ = '+str('%0.2f' % med_res) + '\n' \
           + r'$\tau$ = '+str('%0.2f' % std_res)
plt.text(xloc2, yloc2, stat_txt, va='bottom', ha ='right', fontsize=12, bbox=props)

plt.xlabel('$\mathregular{M_{Brune}}$', fontsize=16)
#plt.ylabel('$\mathregular{M_{Brune}}$', fontsize=16)#plt.ylabel('Brune $\mathregular{M_W}$', fontsize=18)
plt.grid(which='both')
plt.xlim([2.75,6.75])
plt.ylim([-0.8,0.8])
plt.text(2.3,0.98, '(d)', va='top', ha ='left', fontsize=16)
#mag_res = ml_array - mw_array

# make colourbar
plt.gcf().subplots_adjust(right=0.85)
cax = fig.add_axes([0.87,0.2,0.02,0.6]) # setup colorbar axes.
norm = colors.Normalize(vmin=-0.7, vmax=1.7)
cb = colorbar.ColorbarBase(cax, cmap=cm, norm=norm, orientation='vertical')

logticks = cb.get_ticks()[0:]
labels = [str('%0.1f' % 10**x) for x in logticks]
labels[0] = str('%0.2f' % 10**logticks[0])
labels[1] = str('%0.2f' % 10**logticks[1])
labels[2] = str('%0.2f' % 10**logticks[2])
labels[3] = str('%0.2f' % 10**logticks[3])
labels[4] = str('%0.2f' % 10**logticks[4])

cb.set_ticks(logticks)
cb.set_ticklabels(labels)

cb.set_label('Brune Stress Drop (MPa)', rotation=270, fontsize=16, labelpad=16)

# now show figs 
plt.savefig('mw_res_zone.png', fmt='png', dpi=300, bbox_inches='tight')       
plt.show()    

##########################################################################################
# plot regional residual with ML
mw_pred = c0 * ml_array**2 + c1 * ml_array + c2 # NSHA23
mw_res = mw_array - mw_pred
xloc = 2.87
yloc = 0.736        
xloc2 = 6.63
yloc2 = -0.736

fig = plt.figure(10, figsize=(14, 8.))

#plt all data
plt.subplot(221)
plt.plot([2,7], [0,0], 'k--', lw=0.5)
sc = plt.scatter(ml_array, mw_res, c=log10(sd_array), vmin=-0.7, vmax=1.7, s=36, cmap=cm, label='Data')
#sc = plt.scatter(ml_array, mw_res, c=log10(sd_array), vmin=-0.7, vmax=1.7, s=36, cmap=cm, label='Data')
plt.text(xloc, yloc, 'All Data', va='top', ha ='left', fontsize=12, bbox=props)
plt.grid(which='both')
plt.ylabel('$\mathregular{M_{Brune}}$ - $\mathregular{M_{W(Conv)}}$', fontsize=16)
plt.xlim([2.75,7])
plt.ylim([-0.8,0.8])
plt.text(2.3,0.98, '(a)', va='top', ha ='left', fontsize=16)

# get mean & std
mean_res = nanmean(mw_res)
med_res = nanmedian(mw_res[idx])
std_res = nanstd(mw_res)
stat_txt = r'$\mu$ = '+str('%0.2f' % mean_res) + '\n' \
           + r'$\tilde{x}$ = '+str('%0.2f' % med_res) + '\n' \
           + r'$\tau$ = '+str('%0.2f' % std_res)
plt.text(xloc2, yloc2, stat_txt, va='bottom', ha ='right', fontsize=12, bbox=props)


#plt EA
plt.subplot(222)
plt.plot([2,7], [0,0], 'k--', lw=0.5)
idx = where(mzone_array == 'EA')[0]
sc = plt.scatter(ml_array[idx], mw_res[idx], c=log10(sd_array[idx]), vmin=-0.7, vmax=1.7, s=36, cmap=cm, label='Data')
#sc = plt.scatter(ml_array[idx], mw_res[idx], c=log10(sd_array[idx]), vmin=-0.7, vmax=1.7, s=36, cmap=cm, label='Data')
plt.text(xloc, yloc, 'Eastern Australia', va='top', ha ='left', fontsize=12, bbox=props)
#plt.xlabel('$\mathregular{M_{L(MLM92)}}$', fontsize=18)
#plt.ylabel('Brune $\mathregular{M_W}$', fontsize=18)

# get mean & std
mean_res = nanmean(mw_res[idx])
med_res = nanmedian(mw_res[idx])
std_res = nanstd(mw_res[idx])
stat_txt = r'$\mu$ = '+str('%0.2f' % mean_res) + '\n' \
           + r'$\tilde{x}$ = '+str('%0.2f' % med_res) + '\n' \
           + r'$\tau$ = '+str('%0.2f' % std_res)
plt.text(xloc2, yloc2, stat_txt, va='bottom', ha ='right', fontsize=12, bbox=props)

plt.grid(which='both')
plt.xlim([2.75,7])
plt.ylim([-0.8,0.8])
plt.text(2.3,0.98, '(b)', va='top', ha ='left', fontsize=16)

#plt WCA
plt.subplot(223)
plt.plot([2,7], [0,0], 'k--', lw=0.5)
idx = where(mzone_array == 'WCA')[0]
sc = plt.scatter(ml_array[idx], mw_res[idx], c=log10(sd_array[idx]), vmin=-0.7, vmax=1.7, s=36, cmap=cm, label='Data')
#sc = plt.scatter(ml_array[idx], mw_res[idx], c=log10(sd_array[idx]), vmin=-0.7, vmax=1.7, s=36, cmap=cm, label='Data')
plt.text(xloc, yloc, 'Western and Central Australia', va='top', ha ='left', fontsize=12, bbox=props)

# get mean & std
mean_res = nanmean(mw_res[idx])
med_res = nanmedian(mw_res[idx])
std_res = nanstd(mw_res[idx])
stat_txt = r'$\mu$ = '+str('%0.2f' % mean_res) + '\n' \
           + r'$\tilde{x}$ = '+str('%0.2f' % med_res) + '\n' \
           + r'$\tau$ = '+str('%0.2f' % std_res)
plt.text(xloc2, yloc2, stat_txt, va='bottom', ha ='right', fontsize=12, bbox=props)

plt.xlabel('$\mathregular{M_{L}}$', fontsize=16)
plt.ylabel('$\mathregular{M_{Brune}}$ - $\mathregular{M_{W(Conv)}}$', fontsize=16)
plt.grid(which='both')
plt.xlim([2.75,7])
plt.ylim([-0.8,0.8])
plt.text(2.3,0.98, '(c)', va='top', ha ='left', fontsize=16)

#plt SA
plt.subplot(224)
plt.plot([2,7], [0,0], 'k--', lw=0.5)
idx = where(mzone_array == 'SA')[0]
sc = plt.scatter(ml_array[idx], mw_res[idx], c=log10(sd_array[idx]), vmin=-0.7, vmax=1.7, s=36, cmap=cm, label='Data')
plt.text(xloc, yloc, 'South Australia', va='top', ha ='left', fontsize=12, bbox=props)

# get mean & std
mean_res = nanmean(mw_res[idx])
med_res = nanmedian(mw_res[idx])
std_res = nanstd(mw_res[idx])
stat_txt = r'$\mu$ = '+str('%0.2f' % mean_res) + '\n' \
           + r'$\tilde{x}$ = '+str('%0.2f' % med_res) + '\n' \
           + r'$\tau$ = '+str('%0.2f' % std_res)
plt.text(xloc2, yloc2, stat_txt, va='bottom', ha ='right', fontsize=12, bbox=props)

plt.xlabel('$\mathregular{M_{L}}$', fontsize=16)
#plt.ylabel('$\mathregular{M_{Brune}}$', fontsize=16)#plt.ylabel('Brune $\mathregular{M_W}$', fontsize=18)
plt.grid(which='both')
plt.xlim([2.75,7])
plt.ylim([-0.8,0.8])
plt.text(2.3,0.98, '(d)', va='top', ha ='left', fontsize=16)
#mag_res = ml_array - mw_array

# make colourbar
plt.gcf().subplots_adjust(right=0.85)
cax = fig.add_axes([0.87,0.2,0.02,0.6]) # setup colorbar axes.
norm = colors.Normalize(vmin=-0.7, vmax=1.7)
cb = colorbar.ColorbarBase(cax, cmap=cm, norm=norm, orientation='vertical')

logticks = cb.get_ticks()[0:]
labels = [str('%0.1f' % 10**x) for x in logticks]
labels[0] = str('%0.2f' % 10**logticks[0])
labels[1] = str('%0.2f' % 10**logticks[1])
labels[2] = str('%0.2f' % 10**logticks[2])
labels[3] = str('%0.2f' % 10**logticks[3])
labels[4] = str('%0.2f' % 10**logticks[4])

cb.set_ticks(logticks)
cb.set_ticklabels(labels)

cb.set_label('Brune Stress Drop (MPa)', rotation=270, fontsize=16, labelpad=16)

# now show figs 
plt.savefig('mw_res_zone_with_ml.png', fmt='png', dpi=300, bbox_inches='tight')       
plt.show()    


##########################################################################################
# plot ML-MW diff with SD
fig = plt.figure(4, figsize=(10, 8))
m_res = ml_array - mw_array


plt.plot([.1, 200], [0,0], 'k--')
cm = plt.cm.get_cmap('viridis')
#cm = cmc.batlow
sc = plt.scatter(sd_array, m_res, c=mw_array, vmin=3, vmax=6.6, s=36, cmap=cm, label='Data')
plt.xlabel(r"$\Delta\sigma$ (MPa)", fontsize=18)
plt.ylabel('$\mathregular{M_{L(2800)}}$ - $\mathregular{M_{Brune}}$', fontsize=18)
plt.xscale("log")  
plt.grid(which='both') 
plt.xlim([0.05, 100])

# make colourbar
plt.gcf().subplots_adjust(right=0.85)
cax = fig.add_axes([0.87,0.2,0.02,0.6]) # setup colorbar axes.
norm = colors.Normalize(vmin=3, vmax=6.6)
cb = colorbar.ColorbarBase(cax, cmap=cm, norm=norm, orientation='vertical')
cb.set_label('$\mathregular{M_{Brune}}$', rotation=270, fontsize=18, labelpad=20)

# now show figs 
plt.savefig('mw_ml_res_vs_stressdrop.png', fmt='png', dpi=300, bbox_inches='tight')       
plt.show()    

##########################################################################################
# plot ML-MW diff with SD std
fig = plt.figure(10, figsize=(10, 8))
m_res = ml_array - mw_array


plt.plot([.05, 200], [0,0], 'k--')
cm = plt.cm.get_cmap('plasma')
#cm = cmc.batlow
sc = plt.scatter(sd_array, m_res, c=sd_std_array, vmin=0.1, vmax=0.9, s=36, cmap=cm, label='Data')
plt.xlabel(r"$\Delta\sigma$ (MPa)", fontsize=18)
plt.ylabel('$\mathregular{M_{L(2800)}}$ - $\mathregular{M_{Brune}}$', fontsize=18)
plt.xscale("log")  
plt.grid(which='both') 
plt.xlim([0.05, 100])

# make colourbar
plt.gcf().subplots_adjust(right=0.85)
cax = fig.add_axes([0.87,0.2,0.02,0.6]) # setup colorbar axes.
norm = colors.Normalize(vmin=0.1, vmax=0.9)
cb = colorbar.ColorbarBase(cax, cmap=cm, norm=norm, orientation='vertical')
cb.set_label(r"log $\Delta\sigma$ STD (MPa)", rotation=270, fontsize=18, labelpad=20)

# now show figs 
plt.savefig('mw_ml_res_vs_stressdrop_std.png', fmt='png', dpi=300, bbox_inches='tight')       
plt.show()    

##########################################################################################
# plot SD vs mag
fig = plt.figure(5, figsize=(14, 4))

plt.subplot(121)
#plt.plot([.1, 200], [0,0], 'k--')  '#1f77b4', '#ff7f0e'
#sc = plt.scatter(sd_array, m_res, c=mw_array, vmin=3, vmax=6.5, s=36, cmap=cm, label='Data')
idx = where(isnan(sd_array) == False)[0]
c = linregress(mw_array[idx], log10(sd_array[idx]))
mplt = array([3, 7])
splt = c[0]*mplt + c[1]

surf_rup_sd = dictlist2array(mwdat, 'sd')
surf_rup_mag = dictlist2array(mwdat, 'mw')
surf_rup_flag = dictlist2array(mwdat, 'surf_rup')
sidx = where(array(surf_rup_flag)==1)[0]

print('Correlation Coef: '+str(c[2]))
print('log mean: '+str(nanmean(log10(sd_array))))
print('log std: '+str(nanstd(log10(sd_array))))
print('Correlation Coef: '+str(c[2]))

plt.semilogy(mw_array, sd_array, 'o', c='0.8', ms=6.5, label='Blind Ruptures')
plt.semilogy(surf_rup_mag[sidx], surf_rup_sd[sidx], 'o', c='0.8', ms=6.5, mec='r', mew=0.75, label='Surface Ruptures')
plt.semilogy(mplt, 10**splt, '--', c='k', lw=3)
plt.ylabel(r"$\Delta\sigma$ (MPa)", fontsize=18)
#plt.xlabel('$\mathregular{M_{Brune}}$', fontsize=18)
plt.xlabel('M', fontsize=18, weight="bold")
plt.grid(which='both') 
plt.ylim([0.05, 100])
plt.legend(loc=4, numpoints=1, fontsize=10)

# annotate
xloc2 = 6.92
yloc2 = get_log_xy_locs([0.05, 100], 0.21)

stat_txt = r'$\mu$ = '+str('%0.2f' % 10**nanmean(log10(sd_array[idx]))) + '\n' \
           + r'$\tilde{x}$ = '+str('%0.2f' % nanmedian(sd_array[idx])) + '\n' \
           + r'log $\tau$ = '+str('%0.2f' % nanstd(log10(sd_array[idx]))) + '\n' \
           + '$r^{2}$ = '+str('%0.3f' % c.rvalue**2)
plt.text(xloc2, yloc2, stat_txt, va='bottom', ha ='right', fontsize=12, bbox=props)
plt.xlim([3, 7])

'''
# make colourbar
plt.gcf().subplots_adjust(right=0.85)
cax = fig.add_axes([0.87,0.2,0.02,0.6]) # setup colorbar axes.
norm = colors.Normalize(vmin=3, vmax=6.5)
cb = colorbar.ColorbarBase(cax, cmap=cm, norm=norm, orientation='vertical')
cb.set_label('$\mathregular{M_{Brune}}$', rotation=270, fontsize=18, labelpad=20)
'''

# now plot hist
ax = plt.subplot(143)
bins = arange(-1.3, 2.5, 0.2)
plt.hist(sd_array[idx], bins=10**bins, orientation='horizontal', facecolor='0.8') #, width=0.8)
#plt.yticks([])
#ax.set_yticklabels([])
ax.tick_params(labelleft=False) 
plt.yscale('log')
plt.ylim([0.05, 100])
plt.xlabel('Count', fontsize=18)
plt.tight_layout() 
#plt.subplots_adjust(wspace=0.07)

# now show figs 
plt.savefig('stressdrop_vs_mag.png', fmt='png', dpi=300, bbox_inches='tight')       
plt.show() 

##########################################################################################
# plot SD vs ML and ML
fig = plt.figure(5, figsize=(14, 4))

plt.subplot(121)
#plt.plot([.1, 200], [0,0], 'k--')  '#1f77b4', '#ff7f0e'
#sc = plt.scatter(sd_array, m_res, c=mw_array, vmin=3, vmax=6.5, s=36, cmap=cm, label='Data')
idx = where((isnan(sd_array) == False) & (isnan(ml_array) == False))[0]
c = linregress(ml_array[idx], log10(sd_array[idx]))
mplt = array([3, 7])
splt = c[0]*mplt + c[1]

print('Correlation Coef: '+str(c[2]))
print('log mean: '+str(nanmean(log10(sd_array))))
print('log std: '+str(nanstd(log10(sd_array))))
print('Correlation Coef: '+str(c[2]))

plt.semilogy(ml_array, sd_array, 'o', c='0.8', ms=6.5)
plt.semilogy(mplt, 10**splt, '--', c='k', lw=3)
plt.ylabel(r"$\Delta\sigma$ (MPa)", fontsize=18)
plt.xlabel('$\mathregular{M_{L(2800)}}$', fontsize=18)
plt.grid(which='both') 
plt.ylim([0.05, 100])
#plt.legend(loc=4, numpoints=1, fontsize=10)

# annotate
xloc2 = 6.9
yloc2 = get_log_xy_locs([0.05, 100], 0.04)

stat_txt = r'$r^{2}$ = '+str('%0.3f' % c.rvalue**2)
plt.text(xloc2, yloc2, stat_txt, va='bottom', ha ='right', fontsize=14, bbox=props)
plt.xlim([3, 7])


# now plot hist
ax = plt.subplot(122)
idx = where(isnan(sd_array) == False)[0]
c = linregress(mw_array[idx], log10(sd_array[idx]))
mplt = array([3, 7])
splt = c[0]*mplt + c[1]

print('Correlation Coef: '+str(c[2]))
print('log mean: '+str(nanmean(log10(sd_array))))
print('log std: '+str(nanstd(log10(sd_array))))
print('Correlation Coef: '+str(c[2]))

plt.semilogy(mw_array, sd_array, 'o', c='0.8', ms=6.5)
plt.semilogy(mplt, 10**splt, '--', c='k', lw=3)
#plt.ylabel(r"$\Delta\sigma$ (MPa)", fontsize=18)
plt.xlabel('M', fontsize=18, weight="bold")
plt.grid(which='both') 
plt.ylim([0.05, 100])
#plt.legend(loc=4, numpoints=1, fontsize=10)

# annotate
stat_txt = r'$r^{2}$ = '+str('%0.3f' % c.rvalue**2)
plt.text(xloc2, yloc2, stat_txt, va='bottom', ha ='right', fontsize=14, bbox=props)
plt.xlim([3, 7])
plt.tight_layout() 

# now show figs 
plt.savefig('stressdrop_vs_ml_mw.png', fmt='png', dpi=300, bbox_inches='tight')       
plt.show() 