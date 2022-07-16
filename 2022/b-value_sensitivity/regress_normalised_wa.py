from numpy import array, nan, unique, where, log10, median, nanmean, nanmedian, nanstd, polyfit, savetxt, loadtxt, random, \
                  isnan, mean, hstack, logspace, linspace, arange, sqrt, exp, zeros_like, ones_like, unique, polyval
from scipy.stats import linregress
import matplotlib.pyplot as plt
from misc_tools import get_binned_stats, get_binned_stats_mean, remove_last_cmap_colour
from scipy.odr import Data, Model, ODR, models
import scipy.odr.odrpack as odrpack
import matplotlib as mpl
from misc_tools import dictlist2array, get_mpl2_colourlist, remove_last_cmap_colour
from gmt_tools import cpt2colormap 
from os import getcwd

mpl.style.use('classic')

print('change stress drop')

if getcwd().startswith('/nas'):
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/Paired_10.cpt'
else:
    cptfile = '//Users//trev//Documents//DATA//GMT//cpt//Paired_10.cpt'
    
ncolours = 11
cmap, zvals = cpt2colormap(cptfile, ncolours)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncolours-1)))


###############################################################################
# load simulated W-A data
###############################################################################

wafile = 'smsim/simulated_wa_amplitudes.csv'

data = loadtxt(wafile, delimiter=',', skiprows=1)

events = data[:,0]
mags = data[:,1]
dists = data[:,2]
wa_amps = data[:,3]

'''	
# add random Gaussian to WA amps - this is now being done in the simulate code
rand_gaus = random.normal(loc=0., scale=0.2, size=len(wa_amps))
wa_amps = 10**(log10(wa_amps) + rand_gaus)
'''
	
uevents = unique(events)

###############################################################################
# get median near-source attenuation
###############################################################################
afit = []
umw = []
for ue in uevents:
    idx = where((events == ue) & (dists <= 100.))[0]
    if len(idx) >= 3:
        afit.append(linregress(log10(dists[idx]), log10(wa_amps[idx]))[0])
    
    # get unique mws
    idx = where(events == ue)[0]
    umw.append(mags[idx][0])

afit = array(afit)
umw = array(umw)

idx = where(~isnan(afit))[0]
plt.hist(afit[idx],arange(0., 3.05, 0.1))
plt.xlabel('Intra-Event Geometric Attenuation | Rhyp <= 100 km', fontsize=16)
plt.ylabel('Count', fontsize=16)
 
# get average
alpha = abs(nanmean(afit))  
alpha = abs(nanmedian(afit))
plt.text(-2.75, 14.1, 'G(R) = '+str('%0.2f' % alpha), va='center', fontsize=16) 
#plt.xlim([-3, 1.5])
#plt.ylim([0, 15])
plt.savefig('av_geo_spread_hist.png', fmt='png', bbox_inches='tight')

#plt.show()

###############################################################################
# get normalised data
###############################################################################
lognormamp = []
rhyp = []
for ue in uevents:
    # get src amp
    idx100 = where((dists <= 100.) & (events == ue))[0]
    if len(idx100) >= 1:
        srcamp = mean(alpha * log10(dists[idx100]) + log10(wa_amps[idx100]))
        
        # get norm amps
        idx = where(events == ue)[0]
        tempamps = log10(wa_amps[idx]) - srcamp
        lognormamp = hstack((lognormamp, tempamps))
        
        # get distances
        rhyp = hstack((rhyp, dists[idx]))
        
###############################################################################
# plot data
###############################################################################
        
# now window data at 17 km and get scaling factor
log10rhyp = log10(rhyp)
log17 = log10(17.) # from HB 1987
idx = where((log10rhyp >= log17-0.05) & (log10rhyp <= log17+0.05))[0]

# 100 km - override
idx = where((log10rhyp >= 1.95) & (log10rhyp <= 2.05))[0]
scale_fact = median(lognormamp[idx])

# adjust with scale factor
lognormamp -= scale_fact
#lognormamp -= scale_fact-1 #17 km

plt.figure(2, figsize=(8,8))
plt.plot(log10(rhyp), lognormamp, '+', color='0.6')
plt.xlim([log10(3), 3])
plt.ylim([-1.5, 2.75])

# bin data
bins = arange(0., 2.8, 0.1)
medbin, stdbin, meanx, outbins, nperbin = get_binned_stats(bins, log10(rhyp), lognormamp)
idx = where(~isnan(medbin))[0]
plt.errorbar(meanx, medbin, yerr=stdbin, color='r', fmt='s', mec='r')

plt.xlabel('log Hypocentral Distance (km)', fontsize=16)
plt.ylabel('Normalised log Wood-Anderson (mm)', fontsize=16)
plt.grid(which='major')

###############################################################################
# fit linear WA attenuation
###############################################################################

# fit linear
lin_fit = linregress(meanx, medbin) # regress binned data
#lin_fit = linregress(log10(rhyp), lognormamp) # regress all data
xfit = array([0.5, 2.85])
yfit = lin_fit[0]*xfit + lin_fit[1]
#plt.plot(xfit, yfit, 'k-', lw=2.0)

###############################################################################
# fit tri-linear WA attenuation
###############################################################################
# fit tri-linear
R1 = 70.
R2 = 150.

meanx = log10(rhyp) # use all data
medbin = lognormamp

# get dists LE 90
didx = where(10**meanx <= R1)[0]
tri1 = linregress(meanx[didx], medbin[didx])
#tri1 = linregress(log10(rhyp[didx]), lognormamp[didx]) # regress all data

b1 = tri1[0]
c0 = tri1[1]

### need to fit each segment ###

def fit_trilinear(c, x):
    from numpy import zeros_like, where
    # set hinge distances (from Allen, 2012)
    #R1 = 90.
    #R2 = 150.
    
    ans = zeros_like(x)
    
    idx1 = where(x <= log10(R1))[0]
    idx2 = where((x > log10(R1)) & (x <= log10(R2)))[0]
    idx3 = where(x > log10(R2))[0]
        
    ans[idx1] = c0 + b1 * x[idx1]
    
    ans[idx2] = c0 + b1 * log10(R1) \
                + c[0] * log10(10**x[idx2] / R1)
    
    ans[idx3] = c0 + b1 * log10(R1) \
                + c[0] * log10(R2 / R1) \
                + c[1] * log10(10**x[idx3] / R2)
    
    return ans

'''
data = odrpack.RealData(meanx, medbin)
tri = odrpack.Model(fit_trilinear)
odr = odrpack.ODR(data, tri, beta0=[-1.0, 0.0, -2.0])
odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq
out = odr.run()
out.pprint()
'''

def fit_b2(c, x):
    # c0 defined above
    
    from numpy import zeros_like, where
    
    ans = zeros_like(x)
    
    idx1 = where(x <= log10(R1))[0]
    idx2 = where((x > log10(R1)) & (x <= log10(R2)))[0]
        
    ans[idx1] = c0 + b1 * x[idx1]
    
    ans[idx2] = c0 + b1 * log10(R1) \
                + c[0] * log10(10**x[idx2] / R1)
 
    return ans

# only fit data LT R2
didx = where(10**meanx <= R2)[0]

data = odrpack.RealData(meanx[didx], medbin[didx]) # regress binned data
#data = odrpack.RealData(log10(rhyp[didx]), lognormamp[didx]) # regress all data

tri2 = odrpack.Model(fit_b2)
odr = odrpack.ODR(data, tri2, beta0=[0.0])
odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq
out = odr.run()

b2 = out.beta[0] # 90-150 km
#b3 = out.beta[1] # GT 150 km

def fit_b3(c, x):
    from numpy import zeros_like, where
    # c0 defined above
    
    ans = zeros_like(x)
    
    idx1 = where(x <= log10(R1))[0]
    idx2 = where((x > log10(R1)) & (x <= log10(R2)))[0]
    idx3 = where(x > log10(R2))[0]
        
    ans[idx1] = c0 + b1 * x[idx1]
    
    ans[idx2] = c0 + b1 * log10(R1) \
              + b2 * log10(10**x[idx2] / R1)
    
    ans[idx3] = c0 + b1 * log10(R1) \
              + b2 * log10(R2 / R1) \
              + c[0] * log10(10**x[idx3] / R2)
    
    return ans

# fit all data 
data = odrpack.RealData(meanx, medbin) # regress binned data
#data = odrpack.RealData(log10(rhyp), lognormamp) # regress all data

tri3 = odrpack.Model(fit_b3)
odr = odrpack.ODR(data, tri3, beta0=[-1.7])
odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq
out = odr.run()

b3 = out.beta[0] # 90-150 km

###############################################################################
# plot tri-linear WA attenuation
###############################################################################

# plot model
x1 = log10([5, R1])
y1 = b1*x1 + c0
plt.plot(x1, y1, '-', color='limegreen', lw=2.0)

x2 = log10([R1, R2])
y2 = c0 + b1 * log10(R1) \
     + b2 * log10(10**x2 / R1)
plt.plot(x2, y2, '-', color='limegreen', lw=2.0)

x3 = log10([R2, 800])
y3 = c0 + b1 * log10(R1) \
     + b2 * log10(R2 / R1) \
     + b3 * log10(10**x3 / R2)
plt.plot(x3, y3, '-', color='limegreen', lw=2.0)

###############################################################################
# fit standard ML eqn
###############################################################################

def fit_ml_curve(c, x):
    '''
    –log A0 = n log(r/100) + K(r – 100) + 3.0
    '''
    
    from numpy import sqrt, log10
    
    ans = c[0] * log10(10**x/100.) + c[1]*(10**x - 100) + c[2]
    #ans = c[0] * log10(10**x/17.) + c[1]*(10**x - 17) + c[2]
    
    return ans
    
def fit_ml_curve_3(c, x):
    '''
    –log A0 = n log(r/100) + K(r – 100) + 3.0
    '''
    
    from numpy import sqrt, log10
    
    ans = c[0] * log10(10**x/100.) + c[1]*(10**x - 100.)
    #ans = c[0] * log10(10**x/17.) + c[1]*(10**x - 17) + c[2]
    
    return ans

# fit all data
data = odrpack.RealData(meanx, medbin)

''' fit standard ML eqn with variable DC'''
afit = odrpack.Model(fit_ml_curve)
odr = odrpack.ODR(data, afit, beta0=[1.3, 0.0005, 1.0])

odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
out = odr.run()
c = out.beta

# plot fitted curve
'''
xplt = logspace(0, log10(800), 100)
yplt = c[0] * log10(xplt/100.) + c[1]*(xplt - 100.) + c[2]
#yplt = c[0] * log10(xplt/17.) + c[1]*(xplt - 17.) + c[2]
plt.plot(log10(xplt), yplt, 'r-', lw=2.0)
'''
''' fit standard ML eqn with fixed DC'''
afit = odrpack.Model(fit_ml_curve_3)
odr = odrpack.ODR(data, afit, beta0=[1.3, 0.0005])

odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
out = odr.run()
c = out.beta

# plot fitted curve
xplt = logspace(0, log10(800), 100)
yplt = c[0] * log10(xplt/100.) + c[1]*(xplt - 100.)
#yplt = c[0] * log10(xplt/17.) + c[1]*(xplt - 17.) + 2.
plt.plot(log10(xplt), yplt, 'r-', lw=2.0)

###############################################################################
# plt HB87 reference
###############################################################################

# plot HB84
hb87_curve = -1*(log10(xplt) + 0.00301 * xplt + 0.7) + 3
plt.plot(log10(xplt), hb87_curve, '-', c='dodgerblue', lw=2.0)

#plt.plot([0.0, log17],[1., 1.], 'k--', lw=1.5)
#plt.plot([log17, log17],[-1.5, 1.], 'k--', lw=1.5)

plt.plot([0.0, 2],[0., 0.], 'k--', lw=1.5)
plt.plot([2, 2],[-1.5, 0.], 'k--', lw=1.5)

plt.savefig('norm_sim_logamps.png', bbox_inches='tight')

###############################################################################
# calc SED magnitude
###############################################################################
# ref: Appendix H ECOS-09: Swiss instrumental local magnitudes - Nicholas Deichmann, SED, 2009/12/17
def calc_SED84(comp, logA, rhyp):
    '''
    Assumes W-A amplification factor of 2800
    '''    
    from numpy import nan, log10
    
    # convert logA assuming W-A amplification of 2800
    logA2800 = log10(2800. * (10**logA) / 2080.)
    
    # pre-allocate value
    SED84 = nan
    Ce = 0.1
    if comp == 0: # assume vertical only
        if rhyp <= 60.:
            Cd = 0.0180 * rhyp + 1.77
        else:
            Cd = 0.0038 * rhyp + 2.62
        
        SED84 = logA2800 + Cd + Ce
            
    magstr = 'SED84:\t' + str("%0.2f" % SED84)

    return SED84

###############################################################################
# now calculate MLs using new functions
###############################################################################

fig = plt.figure(3, figsize=(7,14))
ax = plt.subplot(2,1,1)
tri_ml = []
tri_std = []
com_ml = []
com_std = []
mlm92_ml = []
bj84_ml = []
w20_ml = []
sed84_ml = []

rhyp = dists # need to do this as one um resulted in nan afit

ml_fit = arange(1., 6.05, 0.05)

#tri_c = 2 + b1*log10(17)
tri_c = 0
#com_c = 2 + (c[0] * log10(17./17.) + c[1]*(17 - 17.) + c[2])

com_c = c[0] * log10(100./100.) + c[1]*(100 - 100.) + 3.
for ue in uevents:
    # get src amp
    idx = where(events == ue)[0]
    
    # get tri-linear
    tmplogAr = []
    for r in rhyp[idx]:
        if r <= R1:
            logAr = b1*log10(r)
        elif r > R1 and r <= R2:
            logAr = b1*log10(R1) + b2*log10(r/R1)
        elif r > R2:
            logAr = b1*log10(R1) + b2*log10(R2/R1) + b3*log10(r/R2)
        
        tmplogAr.append(logAr)
    
    logA100 = b1*log10(R1) + b2*log10(100./R1)
    tri_c = logA100 + 3.0
    
    tmplogAr = array(tmplogAr)
    tmp_tri_ml = nanmedian(log10(wa_amps[idx]) - tmplogAr + tri_c) + 0.18
    tri_ml.append(tmp_tri_ml)
    tri_std.append(nanstd(log10(wa_amps[idx]) - tmplogAr + tri_c))
    
    
    # get standard functional form
    #logAr = c[0] * log10(rhyp[idx]/17.) + c[1]*(rhyp[idx] - 17.) + c[2]
    logAr = c[0] * log10(rhyp[idx]/100.) + c[1]*(rhyp[idx] - 100.) # fixed DC
    tmp_com_ml = nanmedian(log10(wa_amps[idx]) - logAr + 3.0) + 0.18 # when V amp == 1 mm, actually ML 3.18 on H
    com_ml.append(tmp_com_ml)
    com_std.append(nanstd(log10(wa_amps[idx]) - logAr + 3. + 0.18))
    
    # calc using BJ84
    #log10(rhyp) + 0.00301 * rhyp + 0.70
    bj84_ml.append(nanmedian(log10(wa_amps[idx]) + log10(rhyp[idx]) + 0.00301 * rhyp[idx] + 0.7))
    
    # calc using MLM92
    mlm92_ml.append(nanmedian(log10(wa_amps[idx]) + 1.34 * log10(rhyp[idx] / 100.) + 0.00055 * (rhyp[idx] - 100) + 3.13))
    
    # calc Walter et al (2020)
    didx = where(rhyp[idx] <= 160)[0]
    w20_ml.append(nanmedian(log10(wa_amps[idx][didx]) + 2.01 * log10(rhyp[idx][didx]) - 0.0057 * rhyp[idx][didx] - 0.45) + 0.18) # add 0.18 to covert V to H
    
    # get SED84
    Ce = 0.1 # BB constant?
    logA2800 = log10(2800. * (wa_amps[idx]) / 2080.)
    
    Cd = 0.0180 * rhyp[idx] + 1.77
    didx = rhyp[idx] > 60
    Cd[didx] = 0.0038 * rhyp[idx][didx] + 2.62
    
    sed84_ml.append(nanmedian(logA2800 + Cd + Ce) + 0.18) # add 0.18 to cov
    
    
tri_p = polyfit(tri_ml, umw, 2)
tri_fit = polyval(tri_p, ml_fit)
com_p = polyfit(com_ml, umw, 2)
com_fit = polyval(com_p, ml_fit)  
bj84_p = polyfit(bj84_ml, umw, 2)
bj84_fit = polyval(bj84_p, ml_fit) 
mlm92_p = polyfit(mlm92_ml, umw, 2)
mlm92_fit = polyval(mlm92_p, ml_fit) 
w20_p = polyfit(w20_ml, umw, 2)
w20_fit = polyval(w20_p, ml_fit)
sed84_p = polyfit(sed84_ml, umw, 2)
sed84_fit = polyval(sed84_p, ml_fit) 


print('!!!!! add swiss & OK models !!!!!')

plt.plot([1, 6], [1, 6], 'k--')
plt.plot(tri_ml, umw, 'o', c=cs[0])
plt.plot(com_ml, umw, 'o',  c=cs[2])
plt.plot(bj84_ml, umw, 'o',  c=cs[4])
plt.plot(mlm92_ml, umw, 'o',  c=cs[6])
'''
plt.plot(w20_ml, umw, 'o',  c=cs[8])
plt.plot(sed84_ml, umw, 'o',  c='0.75')
'''
plt.plot(ml_fit, tri_fit, '-', lw=2.0, c=cs[1], label='trilinear')
plt.plot(ml_fit, com_fit, '-', lw=2.0, c=cs[3], label='standard')
plt.plot(ml_fit, bj84_fit, '-', lw=2.0, c=cs[5], label='BJ84')
plt.plot(ml_fit, mlm92_fit, '-', lw=2.0, c=cs[7], label='MLM92')
'''
plt.plot(ml_fit, w20_fit, '-', lw=2.0, c=cs[9], label='Wea20')
plt.plot(ml_fit, sed84_fit, '-', lw=2.0, c='k', label='SED84')
'''
plt.legend(loc=2, numpoints=1)
plt.xlim([1, 6])
plt.ylim([1, 6])
plt.grid(which='both')
plt.ylabel('MW')

######################################################################################
# plt difference
ax = plt.subplot(4,1,3)
plt.plot(ml_fit, tri_fit-ml_fit, '-', lw=2.0, c=cs[1], label='trilinear')
plt.plot(ml_fit, com_fit-ml_fit, '-', lw=2.0, c=cs[3], label='standard')
plt.plot(ml_fit, bj84_fit-ml_fit, '-', lw=2.0, c=cs[5], label='BJ84')
plt.plot(ml_fit, mlm92_fit-ml_fit, '-', lw=2.0, c=cs[7], label='MLM92')
plt.xlabel('ML')
plt.grid(which='both')


plt.savefig('simulated_ml_mw_relationships.png', fmt='png', dpi=300, 	bbox_inches='tight')
plt.show()