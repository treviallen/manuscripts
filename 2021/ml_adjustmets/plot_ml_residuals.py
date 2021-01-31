# Allen 2016 (unpub) - SE Aust
def calc_A16(logA, rhyp, comp):
    from numpy import loadtxt, log10
    coeffile = 'A16_coeffs.txt'
    dat = loadtxt(coeffile)
    b1 = dat[0]
    b2 = dat[1]
    b3 = dat[2]
    r1 = dat[3]
    r2 = dat[4]
    c0 = dat[5] - 3.0 # norm 
    c1 = dat[6] # h/v correction
    c1l = dat[7]
    c2l = dat[8]
    
    if rhyp <= r1:
        logASEA = b1*log10(rhyp) + c0
    elif rhyp > r1 and rhyp <= r2:
        logASEA = b1*log10(r1) + b2*log10(rhyp/r1) + c0
    elif rhyp > r2:
        logASEA = b1*log10(r1) + b2*log10(r2/r1) + b3*log10(rhyp/r2) + c0

    A16 = logA - logASEA + c1 # constant H/V
    A16 = logA - logASEA + (c1l * log10(rhyp) + c2l) # linear H/V
    
    #print('A16:\t' + str("%0.1f" % A16), logA - logASEA, dat[5],dat[6]
    
    return A16

from datetime import datetime
from numpy import array, arange, hstack, log10, unique, nanmean, nanstd, isnan, ones_like, where, mean, std, random, argsort
from scipy.stats import trim_mean
from misc_tools import get_binned_stats, get_binned_stats_mean
from calculate_magnitudes import *
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('classic')

# parse pre-existing data to find event time
wavals = []
dtvect = []
hvrat = []
hmax = []
lines = open('/Users/trev/Dropbox/Magnitudes/2016_working/seaust_w-a_amps_post2009.csv').readlines()
for i, line in enumerate(lines):
    if line.startswith('#') == False:
        dat = line.strip().split(',')
        
        # make dt object
        dt = datetime.strptime(dat[0], '%Y-%m-%d %H:%M:%S')
        dtvect.append(dt)
        
        # old data - not logged
        if i <= 1169:
            adict = {'datetime':dt, 'datestr':dat[0], 'lon':dat[1], 'lat':dat[2], \
                     'dep':dat[3], 'mag':dat[4], 'stn':dat[5], 'rhyp':dat[7], 'repi':dat[6], \
                     'logAe':log10(float(dat[8])), 'logAn':log10(float(dat[8])), 'logAz':log10(float(dat[10]))}
            
        # new data logged
        else:
            adict = {'datetime':dt, 'datestr':dat[0], 'lon':dat[1], 'lat':dat[2], \
                     'dep':dat[3], 'mag':dat[4], 'stn':dat[5], 'rhyp':dat[7], 'repi':dat[6], \
                     'logAe':float(dat[8]), 'logAn':float(dat[8]), 'logAz':float(dat[10])}
        
        hmax = max([adict['logAn'], adict['logAe']])
        adict['logAHmax'] = hmax
                         	
        wavals.append(adict)

# get unique events
dtvect = array(dtvect)
udt = unique(dtvect)

# loop thru events
evdict = []

for u in udt:
    tev = {'datetime':u, 'stns':[], 'e':[], 'n':[], 'z':[], 'hmax':[], 'rhyp':[], 'repi':[]}
    stns = []
    e = []
    n = []
    z = []
    hmax = []
    repi = []
    rhyp = []
    
    for waval in wavals:
        if waval['datetime'] == u:
            tev['evlo'] = waval['lon']
            tev['evla'] = waval['lat']
            tev['evdp'] = waval['dep']
            if float(waval['rhyp']) > 500. and float(waval['mag']) < 4.0:
                # ignore data
                a = 1
            elif float(waval['rhyp']) > 250. and float(waval['mag']) < 3.0:
                # ignore data
                a = 1
            else:
                # add data
                stns.append(waval['stn'])
                e.append(waval['logAe'])
                n.append(waval['logAn'])
                z.append(waval['logAz'])
                hmax.append(waval['logAHmax'])
                rhyp.append(float(waval['rhyp']))
                repi.append(float(waval['repi']))
    tev['stns'] = array(stns)
    tev['e'] = array(e)
    tev['n'] = array(n)
    tev['z'] = array(z)
    tev['hmax'] = array(hmax)
    tev['rhyp'] = array(rhyp)
    tev['repi'] = array(repi)

    evdict.append(tev)

###############################################################################
# get mag residuals
###############################################################################

# loop thru events
r35res   = []
bj84res = []
sed84res = []
hb87res  = []
gg91res  = []
mlm92res = []
wgw96res = []
a16res   = []
r35evstd   = []
bj84evstd = []
sed84evstd = []
hb87evstd  = []
gg91evstd  = []
mlm92evstd = []
wgw96evstd = []
a16evstd   = []
rhypres  = []
lores = []
lares = []
seamags  = []

mlm92meanmags = []
a16meanmags = []
r35meanmags = []
wgw96meanmags = []
bj84meanmags = []
sed84meanmags = []
hb87meanmags = []

mlm92sampstd = []
a16sampstd = []
r35sampstd = []
wgw96sampstd = []
bj84sampstd = []
hb87sampstd = []
sed84sampstd = []

nrecs = []
maxDist = []
minDist = []

###############################################################
# function to sample errors
###############################################################

def sample_mag_stds(stnmags, meanmag):
    if meanmag >= 5.0:
        nsamp = 4
    elif meanmag >= 4.0:
        nsamp = 3
    elif meanmag >= 3.0:
        nsamp = 3
    else:
        nsamp = 3
    
    sampled_mags = []
    for i in range(0,1000):
        sampled_mags.append(nanmean(random.choice(stnmags, nsamp)))
    
    return std(sampled_mags)

###############################################################
# loop thru events

# make header
lines = 'DATETIME,STN,RHYP,R35,BJ84,SED84,HB87,GG91,MLM92,WGW96,A16,RES\n'

for ev in evdict:
    # get sort idx
    si=argsort(ev['rhyp'])
    
    # calc R35
    mag = []
    comp = 1
    for amp, dist in zip(ev['hmax'], ev['repi']):
        mag.append(calc_R35(comp, amp, dist))
    r35mag = array(mag)
    
    idx = where(isnan(r35mag) == False)[0]
    magres = array(mag) - trim_mean(r35mag[idx], 0.125)
    r35evstd.append(nanstd(magres))
    r35res = hstack((r35res, magres))
    r35meanmags.append(trim_mean(r35mag[idx], 0.125))
    r35sampstd.append(sample_mag_stds(r35mag[si], trim_mean(r35mag[idx], 0.125)))
    
    # calc JB84
    mag = []
    comp = 1
    for amp, dist in zip(ev['hmax'], ev['repi']):
        mag.append(calc_BJ84(comp, amp, dist))
    bj84mag = array(mag)
    
    idx = where(isnan(bj84mag) == False)[0]
    magres = array(mag) - trim_mean(bj84mag[idx], 0.125)
    bj84evstd.append(nanstd(magres))
    bj84res = hstack((bj84res, magres))
    bj84meanmags.append(trim_mean(bj84mag[idx], 0.125))
    bj84sampstd.append(sample_mag_stds(bj84mag[si], trim_mean(bj84mag[idx], 0.125)))
    
    # calc SED84
    mag = []
    comp = 0
    for amp, dist in zip(ev['z'], ev['rhyp']):
        mag.append(calc_SED84(comp, amp, dist))
    sed84mag = array(mag)
    
    idx = where(isnan(sed84mag) == False)[0]
    magres = array(mag) - trim_mean(sed84mag[idx], 0.125)
    sed84evstd.append(nanstd(magres))
    sed84res = hstack((sed84res, magres))
    sed84meanmags.append(trim_mean(sed84mag[idx], 0.125))
    sed84sampstd.append(sample_mag_stds(sed84mag[si], trim_mean(sed84mag[idx], 0.125)))
    
    # calc HB87
    mag = []
    comp = 1
    for amp, dist in zip(ev['hmax'], ev['rhyp']):
        mag.append(calc_HB87(comp, amp, dist))
    hb87mag = array(mag)
    
    idx = where(isnan(hb87mag) == False)[0]
    magres = array(mag) - trim_mean(hb87mag[idx], 0.125)
    hb87evstd.append(nanstd(magres))
    hb87res = hstack((hb87res, magres))
    
    # calc GG91
    mag = []
    comp = 0
    for amp, dist in zip(ev['z'], ev['rhyp']):
        mag.append(calc_GG91(comp, amp, dist))
    gg91mag = array(mag)
    
    magres = array(mag) - trim_mean(mag, 0.125)
    gg91evstd.append(nanstd(magres))
    gg91res = hstack((gg91res, magres))
    
    # calc MLM92
    mag = []
    comp = 0
    mxDist = 0
    mnDist = 99999
    for amp, dist in zip(ev['z'], ev['rhyp']):
        mag.append(calc_MLM92(comp, amp, dist))
    mlm92mag = array(mag)
    
    magres = array(mag) - trim_mean(mag, 0.125)
    mlm92evstd.append(nanstd(magres))
    mlm92res = hstack((mlm92res, magres))
    mlm92meanmags.append(trim_mean(mag, 0.125))
    mlm92sampstd.append(sample_mag_stds(mlm92mag[si], trim_mean(mag, 0.125)))
    
    # make official mag array
    seamags = hstack((seamags, ones_like(magres)*trim_mean(mag, 0.125)))
    rhypres = hstack((rhypres, ev['rhyp']))
    lores = hstack((lores, ones_like(magres)*float(ev['evlo'])))
    lares = hstack((lares, ones_like(magres)*float(ev['evla'])))
    
    # calc WGW96
    mag = []
    comp = 0
    for amp, dist in zip(ev['z'], ev['rhyp']):
        mag.append(calc_WGW96(comp, amp, dist))
    wgw96mag = array(mag)
    
    magres = array(mag) - trim_mean(mag, 0.125)
    wgw96evstd.append(nanstd(magres))
    wgw96res = hstack((wgw96res, magres))
    wgw96meanmags.append(trim_mean(mag, 0.125))
    wgw96sampstd.append(sample_mag_stds(wgw96mag[si], trim_mean(mag, 0.125)))
    
    # calc A16 - still need to check!
    mag = []
    for amp, dist in zip(ev['z'], ev['rhyp']):
        mag.append(calc_A16(amp, dist, 0))
    a16mag = array(mag)
    
    magres = array(mag) - trim_mean(mag, 0.125)
    a16evstd.append(nanstd(magres))
    a16res = hstack((a16res, magres)) 
    a16meanmags.append(trim_mean(mag, 0.125))
    a16sampstd.append(sample_mag_stds(a16mag[si], trim_mean(mag, 0.125)))
    
    # get stats
    maxDist.append(max(ev['rhyp']))
    minDist.append(min(ev['rhyp']))
    nrecs.append(len(ev['rhyp']))
    
    # make lines
    dt = ev['datetime']
    for i in range (0, len(ev['rhyp'])):
        lines += ','.join((dt.strftime('%Y%m%d%H%M'), ev['stns'][i], str(ev['rhyp'][i]), \
                           str('%0.2f' % r35mag[i]), str('%0.2f' % bj84mag[i]), \
                           str('%0.2f' % sed84mag[i]), str('%0.2f' % hb87mag[i]), \
                           str('%0.2f' % gg91mag[i]), str('%0.2f' % mlm92mag[i]), \
                           str('%0.2f' % wgw96mag[i]), str('%0.2f' % a16mag[i]), \
                           str('%0.2f' % (nanmean(mlm92mag) - mlm92mag[i])))) + '\n'

r35evstd   = array(r35evstd)
hb87evstd  = array(hb87evstd)
sed84evstd  = array(sed84evstd)
gg91evstd  = array(gg91evstd)
mlm92evstd = array(mlm92evstd)
wgw96evstd = array(wgw96evstd)
a16evstd   = array(a16evstd)

r35meanmags = array(r35meanmags)
bj84meanmags = array(bj84meanmags)
sed84meanmags = array(sed84meanmags)
hb87meanmags = array(hb87meanmags)
wgw96meanmags = array(wgw96meanmags)
a16meanmags = array(a16meanmags)
mlm92meanmags = array(mlm92meanmags)

midx = where(mlm92meanmags >= 4.)[0]
print('MLM92 mean std', nanmean(mlm92evstd[midx]))
print('A16   mean std', nanmean(a16evstd[midx]))

# export lines
f = open('sea_stn_mags.csv', 'w')
f.write(lines)
f.close()

###############################################################################
# plot spread of sampled stds
###############################################################################

fig = plt.figure(6, figsize=(6,6))
plt.plot([0,0.25],[0,0.25],'k--')
plt.plot(mlm92sampstd, a16sampstd, 'r+', ms=10, lw=1.5)
plt.xlabel('Sampled MLM92 STDs')
plt.ylabel('Sampled A16 STDs')

#plt.show()

testmag = 3.5
idx = where((array(a16sampstd) > array(mlm92sampstd)) & (a16meanmags > testmag))[0]
midx = where(a16meanmags > testmag)[0]
stdrat = len(idx)*1. / len(a16meanmags[midx])

print('Std Rat:', stdrat)

###############################################################################
# export event info
###############################################################################

lines = 'DATETIME,R35,R35STD,BJ84,BJ84STD,SED84,SED84STD,MLM92,MLM92STD,WGW96,WGW96STD,A16,A16STD,MIN_RHYP,MAX_RHYP,NSTA\n'
for i, ev in enumerate(evdict):
    datestr = datetime.strftime(ev['datetime'],'%Y%m%d%H%M')
    lines += ','.join((datestr, str('%0.2f' % r35meanmags[i]), str('%0.2f' % r35evstd[i]), \
                       str('%0.2f' % bj84meanmags[i]), str('%0.2f' % bj84evstd[i]), \
                       str('%0.2f' % sed84meanmags[i]), str('%0.2f' % sed84evstd[i]), \
                       str('%0.2f' % mlm92meanmags[i]), str('%0.2f' % mlm92evstd[i]), \
                       str('%0.2f' % wgw96meanmags[i]), str('%0.2f' % wgw96evstd[i]), \
                       str('%0.2f' % a16meanmags[i]), str('%0.2f' % a16evstd[i]), \
                       str('%0.2f' % minDist[i]), str('%0.2f' % maxDist[i]),
                       str('%0.2f' % nrecs[i]))) + '\n'

# export lines
f = open('sea_ev_mags.csv', 'w')
f.write(lines)
f.close()

###############################################################################
# now plot
###############################################################################

def plt_stats(plt, rhyp, res):
    from numpy import nanmedian, nanstd, nanvar
    meanres = nanmedian(res)
    stdres  = nanstd(res)
    strng = 'Med = ' + str('%0.2f' % meanres) + '\n' \
          + 'Std  = ' + str('%0.2f' % stdres)
    plt.text(784, -0.735, strng, fontsize=12, va='bottom',ha='right')

mmin = [2.0, 4.0]
mmax = [6.0, 6.0]
bins = arange(15, 801, 30)
i = 0

letters = ['(a)', '(b)', '(c)', '(d)']

plt.figure(i, figsize=(18, 7))
for mi, ma in zip(mmin, mmax):
    plt.tick_params(labelsize=15)
    midx = where((seamags > mi) & (seamags <= ma))[0]
    
    # get gippsland data
    gidx = where((seamags > mi) & (seamags <= ma) & (lares < -37.8))[0]
    #plt.suptitle('ML '+str(mi)+'-'+str(ma), fontsize=22)
    
    # get MLM92
    plt.subplot(2,2,2*i+2)
    plt.plot(rhypres[midx], mlm92res[midx], '+', color='0.6')
    #plt.plot(rhypres[gidx], mlm92res[gidx], '+', color='seagreen')
    medbin, stdbin, meanx, outbins, nperbin = get_binned_stats(bins, rhypres[midx], mlm92res[midx])
    plt.errorbar(outbins, medbin, yerr=stdbin, color='r', fmt='s', mec='r')
    
    #medbin, stdbin, meanx, outbins, nperbin = get_binned_stats(bins, rhypres[gidx], mlm92res[gidx])
    #plt.errorbar(outbins, medbin, yerr=stdbin, color='b', fmt='s', mfc='none', mec='b')


    plt.plot([0, 800], [0, 0], 'k--')
    plt.xlim([0, 800])
    plt.ylim([-0.75, 0.75])
    if i == 0:
        plt.title('Michael-Lieba & Malafant (1992)', fontsize=16)
        # r'$\mathregular{R_{JB}}$
    plt.ylabel(r'MLM92 Stn $\mathregular{M_L}$ - Mean $\mathregular{M_L}$', fontsize=13)
    if i == 1:
        plt.xlabel('Hypocentral Distance (km)', fontsize=13)
    plt_stats(plt, rhypres[midx], mlm92res[midx])
    plt.text(15, 0.72, letters[2*i+1], fontsize=18, va='top', ha='left')
    '''
    # get WGW96
    plt.subplot(324)
    plt.plot(rhypres[midx], wgw96res[midx], '+', color='0.6')
    medbin, stdbin, meanx, outbins, nperbin = get_binned_stats(bins, rhypres[midx], wgw96res[midx])
    plt.errorbar(outbins, medbin, yerr=stdbin, color='r', fmt='s', mec='r')

    plt.plot([0, 800], [0, 0], 'k--')
    plt.xlim([0, 800])
    plt.ylim([-0.75, 0.75])
    plt.title('Wilkie et al. (1996)', fontsize=16)
    plt.ylabel('Stn ML - Mean ML', fontsize=14)
    plt.xlabel('Hypocentral Distance (km)', fontsize=14)
    plt_stats(plt, rhypres[midx], wgw96res[midx])
    '''
    '''
    # get A16
    plt.subplot(326)
    plt.plot(rhypres[midx], a16res[midx], '+', color='0.6')
    #plt.plot(rhypres[gidx], mlm92res[gidx], '+', color='seagreen')
    medbin, stdbin, meanx, outbins, nperbin = get_binned_stats(bins, rhypres[midx], a16res[midx])
    plt.errorbar(outbins, medbin, yerr=stdbin, color='r', fmt='s', mec='r')
    
    medbin, stdbin, meanx, outbins, nperbin = get_binned_stats(bins, rhypres[gidx], a16res[gidx])
    #plt.errorbar(outbins, medbin, yerr=stdbin, color='b', fmt='s', mfc='none', mec='b')

    plt.plot([0, 800], [0, 0], 'k--')
    plt.xlim([0, 800])
    plt.ylim([-0.75, 0.75])
    plt.title('Allen (2017)')
    plt.ylabel('Stn ML - Mean ML') 
    plt_stats(plt, rhypres[midx], a16res[midx])
    '''
    '''
    # get GG91
    plt.subplot(324)
    plt.plot(rhypres[midx], gg91res[midx], '+', color='0.6')
    medbin, stdbin, meanx, outbins, nperbin = get_binned_stats(bins, rhypres[midx], gg91res[midx])
    plt.errorbar(outbins, medbin, yerr=stdbin, color='r', fmt='s', mec='r')

    plt.plot([0, 800], [0, 0], 'k--')
    plt.xlim([0, 800])
    plt.ylim([-0.75, 0.75])
    plt.title('Gaull & Gregson (1991)')
    plt.ylabel('Stn ML - Mean ML')
    plt_stats(plt, rhypres[midx], gg91res[midx])
    '''
    # get R35
    plt.subplot(2,2,2*i+1)
    plt.plot(rhypres[midx], r35res[midx], '+', color='0.6')
    medbin, stdbin, meanx, outbins, nperbin = get_binned_stats(bins, rhypres[midx], r35res[midx])
    plt.errorbar(outbins, medbin, yerr=stdbin, color='r', fmt='s', mec='r')

    plt.plot([0, 800], [0, 0], 'k--')
    plt.xlim([0, 800])
    plt.ylim([-0.75, 0.75])
    if i == 0:
        plt.title('Richter (1935)', fontsize=16)
    plt.ylabel(r'R35 Stn $\mathregular{M_L}$ - Mean $\mathregular{M_L}$', fontsize=13)
    if i == 1:
        plt.xlabel('Hypocentral Distance (km)', fontsize=13)
    plt_stats(plt, rhypres[midx], r35res[midx])
    plt.text(15, 0.72, letters[2*i], fontsize=18, va='top', ha='left')
    
    # get HB87
    '''
    plt.subplot(326)
    plt.plot(rhypres[midx], hb87res[midx], '+', color='0.6')
    medbin, stdbin, meanx, outbins, nperbin = get_binned_stats(bins, rhypres[midx], hb87res[midx])
    plt.errorbar(outbins, medbin, yerr=stdbin, color='r', fmt='s', mec='r')

    plt.plot([0, 800], [0, 0], 'k--')
    plt.xlim([0, 800])
    plt.ylim([-0.75, 0.75])
    plt.title('Hutton & Boore (1987)')
    plt.ylabel('Stn ML - Mean ML')
    plt_stats(plt, rhypres[midx], hb87res[midx])
    '''
    '''
    # get BJ84
    plt.subplot(322)
    plt.plot(rhypres[midx], bj84res[midx], '+', color='0.6')
    medbin, stdbin, meanx, outbins, nperbin = get_binned_stats(bins, rhypres[midx], bj84res[midx])
    plt.errorbar(outbins, medbin, yerr=stdbin, color='r', fmt='s', mec='r')

    plt.plot([0, 800], [0, 0], 'k--')
    plt.xlim([0, 800])
    plt.ylim([-0.75, 0.75])
    plt.title('Bakun & Joyner (1984)', fontsize=16)
    plt.ylabel('Stn ML - Mean ML', fontsize=14)
    plt_stats(plt, rhypres[midx], bj84res[midx])
    '''
    '''
    # get Kradolfer
    plt.subplot(325)
    plt.plot(rhypres[midx], sed84res[midx], '+', color='0.6')
    medbin, stdbin, meanx, outbins, nperbin = get_binned_stats(bins, rhypres[midx], sed84res[midx])
    plt.errorbar(outbins, medbin, yerr=stdbin, color='r', fmt='s', mec='r')

    plt.plot([0, 800], [0, 0], 'k--')
    plt.xlim([0, 800])
    plt.ylim([-0.75, 0.75])
    plt.title('Kradolfer (1984)', fontsize=16)
    plt.ylabel('Stn ML - Mean ML', fontsize=14)
    plt.xlabel('Hypocentral Distance (km)', fontsize=14)
    plt_stats(plt, rhypres[midx], sed84res[midx])
    '''
    i += 1
    # save fig
pngfile = 'intra_eqn_res.png'
epsfile = 'intra_eqn_res.eps'
plt.savefig(pngfile, fmt='png', bbox_inches='tight')
plt.savefig(pngfile, fmt='eps', bbox_inches='tight')
"""    
###############################################################################
# now plot MLM92 vs A16
###############################################################################

plt.figure(5, figsize=(19, 13))
plt.subplot(231)
plt.plot(mlm92meanmags, bj84meanmags, 'r+', ms=12, mew=2.0)
plt.plot([1.0, 6.0], [1.0, 6.0], 'k--')
plt.xlabel('ML(MLM92)', fontsize=16)
plt.ylabel('ML(BJ84)', fontsize=16)
plt.xlim([1.0, 6.0])
plt.ylim([1.0, 6.0])

plt.subplot(437)
magdiff = mlm92meanmags - bj84meanmags
plt.plot(mlm92meanmags, magdiff, 'r+', ms=10, mew=1.0)
plt.plot([1.0, 6.0], [0.0, 0.0], 'k--')
plt.xlabel('ML(MLM92)', fontsize=16)
plt.ylabel('ML(MLM92) - ML(BJ84)', fontsize=16)
plt.xlim([1.0, 6.0])
plt.ylim([-0.3, 0.3])

plt.subplot(232)
plt.plot(mlm92meanmags, a16meanmags, 'r+', ms=12, mew=2.0)
plt.plot([1.0, 6.0], [1.0, 6.0], 'k--')
plt.xlabel('ML(MLM92)', fontsize=16)
plt.ylabel('ML(A16)', fontsize=16)
plt.xlim([1.0, 6.0])
plt.ylim([1.0, 6.0])

plt.subplot(438)
magdiff = mlm92meanmags - a16meanmags
plt.plot(mlm92meanmags, magdiff, 'r+', ms=10, mew=1.0)
plt.plot([1.0, 6.0], [0.0, 0.0], 'k--')
plt.xlabel('ML(MLM92)', fontsize=16)
plt.ylabel('ML(MLM92) - ML(A16)', fontsize=16)
plt.xlim([1.0, 6.0])
plt.ylim([-0.3, 0.3])

plt.subplot(233)
plt.plot(sed84meanmags, a16meanmags, 'r+', ms=12, mew=2.0)
plt.plot([1.0, 6.0], [1.0, 6.0], 'k--')
plt.xlabel('ML(SED)', fontsize=16)
plt.ylabel('ML(A16)', fontsize=16)
plt.xlim([1.0, 6.0])
plt.ylim([1.0, 6.0])

plt.subplot(439)
magdiff = sed84meanmags - a16meanmags
plt.plot(sed84meanmags, magdiff, 'r+', ms=10, mew=1.0)
plt.plot([1.0, 6.0], [0.0, 0.0], 'k--')
plt.xlabel('ML(SED)', fontsize=16)
plt.ylabel('ML(SED) - ML(A16)', fontsize=16)
plt.xlim([1.0, 6.0])
plt.ylim([-0.2, 0.4])

plt.savefig('mlm92_a16_cmp.png', fmt='png', bbox_inches='tight')
    
#plt.show()
    
###############################################################################
# now event stds
###############################################################################

fig = plt.figure(8)
plt.plot(r35evstd, 'co')   
plt.plot(hb87evstd, '^', color='limegreen')  
plt.plot(gg91evstd, 'bv')  
plt.plot(mlm92evstd, 's', color='orange')
plt.plot(wgw96evstd, 'rd') 
plt.plot(a16evstd, 'kh')

plt.show()
"""