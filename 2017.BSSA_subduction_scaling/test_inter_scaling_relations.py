'''
Code to read data and scaling curves and calculate residuals
'''

def get_binned_stats(bins, binsize, xdat, yres):
    from numpy import array, nanstd, where
    from scipy.stats import nanmean

    medres = []
    stdres = []

    yres = array(yres)

    halfbin = binsize / 2.0

    for bin in bins:
        index = array(where((xdat >= bin-halfbin) & (xdat < bin+halfbin))[0])
        medres.append(nanmean(yres[index]))
        stdres.append(nanstd(yres[index]))

    return array(medres), array(stdres)

from numpy import arange, array, log10, mean, median, std, nanmax, unique, where, interp, \
                                 vstack, percentile, nanmean, nanmedian, nanstd, nan, isnan, exp
import matplotlib.pyplot as plt
from misc_tools import plttext, weighted_avg_and_std
from scipy.stats import ttest_ind
import pandas
from pandas import DataFrame, Series

# read data file
print 'read data file...'
lines = open('summary_trimmed_table.txt', 'rb').readlines()[1:]

evdict = []
for line in lines:
    ev = {}
    dat = line.strip().split('\t')
    ev['id'] = dat[0]
    ev['region'] = dat[1]
    #ev['evloc'] = dat[2]
    ev['evlon'] = float(dat[2])
    ev['evlat'] = float(dat[3])
    ev['evdep'] = float(dat[4])
    ev['evmag'] = float(dat[5])
    ev['evtype'] = dat[6]
    ev['flen'] = float(dat[7])
    ev['fwid'] = float(dat[8])
    ev['farea'] = float(dat[9])
    ev['fmxs'] = float(dat[10])
    ev['favs'] = float(dat[11])
    ev['dip'] = float(dat[12])

    evdict.append(ev)

# read linera len model coeffs
lines = open('len1_coeffs.txt', 'rb').readlines()
lclin1 = float(lines[0].strip())
lclin2 = float(lines[1].strip())

# read len model coeffs
lines = open('len_coeffs.txt', 'rb').readlines()
lc1 = float(lines[0].strip())
lc2 = float(lines[1].strip())
lc3 = float(lines[2].strip())
lc4 = float(lines[3].strip())

# read linera wid model coeffs
lines = open('wid1_coeffs.txt', 'rb').readlines()
wclin1 = float(lines[0].strip())
wclin2 = float(lines[1].strip())

# read bi-linear width coefficients
lines = open('bilin_wid_coeffs.txt', 'rb').readlines()
wc1 = float(lines[0].strip())
wc2 = float(lines[1].strip())
wc3 = float(lines[2].strip())

# read linera wid model coeffs
lines = open('area1_coeffs.txt', 'rb').readlines()
aclin1 = float(lines[0].strip())
aclin2 = float(lines[1].strip())

# read bi-linear width coefficients
lines = open('area_coeffs.txt', 'rb').readlines()
ac1 = float(lines[0].strip())
ac2 = float(lines[1].strip())
ac3 = float(lines[2].strip())
ac4 = float(lines[3].strip())

# read lw coeffs
lines = open('LW_coeffs.txt', 'rb').readlines()
lw1 = float(lines[0].strip())
lw2 = float(lines[1].strip())
lw3 = float(lines[2].strip())

# read max slip model coeffs
lines = open('mxs_coeffs.txt', 'rb').readlines()
ms1 = float(lines[0].strip())
ms2 = float(lines[1].strip())

# read max slip model coeffs
lines = open('avs_coeffs.txt', 'rb').readlines()
as1 = float(lines[0].strip())
as2 = float(lines[1].strip())

'''
get len and wid residuals
'''
# interface
rlen = []
rwid = []
rarea = []
rlenlin = []
rwidlin = []
rarealin = []
mwrlen = []
mwrwid = []
mwrarea = []
mwrlenlin = []
mwrwidlin = []
mwrarealin = []
mag = []
dip = []
reg = []
rmxs = []
ravs = []
mwrmxs = []
mwravs = []
dep = []
rlwwid = []
rlwlen = []
evid = []

# insalb
intra_rlen = []
intra_rwid = []
intra_rarea = []
intra_mag = []
intra_dip = []
intra_reg = []
intra_ztor = []
intra_dep = []
intra_rlw = []
intra_evid = []

# outer rise
outer_rlen = []
outer_rwid = []
outer_rarea = []
outer_mag = []
outer_dip = []
outer_reg = []
outer_ztor = []
outer_dep = []
outer_rlw = []
outer_evid = []

# transform
trans_rlen = []
trans_rwid = []
trans_rarea = []
trans_mag = []
trans_dip = []
trans_reg = []
trans_ztor = []
trans_dep = []
trans_rlw = []
trans_evid = []

stattxt = ''

for ev in evdict:
    #########################
    # linear len
    ev['clenlin'] = 10**(lclin1 + lclin2 * ev['evmag'])
    ev['clenlinmw'] = (log10(ev['flen']) - lclin1) / lclin2
    
    # bi-linear len
    ev['clen'] = 10**(lc2 + lc1 * ev['evmag'])
    ylen = lc2 + lc1 * lc4
    if ev['evmag'] > lc4:
        ev['clen'] = 10**(lc3 * (ev['evmag']-lc4) + ylen)
        
    if ev['evmag'] > lc4:
        ev['clen'] = 10**(lc3 * (ev['evmag'] - lc4) + ylen)
    
    if ev['flen'] < 10**(lc2 + lc1 * lc4) and ev['evmag'] < lc4:
        ev['clenmw'] = (log10(ev['flen']) - lc2) / lc1
    
    elif ev['flen'] < 10**(lc3 * (9.5-lc4) + ylen):
        ev['clenmw'] = ((log10(ev['flen']) - ylen) / lc3) + lc4
    
    else:
        ev['clenmw'] = nan
    #########################
    # linear wid
    ev['cwidlin'] = 10**(wclin1 + wclin2 * ev['evmag'])
    ev['cwidlinmw'] = (log10(ev['fwid']) - wclin1) / wclin2
    
    # bi-linear width
    ev['cwid'] = 10**min([wc2 + wc1*ev['evmag'], wc2 + wc1*wc3])
    if ev['fwid'] < 10**(wc2 + wc1*wc3) and ev['evmag'] < wc3:
        ev['cwidmw'] = (log10(ev['fwid']) - wc2) / wc1
    else:
        ev['cwidmw'] = nan
    
    #########################
    # linear area
    ev['carealin'] = 10**(aclin1 + aclin2 * ev['evmag'])
    ev['carealinmw'] = (log10(ev['farea']) - aclin1) / aclin2
    
    # bi-linear area
    ev['carea'] = 10**(ac2 + ac1 * ev['evmag'])
    yarea = ac2 + ac1 * ac4
    if ev['evmag'] > ac4:
        ev['carea'] = 10**(ac3 * (ev['evmag'] - ac4) + yarea)
    
    # inverse
    if ev['farea'] < 10**(ac2 + ac1 * ac4) and ev['evmag'] < ac4:
        ev['careamw'] = (log10(ev['farea']) - ac2) / ac1
    
    elif ev['farea'] < 10**(ac3 * (9.5-ac4) + yarea):
        ev['careamw'] = ((log10(ev['farea']) - yarea) / ac3) + ac4
    
    else:
        ev['careamw'] = nan
        
    #########################
    # calc len-width
    ev['clwwid'] = 10**(lw2 + lw2 * log10(ev['flen']))
    if ev['clwwid'] > 10**lw3:
        ev['clwwid'] = 10**lw3
        
    # bi-linear L-W
    if ev['fwid'] < 10**lw3 and ev['flen'] < 10**((log10(ev['fwid']) - lw2) / lw1):
        ev['clwlen'] = 10**((log10(ev['fwid']) - lw2) / lw1)
    else:
        ev['clwlen'] = nan    
    #########################   
    # calc max slip
    ev['cmxs'] = 10**(ms1 + ms2 * ev['evmag'])
    ev['cmxsmw'] = (log10(ev['fmxs']) - ms1) / ms2
    #########################
    # calc max slip
    ev['cavs'] = 10**(as1 + as2 * ev['evmag'])
    ev['cavsmw'] = (log10(ev['favs']) - as1) / as2
    
    #########################
    # get log residuals
    ev['rlen'] = log10(ev['flen'] / ev['clen'])
    ev['rlenlin'] = log10(ev['flen'] / ev['clenlin'])
    ev['rwid'] = log10(ev['fwid'] / ev['cwid'])
    ev['rwidlin'] = log10(ev['fwid'] / ev['cwidlin'])
    ev['rarea'] = log10(ev['farea'] / ev['carea'])
    ev['rarealin'] = log10(ev['farea'] / ev['carealin'])
    ev['rmxs'] = log10(ev['fmxs'] / ev['cmxs'])
    ev['ravs'] = log10(ev['favs'] / ev['cavs'])
    ev['rlwwid'] = log10(ev['fwid'] / ev['clwwid'])
    ev['rlenlinmw'] = ev['evmag'] - ev['clenlinmw']
    ev['rarealinmw'] = ev['evmag'] - ev['carealinmw']
    ev['rwidlinmw'] = ev['evmag'] - ev['cwidlinmw']
    ev['rlenmw'] = ev['evmag'] - ev['clenmw']
    ev['rareamw'] = ev['evmag'] - ev['careamw']
    ev['rwidmw'] = ev['evmag'] - ev['cwidmw']
    ev['rmxsmw'] = ev['evmag'] - ev['cmxsmw']
    ev['ravsmw'] = ev['evmag'] - ev['cavsmw']
    ev['rlwlen'] = log10(ev['flen'] / ev['clwlen'])
    
    # get only interface events
    if ev['evtype'] == 'i' or ev['evtype'] == 'h':
        rlen.append(ev['rlenlin']) # not using L2 anymore
        rlenlin.append(ev['rlenlin'])
        rwid.append(ev['rwid'])
        rwidlin.append(ev['rwidlin'])
        rarea.append(ev['rarea'])
        rarealin.append(ev['rarealin'])
        mag.append(ev['evmag'])
        dip.append(ev['dip'])
        reg.append(ev['region'])
        rmxs.append(ev['rmxs'])
        ravs.append(ev['ravs'])
        dep.append(ev['evdep'])
        #rmxs.append(ev['rmxs'])
        #ravs.append(ev['ravs'])
        rlwwid.append(ev['rlwwid'])
        evid.append(ev['id'])
        
        mwrlen.append(ev['rlenmw'])
        mwrlenlin.append(ev['rlenlinmw'])
        mwrwid.append(ev['rwidmw'])
        mwrwidlin.append(ev['rwidlinmw'])
        mwrarea.append(ev['rareamw'])
        mwrarealin.append(ev['rarealinmw'])
        mwrmxs.append(ev['rmxsmw'])
        mwravs.append(ev['ravsmw'])
        rlwlen.append(ev['rlwlen'])

    elif ev['evtype'] == 's':
        intra_rlen.append(ev['rlenlin'])
        intra_rwid.append(ev['rwid'])
        intra_rarea.append(ev['rarea'])
        intra_mag.append(ev['evmag'])
        intra_dip.append(ev['dip'])
        intra_reg.append(ev['region'])
#        intra_ztor.append(ev['fztor'])
        intra_dep.append(ev['evdep'])
        intra_rlw.append(ev['rlwwid'])
        intra_evid.append(ev['id'])

    elif ev['evtype'] == 'o':
        outer_rlen.append(ev['rlenlin'])
        outer_rwid.append(ev['rwid'])
        outer_rarea.append(ev['rarea'])
        outer_mag.append(ev['evmag'])
        outer_dip.append(ev['dip'])
        outer_reg.append(ev['region'])
#        outer_ztor.append(ev['fztor'])
        outer_dep.append(ev['evdep'])
        outer_rlw.append(ev['rlwwid'])
        outer_evid.append(ev['id'])
    elif ev['evtype'] == 't':
        trans_rlen.append(ev['rlenlin'])
        trans_rwid.append(ev['rwid'])
        trans_rarea.append(ev['rarea'])
        trans_mag.append(ev['evmag'])
        trans_dip.append(ev['dip'])
        trans_reg.append(ev['region'])
#        trans_ztor.append(ev['fztor'])
        trans_dep.append(ev['evdep'])
        trans_rlw.append(ev['rlwwid'])
        trans_evid.append(ev['id'])

# get interface wts for weigted std
im = [7.99, 8.0]
isi = [0.0, 0.2]
sx = interp(array(mag), im, isi)

wts = sx / sum(sx)


#weighted_avg_and_std(values, weights)

# print residual sigma
sigmatxt  = 'len1,' + str(nanstd(rlenlin))+',' + str(weighted_avg_and_std(rlenlin, wts)[1])+','+ str(nanstd(mwrlenlin))+ '\n'
sigmatxt += 'len2,' + str(nanstd(rlen))+','+ str(weighted_avg_and_std(rlenlin, wts)[1])+','+ str(nanstd(mwrlen))+ '\n'
sigmatxt += 'wid1,' + str(nanstd(rwidlin))+','+ str(weighted_avg_and_std(rwidlin, wts)[1])+','+ str(nanstd(mwrwidlin))+ '\n'
sigmatxt += 'wid2,' + str(nanstd(rwid))+','+ str(weighted_avg_and_std(rwid, wts)[1])+','+ str(nanstd(mwrwid))+ '\n'
sigmatxt += 'area1,' + str(nanstd(rarealin))+','+ str(weighted_avg_and_std(rarealin, wts)[1])+','+ str(nanstd(mwrarealin))+ '\n'
sigmatxt += 'area2,' + str(nanstd(rarea))+','+ str(weighted_avg_and_std(rarea, wts)[1])+','+ str(nanstd(mwrarea))+ '\n'
rmxs = array(rmxs)
idx = isnan(rmxs)==False  
idxwts = sx[idx] / sum(sx[idx])
sigmatxt += 'mxs,' + str(nanstd(rmxs))+','+ str(weighted_avg_and_std(rmxs[idx], idxwts)[1])+','+ str(nanstd(mwrmxs))+ '\n'
ravs = array(ravs)
idx = isnan(ravs)==False  
idxwts = sx[idx] / sum(sx[idx])
sigmatxt += 'avs,' + str(nanstd(ravs))+','+ str(weighted_avg_and_std(ravs[idx], idxwts)[1])+','+ str(nanstd(mwravs))+ '\n'
sigmatxt += 'lww-l,' + str(nanstd(rlwwid))+','+ str(weighted_avg_and_std(rlwwid, wts)[1])+','+ str(nanstd(rlwlen))+ '\n'

f = open('interface_sigmas.csv', 'wb')
f.write(sigmatxt)
f.close()

'''
plot len residuals
'''
# plot res vs MW
fig = plt.figure(1, figsize=(10,9))
ax = plt.subplot(2, 2, 1)
plt.plot([6.5, 10], [0, 0], 'k--')
plt.plot(mag, rlen, '+', color='0.5', markersize=8)
magbin = arange(6.75, 9.8, 0.5)
medres, stdres = get_binned_stats(magbin, 0.5, mag, rlen)
plt.errorbar(magbin, medres, stdres, fmt='rs', mec='r', mfc='r')

plt.xlabel('Magnitude (MW)')
plt.ylabel('log10(Obs / Pred)')
plt.suptitle('Fault Length Residuals')
plt.ylim([-0.75, 0.75])
#ylims = log10(ax.get_ylim())
#ymax = 10**(0.92*(ylims[1]-ylims[0]) + ylims[0])


# plot res vs dip
ax = plt.subplot(2, 2, 2)
plt.plot(dip, rlen, '+', color='0.5', markersize=8)
dipbin = arange(2.5, 45., 5.)
medres, stdres = get_binned_stats(dipbin, 5., dip, rlen)
plt.errorbar(dipbin, medres, stdres, fmt='rs', mec='r', mfc='r')

plt.plot([0, 45], [0, 0], 'k--')
plt.xlabel('Centroid Dip (degrees)')
plt.ylabel('log10(Obs / Pred)')
plt.ylim([-0.75, 0.75])

# plot box/whisker plots per subduction
regset = unique(reg)[1:]
regset = list(unique(reg)[1:])
regset.append('')

data = []
tmplab = []
stattxt = 'Rupture Length by Region\n'
rlen = array(rlen)
rlenlin = array(rlenlin)
for sz in regset:
    # get region index
    index1 = array(where(array(reg) == sz))[0]
    index2 = array(where(array(reg) != sz))[0]
    if sz == '':
        sz = 'OTH'
        
    if len(index1) > 1:
        data.append(rlen[index1])
        tmplab.append(sz+' '+str(len(index1)))

        # do t-test
        '''
        reject null hypothesis (i.e. that medians are the same) if p < 0.05
        '''
        print sz,'len t-test:',ttest_ind(rlen[index2], rlen[index1], equal_var=False), median(rlen[index1]) # Welch's t test
        stattxt += ','.join((sz, str('%0.2f' % median(rlen[index1])), str('%0.3f' % std(rlen[index1])), 
                            str('%0.3f' % ttest_ind(rlen[index2], rlen[index1], equal_var=False)[1]))) + '\n'
                            
stattxt += ','.join(('ALL', str('%0.2f' % median(rlen)), str('%0.3f' % std(rlen)), 
                     str('%0.3f' % ttest_ind(rlen, rlen, equal_var=False)[1]))) + '\n'
print '\n'
stattxt += '\n'

ax = plt.subplot(2, 2, 3)
bp = plt.boxplot(data)
plt.setp(bp['boxes'], color='black')
plt.setp(bp['whiskers'], color='black')
plt.setp(bp['medians'], color='black', lw=2)
plt.plot([0, len(tmplab)+1], [0, 0], 'k--')
plt.xticks(range(1,len(tmplab)+1),tmplab, rotation=45)

plt.xlabel('Subduction Zone')
plt.ylabel('log10(Obs / Pred)')
plt.ylim([-0.75, 0.75])

'''
# plot res vs ztor
ax = plt.subplot(2, 2, 4)
plt.plot([-2, 80], [0, 0], 'k--')
plt.plot(ztor, rlen, '+', color='0.5', markersize=8)
zbin = arange(2.5, 40, 5.)
medres, stdres = get_binned_stats(zbin, 5., ztor, rlen)
plt.errorbar(zbin, medres, stdres, fmt='rs', mec='r', mfc='r')
plt.xlabel('Depth to top of Rupture (km)')
plt.ylabel('log10(Obs / Pred)')
plt.ylim([-0.75, 0.75])
plt.xlim([-2, 35])
'''
plt.savefig('length_res.pdf',dpi=300,format='pdf')

'''
plot wid residuals
'''
# plot res vs MW
fig = plt.figure(2, figsize=(10,9))
ax = plt.subplot(2, 2, 1)
plt.plot([6.5, 10], [0, 0], 'k--')
plt.plot(mag, rwid, '+', color='0.5', markersize=8)
medres, stdres = get_binned_stats(magbin, 0.5, mag, rwid)
plt.errorbar(magbin, medres, stdres, fmt='rs', mec='r', mfc='r')
plt.xlabel('Magnitude (MW)')
plt.ylabel('log10(Obs / Pred)')
plt.suptitle('Fault Width Residuals')
plt.ylim([-0.75, 0.75])

# plot res vs dip
ax = plt.subplot(2, 2, 2)
plt.plot([0, 45], [0, 0], 'k--')
plt.plot(dip, rwid, '+', color='0.5', markersize=8)
medres, stdres = get_binned_stats(dipbin, 5., dip, rwid)
plt.errorbar(dipbin, medres, stdres, fmt='rs', mec='r', mfc='r')
plt.xlabel('Centroid Dip (degrees)')
plt.ylabel('log10(Obs / Pred)')
plt.ylim([-0.75, 0.75])

# plot box/whisker plots of residuals per subduction
data = []
tmplab = []
stattxt += 'Rupture Width by Region\n'
rwid = array(rwid)
for sz in regset:
    # get region index
    index1 = array(where(array(reg) == sz))[0]
    index2 = array(where(array(reg) != sz))[0]
    if sz == '':
        sz = 'OTH'
    
    if len(index1) > 1:
        data.append(rwid[index1])
        tmplab.append(sz+' '+str(len(index1)))

        # do t-test
        print sz,'wid t-test:',ttest_ind(rwid[index2], rwid[index1], equal_var=False), median(rwid[index1]) # Welch's t test
        stattxt += ','.join((sz, str('%0.2f' % median(rwid[index1])), str('%0.3f' % std(rwid[index1])), 
                            str('%0.3f' % ttest_ind(rwid[index2], rwid[index1], equal_var=False)[1]))) + '\n'

stattxt += ','.join(('ALL', str('%0.2f' % median(rwid)), str('%0.3f' % std(rwid)), 
                     str('%0.3f' % ttest_ind(rwid, rwid, equal_var=False)[1]))) + '\n'

print '\n'
stattxt += '\n'

ax = plt.subplot(2, 2, 3)
bp = plt.boxplot(data)
plt.plot([0, len(tmplab)+1], [0, 0], 'k--')
plt.xticks(range(1,len(tmplab)+1),tmplab, ha='right', rotation=45)

plt.xlabel('Subduction Zone')
plt.ylabel('log10(Obs / Pred)')
plt.ylim([-0.75, 0.75])

# plot box/whisker plots of dips per subduction
data = []
tmplab = []
noidx = []
dip = array(dip)
for sz in regset:
    # get region index
    index = array(where((array(reg) == sz) & (dip > 0)))[0]
    if len(index) > 1:
        data.append(dip[index])
        tmplab.append(sz+' '+str(len(index)))

ax = plt.subplot(2, 2, 4)
bp = plt.boxplot(data)
plt.plot([0, len(tmplab)+1], [median(dip), median(dip)], 'k--')
plt.xticks(range(1,len(tmplab)+1),tmplab, ha='right', rotation=45)

plt.xlabel('Subduction Zone')
plt.ylabel('Centroid Dip')
plt.ylim([0, 50])
plt.savefig('width_res.pdf',dpi=200,format='pdf')


'''
plot area residuals
'''
# plot res vs MW
fig = plt.figure(3, figsize=(10,9))
ax = plt.subplot(2, 2, 1)
plt.plot([6.5, 10], [0, 0], 'k--')
plt.plot(mag, rarea, '+', c='0.5', markersize=8)
medres, stdres = get_binned_stats(magbin, 0.5, mag, rarea)
plt.errorbar(magbin, medres, stdres, fmt='rs', mec='r', mfc='r')
plt.xlabel('Magnitude (MW)')
plt.ylabel('log10(Obs / Pred)')
plt.suptitle('Rupture Area Residuals')
plt.ylim([-0.75, 0.75])

# plot res vs dip
ax = plt.subplot(2, 2, 2)
plt.plot([0, 45], [0, 0], 'k--')
plt.plot(dip, rarea, '+', color='0.5', markersize=8)
medres, stdres = get_binned_stats(dipbin, 5., dip, rarea)
plt.errorbar(dipbin, medres, stdres, fmt='rs', mec='r', mfc='r')
plt.xlabel('Centroid Dip (degrees)')
plt.ylabel('log10(Obs / Pred)')
plt.ylim([-0.75, 0.75])

# plot box/whisker plots of residuals per subduction
data = []
tmplab = []
stattxt += 'Rupture Area by Region\n'
rarea = array(rarea)
for sz in regset:
    # get region index
    index1 = array(where(array(reg) == sz))[0]
    index2 = array(where(array(reg) != sz))[0]
    if sz == '':
        sz = 'OTH'
        
    if len(index1) > 1:
        data.append(rarea[index1])
        tmplab.append(sz+' '+str(len(index1)))
        # do t-test
        print sz,'area t-test:',ttest_ind(rarea[index2], rarea[index1], equal_var=False), median(rarea[index1]) # Welch's t test
        stattxt += ','.join((sz, str('%0.2f' % median(rarea[index1])), str('%0.3f' % std(rarea[index1])), 
                            str('%0.3f' % ttest_ind(rarea[index2], rarea[index1], equal_var=False)[1]))) + '\n'

stattxt += ','.join(('ALL', str('%0.2f' % median(rarea)), str('%0.3f' % std(rarea)), 
                     str('%0.3f' % ttest_ind(rarea, rarea, equal_var=False)[1]))) + '\n'

print '\n'
stattxt += '\n'

ax = plt.subplot(2, 2, 3)
bp = plt.boxplot(data)
plt.plot([0, len(tmplab)+1], [0, 0], 'k--')
plt.xticks(range(1,len(tmplab)+1),tmplab, ha='right', rotation=45)

plt.xlabel('Subduction Zone')
plt.ylabel('log10(Obs / Pred)')
plt.ylim([-0.75, 0.75])

'''
# plot hist of zbor
ax = plt.subplot(2, 2, 4)
plt.hist(zbor,bins=arange(7.5, 101, 5))
plt.xlabel('Depth to Bottom of Rupture')
plt.ylabel('Count')
perctxt = '75th Percentile = ' + str('%0.1f' % percentile(zbor, 75.))
plttext(plt ,55, 9.25, perctxt, fsize=11)
'''
plt.savefig('area_res.png',dpi=300,format='png')

##############################################################################
'''
plot inter vs intra vs outer, trans
'''
# do bw plots for inter/intra/outerr
fig = plt.figure(4, figsize=(15,9))

szset = ('Interface '+str(len(rlen)), 'Intraslab '+str(len(intra_rlen)), \
         'Outer-Rise '+str(len(outer_rlen)),'Strike-Slip '+str(len(trans_rlen)))
rldata = array([rlenlin, intra_rlen, outer_rlen, trans_rlen])
rwdata = array([rwid, intra_rwid, outer_rwid, trans_rwid])
radata = array([rarea, intra_rarea, outer_rarea, trans_rarea])


print 'Interface','len t-test:',ttest_ind(rlen, rlen, equal_var=False), median(rlen) # Welch's t test
print 'In-Slab','len t-test:',ttest_ind(rlen, intra_rlen, equal_var=False), median(intra_rlen) # Welch's t test
print 'Outer-Rise','len t-test:',ttest_ind(rlen, outer_rlen, equal_var=False), median(outer_rlen) # Welch's t test
print 'Strike-Slip','len t-test:',ttest_ind(rlen, trans_rlen, equal_var=False), median(trans_rlen) # Welch's t test
stattxt += 'Rupture Length by Event Type\n'
stattxt += ','.join(('Interface', str('%0.2f' % median(rlen)), str('%0.3f' % std(rlen)), 
                    str('%0.3f' % ttest_ind(rlen, rlen, equal_var=False)[1]))) + '\n'
stattxt += ','.join(('Intraslab', str('%0.2f' % median(intra_rlen)), str('%0.3f' % std(intra_rlen)), 
                    str('%0.3f' % ttest_ind(rlen, intra_rlen, equal_var=False)[1]))) + '\n'
stattxt += ','.join(('Outer-Rise', str('%0.2f' % median(outer_rlen)), str('%0.3f' % std(outer_rlen)), 
                    str('%0.3f' % ttest_ind(rlen, outer_rlen, equal_var=False)[1]))) + '\n'
stattxt += ','.join(('Strike-Slip', str('%0.2f' % median(trans_rlen)), str('%0.3f' % std(trans_rlen)), 
                    str('%0.3f' % ttest_ind(rlen, trans_rlen, equal_var=False)[1]))) + '\n\n'
print '\n'
print 'Interface','wid t-test:',ttest_ind(rwid, rwid, equal_var=False), median(rwid) # Welch's t test
print 'In-Slab','wid t-test:',ttest_ind(rwid, intra_rwid, equal_var=False), median(intra_rwid) # Welch's t test
print 'Outer-Rise','wid t-test:',ttest_ind(rwid, outer_rwid, equal_var=False), median(outer_rwid) # Welch's t test
print 'Strike-Slip','wid t-test:',ttest_ind(rwid, trans_rwid, equal_var=False), median(trans_rwid) # Welch's t test
stattxt += 'Rupture Width by Event Type\n'
stattxt += ','.join(('Interface', str('%0.2f' % median(rwid)), str('%0.3f' % std(rwid)), 
                    str('%0.3f' % ttest_ind(rwid, rwid, equal_var=False)[1]))) + '\n'
stattxt += ','.join(('Intraslab', str('%0.2f' % median(intra_rwid)), str('%0.3f' % std(intra_rwid)), 
                    str('%0.3f' % ttest_ind(rwid, intra_rwid, equal_var=False)[1]))) + '\n'
stattxt += ','.join(('Outer-Rise', str('%0.2f' % median(outer_rwid)), str('%0.3f' % std(outer_rwid)), 
                    str('%0.3f' % ttest_ind(rwid, outer_rwid, equal_var=False)[1]))) + '\n'
stattxt += ','.join(('Strike-Slip', str('%0.2f' % median(trans_rwid)), str('%0.3f' % std(trans_rwid)), 
                    str('%0.3f' % ttest_ind(rwid, trans_rwid, equal_var=False)[1]))) + '\n\n'
print '\n'
print 'Interface','area t-test:',ttest_ind(rarea, rarea, equal_var=False), median(rarea) # Welch's t test
print 'In-Slab','area t-test:',ttest_ind(rarea, intra_rarea, equal_var=False), median(intra_rarea) # Welch's t test
print 'Outer-Rise','area t-test:',ttest_ind(rarea, outer_rarea, equal_var=False), median(outer_rarea) # Welch's t test
print 'Strike-Slip','area t-test:',ttest_ind(rarea, trans_rarea, equal_var=False), median(trans_rarea) # Welch's t test
stattxt += 'Rupture Area by Event Type\n'
stattxt += ','.join(('Interface', str('%0.2f' % median(rarea)), str('%0.3f' % std(rarea)), 
                    str('%0.3f' % ttest_ind(rarea, rarea, equal_var=False)[1]))) + '\n'
stattxt += ','.join(('Intraslab', str('%0.2f' % median(intra_rarea)), str('%0.3f' % std(intra_rarea)), 
                    str('%0.3f' % ttest_ind(rarea, intra_rarea, equal_var=False)[1]))) + '\n'
stattxt += ','.join(('Outer-Rise', str('%0.2f' % median(outer_rarea)), str('%0.3f' % std(outer_rarea)), 
                    str('%0.3f' % ttest_ind(rarea, outer_rarea, equal_var=False)[1]))) + '\n'
stattxt += ','.join(('Strike-Slip', str('%0.2f' % median(trans_rarea)), str('%0.3f' % std(trans_rarea)), 
                    str('%0.3f' % ttest_ind(rarea, trans_rarea, equal_var=False)[1]))) + '\n\n'
print '\n'

print 'Interface','area t-test:',ttest_ind(rlwwid, rlwwid, equal_var=False), median(rlwwid) # Welch's t test
print 'In-Slab','area t-test:',ttest_ind(rlwwid, intra_rlw, equal_var=False), median(intra_rlw) # Welch's t test
print 'Outer-Rise','area t-test:',ttest_ind(rlwwid, outer_rlw, equal_var=False), median(outer_rlw) # Welch's t test
print 'Strike-Slip','area t-test:',ttest_ind(rlwwid, trans_rlw, equal_var=False), median(trans_rlw) # Welch's t test
stattxt += 'Length-Width Residuals by Event Type\n'
stattxt += ','.join(('Interface', str('%0.2f' % median(rlwwid)), str('%0.3f' % std(rlwwid)), 
                    str('%0.3f' % ttest_ind(rlwwid, rlwwid, equal_var=False)[1]))) + '\n'
stattxt += ','.join(('Intraslab', str('%0.2f' % median(intra_rlw)), str('%0.3f' % std(intra_rlw)), 
                    str('%0.3f' % ttest_ind(rlwwid, intra_rlw, equal_var=False)[1]))) + '\n'
stattxt += ','.join(('Outer-Rise', str('%0.2f' % median(outer_rlw)), str('%0.3f' % std(outer_rlw)), 
                    str('%0.3f' % ttest_ind(rlwwid, outer_rlw, equal_var=False)[1]))) + '\n'
stattxt += ','.join(('Strike-Slip', str('%0.2f' % median(trans_rlw)), str('%0.3f' % std(trans_rlw)), 
                    str('%0.3f' % ttest_ind(rlwwid, trans_rlw, equal_var=False)[1]))) + '\n\n'
print '\n'

# len boxplots
ax = plt.subplot(2, 3, 1)
xbox = [0.53, 1.47, 1.47, 0.53, 0.53]
ybox = [.985, .985, -.985, -.985, .985]
plt.plot(xbox, ybox, '-', color='gray', lw=2.5)
bp = plt.boxplot(rldata)
plt.setp(bp['boxes'], color='black')
plt.setp(bp['whiskers'], color='black')
plt.setp(bp['medians'], color='black')
plt.plot([0, len(szset)+1], [0, 0], 'k--')
#plttext(ax, 0.65, 0.85, '(A)', fsize='large')
plt.xticks(range(1,len(szset)+1),szset, ha='right', rotation=45)

#plt.xlabel('Subduction Type')
plt.ylabel(r'log$\mathregular{_{10}}$(Obs / Pred)', fontsize=14)
plt.title('Rupture Length Residuals')
plt.ylim([-1, 1])
plt.text(4., 0.82, '(A)', fontsize=18)


# wid boxplots
ax = plt.subplot(2, 3, 2)
plt.plot(xbox, ybox, '-', color='gray', lw=2.5)
bp = plt.boxplot(rwdata)
plt.setp(bp['boxes'], color='black')
plt.setp(bp['whiskers'], color='black')
plt.setp(bp['medians'], color='black')
plt.plot([0, len(szset)+1], [0, 0], 'k--')
#plttext(ax, 0.65, 0.85, '(B)', fsize='large')
plt.xticks(range(1,len(szset)+1),szset, ha='right', rotation=45)

#plt.xlabel('Subduction Type')
#plt.ylabel('log10(Obs / Pred)')
plt.title('Rupture Width Residuals')
plt.ylim([-1, 1])
plt.text(4., 0.82, '(B)', fontsize=18)

# area boxplots
ax = plt.subplot(2, 3, 3)
plt.plot(xbox, ybox, '-', color='gray', lw=2.5)
bp = plt.boxplot(radata)
plt.setp(bp['boxes'], color='black')
plt.setp(bp['whiskers'], color='black')
plt.setp(bp['medians'], color='black')
plt.plot([0, len(szset)+1], [0, 0], 'k--')
#plttext(ax, 0.65, 0.85, '(C)', fsize='large')
plt.xticks(range(1,len(szset)+1),szset, ha='right', rotation=45)

#plt.xlabel('Subduction Type')
#plt.ylabel('log10(Obs / Pred)')
plt.title('Rupture Area Residuals')
plt.ylim([-1, 1])
plt.text(4., 0.82, '(C)', fontsize=18)

plt.savefig('inter_intra_outer_bwplot.png',dpi=300,format='png', bbox_inches='tight')
plt.savefig('inter_intra_outer_bwplot.pdf',dpi=300,format='pdf', bbox_inches='tight')
plt.show()

# export statistics
f = open('residual_stats.csv', 'wb')
f.write(stattxt)
f.close()

##############################################################################
'''
plot region
'''
fig = plt.figure(5, figsize=(15,9))

# do length
data = []
tmplab = []
rlen = array(rlen)
for sz in regset:
    
    # get region index
    index = array(where(array(reg) == sz))[0]
    if sz == '':
        sz = 'OTH'
    if len(index) > 1:
        data.append(rlen[index])
        tmplab.append(sz+' '+str(len(index)))
        
    else:
        print sz, index
        
ax = plt.subplot(2, 3, 1)
bp = plt.boxplot(data)
plt.setp(bp['boxes'], color='black')
plt.setp(bp['whiskers'], color='black')
plt.setp(bp['medians'], color='black')
plt.plot([0, len(tmplab)+1], [0, 0], 'k--')
plttext(ax, 0.65, 0.51, '(A)', fsize='large')
plt.xticks(range(1,len(tmplab)+1),tmplab, ha='right', rotation=45)

plt.title('Interface Rupture Length Residuals', fontsize=14)
plt.ylabel(r'log$\mathregular{_{10}}$(Obs / Pred)', fontsize=14)
plt.ylim([-0.75, 0.75])
plt.text(0.65, 0.62, '(A)', fontsize=18)


# do width
data = []
tmplab = []
rwid = array(rwid)
for sz in regset:
    # get region index
    index = array(where(array(reg) == sz))[0]
    if sz == '':
        sz = 'OTH'
    if len(index) > 1:
        data.append(rwid[index])
        tmplab.append(sz+' '+str(len(index)))
        
ax = plt.subplot(2, 3, 2)
bp = plt.boxplot(data)
plt.setp(bp['boxes'], color='black')
plt.setp(bp['whiskers'], color='black')
plt.setp(bp['medians'], color='black')
plt.plot([0, len(tmplab)+1], [0, 0], 'k--')
plttext(ax, 0.65, 0.51, '(B)', fsize='large')
plt.xticks(range(1,len(tmplab)+1),tmplab, ha='right', rotation=45)

plt.title('Interface Rupture Width Residuals', fontsize=14)
plt.ylim([-0.75, 0.75])
plt.text(0.65, 0.62, '(B)', fontsize=18)

# do area
data = []
tmplab = []
rarea = array(rarea)
for sz in regset:
    # get region index
    index = array(where(array(reg) == sz))[0]
    if sz == '':
        sz = 'OTH'
    if len(index) > 1:
        data.append(rarea[index])
        tmplab.append(sz+' '+str(len(index)))
        
ax = plt.subplot(2, 3, 3)
bp = plt.boxplot(data)
plt.setp(bp['boxes'], color='black')
plt.setp(bp['whiskers'], color='black')
plt.setp(bp['medians'], color='black')
plt.plot([0, len(tmplab)+1], [0, 0], 'k--')
plttext(ax, 0.65, 0.51, '(C)', fsize='large')
plt.xticks(range(1,len(tmplab)+1),tmplab, ha='right', rotation=45)

plt.title('Interface Rupture Area Residuals', fontsize=14)
plt.ylim([-0.75, 0.75])
plt.text(0.65, 0.62, '(C)', fontsize=18)

plt.savefig('region_bwplot.png',dpi=300,format='png', bbox_inches='tight')
plt.savefig('region_bwplot.pdf',dpi=300,format='pdf', bbox_inches='tight')
plt.show()

##############################################################################
'''
plot region for M >= 7.8
'''
fig = plt.figure(7, figsize=(15,9))
mmax = 7.8

# do length
data = []
tmplab = []
rlen = array(rlen)
for sz in regset:
    
    # get region index
    index = array( where((array(reg) == sz) & (array(mag) >= mmax) ))[0]
    if sz == '':
        sz = 'OTH'
    if len(index) > 1:
        data.append(rlen[index])
        tmplab.append(sz+' '+str(len(index)))
        
    else:
        print sz, index
        
ax = plt.subplot(2, 3, 1)
bp = plt.boxplot(data)
plt.setp(bp['boxes'], color='black')
plt.setp(bp['whiskers'], color='black')
plt.setp(bp['medians'], color='black')
plt.plot([0, len(tmplab)+1], [0, 0], 'k--')
plttext(ax, 0.65, 0.51, '(A)', fsize='large')
plt.xticks(range(1,len(tmplab)+1),tmplab, ha='right', rotation=45)

plt.title('Interface Rupture Length Residuals', fontsize=14)
plt.ylabel(r'log$\mathregular{_{10}}$(Obs / Pred)', fontsize=14)
plt.ylim([-0.6, 0.6])
plt.text(0.65, 0.5, '(A)', fontsize=18)


# do width
data = []
tmplab = []
rwid = array(rwid)
for sz in regset:
    # get region index
    index = array( where((array(reg) == sz) & (array(mag) >= mmax) ))[0]
    if sz == '':
        sz = 'OTH'
    if len(index) > 1:
        data.append(rwid[index])
        tmplab.append(sz+' '+str(len(index)))
        
ax = plt.subplot(2, 3, 2)
bp = plt.boxplot(data)
plt.setp(bp['boxes'], color='black')
plt.setp(bp['whiskers'], color='black')
plt.setp(bp['medians'], color='black')
plt.plot([0, len(tmplab)+1], [0, 0], 'k--')
plttext(ax, 0.65, 0.51, '(B)', fsize='large')
plt.xticks(range(1,len(tmplab)+1),tmplab, ha='right', rotation=45)

plt.title('Interface Rupture Width Residuals', fontsize=14)
plt.ylim([-0.6, 0.6])
plt.text(0.65, 0.5, '(B)', fontsize=18)

# do area
data = []
tmplab = []
rarea = array(rarea)
for sz in regset:
    # get region index
    index = array( where((array(reg) == sz) & (array(mag) >= mmax) ))[0]
    if sz == '':
        sz = 'OTH'
    if len(index) > 1:
        data.append(rarea[index])
        tmplab.append(sz+' '+str(len(index)))
        
ax = plt.subplot(2, 3, 3)
bp = plt.boxplot(data)
plt.setp(bp['boxes'], color='black')
plt.setp(bp['whiskers'], color='black')
plt.setp(bp['medians'], color='black')
plt.plot([0, len(tmplab)+1], [0, 0], 'k--')
plttext(ax, 0.65, 0.51, '(C)', fsize='large')
plt.xticks(range(1,len(tmplab)+1),tmplab, ha='right', rotation=45)

plt.title('Interface Rupture Area Residuals', fontsize=14)
plt.ylim([-0.6, 0.6])
plt.text(0.65, 0.5, '(C)', fontsize=18)

plt.savefig('region_bwplot_m_ge_8.pdf',dpi=300,format='pdf', bbox_inches='tight')
plt.show()

##############################################################################
fig = plt.figure(6, figsize=(5,4.5))

# plot box/whisker plots of dips per subduction
data = []
tmplab = []
noidx = []
dip = array(dip)
for sz in regset:
    # get region index
    index = array(where((array(reg) == sz) & (dip > 0)))[0]
    if len(index) > 1:
        data.append(dip[index])
        tmplab.append(sz+' '+str(len(index)))

#ax = plt.subplot(2, 3, 4)
bp = plt.boxplot(data)
plt.setp(bp['boxes'], color='black')
plt.setp(bp['whiskers'], color='black')
plt.setp(bp['medians'], color='black')
plt.plot([0, len(tmplab)+1], [median(dip), median(dip)], 'k--')
plt.xticks(range(1,len(tmplab)+1),tmplab, ha='right', rotation=45)

plt.ylabel('Average Dip (Degrees)', fontsize=14)
plt.title('Average Interface Dip', fontsize=14)
plt.ylim([0, 50])

plt.savefig('region_bwplot_dip.pdf',dpi=300,format='pdf', bbox_inches='tight')
plt.show()



##############################################################################
'''
calculate AIC
'''
print '\nAkaike Information Criteria\n'

def get_aic(res, nvar, idx):
    from numpy import pi, exp, sqrt, log
    
    resarray = array(res)[idx]
    N = len(resarray)
    SSE = sum(resarray**2)
    s2 = SSE / N
    
    L = ( 1.0/sqrt(2*pi*s2) ) ** N * exp( -SSE/(s2*2.0) )
    
    return 2 * nvar - 2 * log(L)
    
midx0 = where(array(mag) >= 0.)[0]
midx8 = where(array(mag) >= 8.)[0]

# get AIC for lin width
print 'WIDTH\n'
aic2a = get_aic(rwidlin, 2, midx0)
print 'AIC W1 = ', aic2a
aic28 = get_aic(rwidlin, 2, midx8)
print 'AIC W1 = ', aic28

# get AIC for bilin width
aic4a = get_aic(rwid, 4, midx0)
print 'AIC W2 = ', aic4a
aic48 = get_aic(rwid, 4, midx8)
print 'AIC W2 = ', aic48

# test relative likelihood 
rl = exp((aic2a - aic4a)/2.)
print 'relative likelihood (all data)', rl
rl = exp((aic28 - aic48)/2.)
print 'relative likelihood (M >= 8)', rl

# get AIC for lin area
print '\nAREA\n'
aic2a = get_aic(rarealin, 2, midx0)
print 'AIC A1 = ', aic2a
aic28 = get_aic(rarealin, 2, midx8)
print 'AIC A1 = ', aic28

# get AIC for bilin area
aic4a = get_aic(rarea, 4, midx0)
print 'AIC A2 = ', aic4a
aic48 = get_aic(rarea, 4, midx8)
print 'AIC A2 = ', aic48

# test relative likelihood 
rl = exp((aic2a - aic4a)/2.)
print 'relative likelihood (all data)', rl
rl = exp((aic28 - aic48)/2.)
print 'relative likelihood (M >= 8)', rl