import matplotlib.pyplot as plt
from numpy import arange, mean, diff, percentile, array, unique, where, argsort, vstack, unique
from datetime import datetime, timedelta
from catalogue_tools import aki_maximum_likelihood, fit_a_value  
from misc_tools import get_log_xy_locs
import matplotlib 
matplotlib.rc('xtick', labelsize=16) 
matplotlib.rc('ytick', labelsize=16) 

import matplotlib as mpl     
mpl.style.use('classic')


plt.rcParams['pdf.fonttype'] = 42

##########################################################################################
# parse epicentres
##########################################################################################
epifile = 'Moe_Events.csv'
lines = open(epifile).readlines()[2:] # don't read ML 5.4!
evdt = []
evla = []
evlo = []
evdp = []
evml = []
for line in lines:
    dat = line.split(',')
    evdt.append(datetime.strptime(dat[0][0:-2], '%Y-%m-%d %H%M %S'))
    evla.append(float(dat[1]))
    evlo.append(float(dat[2]))
    evdp.append(float(dat[3]))
    evml.append(float(dat[5]))

##########################################################################################
# get MFD
##########################################################################################
fig = plt.figure(1, figsize=(16, 7.5))
ax1 = plt.subplot(1, 2, 2)

mrng = arange(-0.1, 4.5, 0.1)
mbin2 = diff(mrng)[0] / 2.
cum_mag = []
number_obs = []
for m in mrng:
    idx = where(array(evml) >= m)[0]
    cum_mag.append(len(idx))
    
    # get events per bin
    idx = where((array(evml) >= m-mbin2) & (array(evml) < m+mbin2))[0]
    number_obs.append(len(idx))
    
    # get aftershock idx
    if m == 4.4:
        m4idx = idx[0]

uidx = unique(cum_mag[::-1], return_index=True, return_inverse=True)[1]                

cum_diff = [cum_mag[0]]
for i in range(0, len(cum_mag[:-1])):
    cum_diff.append(cum_mag[i] - cum_mag[i+1])
    
nidx = where(array(cum_diff) != 0)[0]

# get b-values
bvals = []
sigs_b = []
for mc in mrng:
    #mc = 0.9
    b_val, sigma_b = aki_maximum_likelihood(mrng, number_obs, mc)
    bvals.append(b_val)
    sigs_b.append(sigma_b)

#plt.semilogy(mrng[nidx], array(cum_mag)[nidx], 'ko', ms=7)

ax1.set_xlabel('Local Magnitude', fontsize=18)
ax1.set_ylabel('Cumulative Number', fontsize=18)

h1 = ax1.semilogy(mrng[::-1][uidx], array(cum_mag)[::-1][uidx], 'ko', ms=7, label='Cumulative')

ax2 = ax1.twinx()
h2 = ax2.errorbar(mrng[:-11], array(bvals)[:-11],yerr=sigs_b[:-11], marker='^', ls='-', c='0.5', ms=7, label='b-value')
ax2.set_ylabel('b-value', fontsize=18)

plt.legend([h1[0], h2[0]], ['Cumulative', 'b-value'], loc=3, fontsize=16, numpoints=1)

#plt.semilogy(mrng[::-1][uidx], array(cum_mag)[::-1][uidx], 'ko', ms=7)
# get xlim
#plt.xlim([datetime(2012,6,15), datetime(2012,12,31)])    

# take b-val at ML 0.9 and get a-value
mc = 1.0
m_upper = mc+1. # only include well behaved data
a_val = fit_a_value(cum_mag, bvals[10], mrng, mc, m_upper)
N = 10**(a_val - bvals[10]*mrng)
ax1.semilogy(mrng[10:], N[10:], 'r-', lw=2.)
ax1.set_xlim([-0.2, 4.5])
ax1.set_ylim([0.5, 1000])

# label b-value
abtxt = 'a-value = '+str('%0.2f' % a_val)+'\nb-value = '+str('%0.2f' % bvals[10])
xtxt = 1. #ax1.get_xlim()[1] * 0.96
ytxt = get_log_xy_locs(ax1.get_ylim(), 0.96)
props = dict(boxstyle='round', facecolor='w', alpha=1)
tb = ax1.text(xtxt, ytxt, abtxt, size=18, ha='left', va='top', weight='normal', bbox=props)

xdiff = diff(ax1.get_xlim())[0]
xtxt = ax1.get_xlim()[0] + xdiff * 0.02
ytxt = get_log_xy_locs(ax1.get_ylim(), 0.98)
ax1.text(xtxt, ytxt, '(b)', size=20, ha='left', verticalalignment='top', weight='normal')

ax1.grid(b=True, which='both',axis='both')

tb.set_zorder(1000)

##########################################################################################
# make subplot A
##########################################################################################

ax = plt.subplot(1, 2, 1)
n_bins = len(evdt)
nbins = 1000
plt.hist(evdt, n_bins, histtype='step', ec='k', cumulative=True, lw=1.5)
plt.xlim([datetime(2012,6,10), evdt[-1]])
plt.grid(True, which='major',axis='both')

# relabel x-axis
ticks = arange(0, 1.1, 1./6.)
months = range(7, 13)
labels = []
xtick = []
for m in months:
    tick = datetime(2012, m, 1)
    labels.append(tick.strftime('%d-%b-%y')[1:])
    xtick.append(tick)

td = timedelta(days=3)    
plt.xticks(xtick, labels, rotation=30, ha='right')
ax.xaxis.labelpad = 5

# annotate
ax.annotate('2012-06-19 ML 5.4', xy=([evdt[0], 0]), xytext=(datetime(2012,7,10), 25), color='k', \
            arrowprops=dict(fc='k', shrink=0.05,width=2, ec='k'), fontsize=15, bbox=props)

ax.annotate('2012-07-20 ML 4.4', xy=([evdt[m4idx], m4idx+1]), xytext=(datetime(2012,8,10), 195), color='k', \
            arrowprops=dict(fc='k', shrink=0.05,width=2, ec='k'), fontsize=15, bbox=props)
                                    
xtxt = get_log_xy_locs(ax.get_xlim(), 0.02)            
ytxt = ax.get_ylim()[1] * 0.98
plt.text(xtxt, ytxt, '(a)', size=20, ha='left', verticalalignment='top', weight='normal')
plt.ylabel('Cumulative Number', fontsize=18)

plt.savefig('moe_cum_rate.png', format='png', bbox_inches='tight', dpi=600)
plt.show()
