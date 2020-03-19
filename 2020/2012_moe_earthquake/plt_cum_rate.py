import matplotlib.pyplot as plt
from numpy import arange, mean, diff, percentile, array, unique, where, argsort, vstack
from datetime import datetime, timedelta
from catalogue_tools import aki_maximum_likelihood, fit_a_value  
from misc_tools import get_log_xy_locs
import matplotlib 
matplotlib.rc('xtick', labelsize=16) 
matplotlib.rc('ytick', labelsize=16) 

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
ax = plt.subplot(1, 2, 1)

mrng = arange(-0.1, 4.3, 0.1)
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
    if m == 4.2:
        m4idx = idx[0]

# get b-value
mc = 0.9
b_val, sigma_b = aki_maximum_likelihood(mrng, number_obs, mc)

# get a-value
m_upper = mc+1.2 # only include well behaced data
a_val = fit_a_value(cum_mag, b_val, mrng, mc, m_upper)
N = 10**(a_val - b_val*mrng)

plt.semilogy(mrng, cum_mag, 'ko', ms=7)
plt.semilogy(mrng[9:], N[9:], 'r-', lw=2.)
# get xlim
#plt.xlim([datetime(2012,6,15), datetime(2012,12,31)])
plt.ylabel('Cumulative Number', fontsize=18)
plt.xlabel('Local Magnitude', fontsize=18)
plt.xlim([-0.2, 4.5])
plt.ylim([0.5, 1000])

# label b-value
abtxt = 'a-value = '+str('%0.2f' % a_val)+'\nb-value = '+str('%0.2f' % b_val)
xtxt = ax.get_xlim()[1] * 0.98
ytxt = get_log_xy_locs(ax.get_ylim(), 0.98)
plt.text(xtxt, ytxt, abtxt, size=18, horizontalalignment='right', verticalalignment='top', weight='normal')

xdiff = diff(ax.get_xlim())
xtxt = ax.get_xlim()[0] + xdiff * 0.02
ytxt = get_log_xy_locs(ax.get_ylim(), 0.98)
plt.text(xtxt, ytxt, 'A', size=20, ha='left', verticalalignment='top', weight='normal')

plt.grid(b=True, which='both',axis='both')

ax = plt.subplot(1, 2, 2)
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
            arrowprops=dict(fc='k', shrink=0.05,width=2, ec='k'), fontsize=15)

ax.annotate('2012-07-20 ML 4.2', xy=([evdt[m4idx], m4idx+1]), xytext=(datetime(2012,8,10), 195), color='k', \
            arrowprops=dict(fc='k', shrink=0.05,width=2, ec='k'), fontsize=15)
                                    
xtxt = get_log_xy_locs(ax.get_xlim(), 0.02)            
ytxt = ax.get_ylim()[1] * 0.98
plt.text(xtxt, ytxt, 'B', size=20, ha='left', verticalalignment='top', weight='normal')


plt.savefig('moe_cum_rate.png', format='png', bbox_inches='tight', dpi=300)
plt.show()
