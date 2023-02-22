from numpy import arange, around, array, random, delete, where, reshape, log10, sqrt, mean, floor, isnan
from os import system, path
from misc_tools import dictlist2array, get_binned_stats
from obspy import Trace
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
from os import getcwd
import scipy.odr.odrpack as odrpack
from scipy.stats import linregress
mpl.style.use('classic')

def fmt_axes_tick_decimals(ax):
    '''
    e.g., ax = plt.subplot(111)
    '''
    
    from matplotlib.ticker import FormatStrFormatter
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))    
    
# load event data
events = pickle.load(open("smsim/simulated_ml_events.pkl", "rb" ))


################################################################################r
# get arrays
################################################################################r

# get arrays for regressing
bj84_2800_array = dictlist2array(events, 'bj84_2800_mean')
mlm92_2800_array = dictlist2array(events, 'mlm92_2800_mean') 

bj84_2080_array = dictlist2array(events, 'bj84_2080_mean')
mlm92_2080_array = dictlist2array(events, 'mlm92_2080_mean') 

mw_array = dictlist2array(events, 'mw') 

################################################################################r
# get mean diff between MLs
################################################################################r

mlm92_diff = mlm92_2800_array - mlm92_2080_array
fig = plt.figure(1, figsize=(6,6))
#bins = arange(0., 0.2, 0.01)
#plt.hist(mlm92_diff, bins)
plt.plot(mlm92_2080_array, mlm92_diff, 'o', ms=6, c='0.5')
plt.ylabel('ML Sensitivity Difference')
plt.xlabel('ML (Sensitivity = 2080)')
plt.grid(which='both')

# get binned data
bins = arange(1.5, 5.5, 0.2)
medamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, mlm92_2080_array, mlm92_diff)
plt.errorbar(medx, medamp, yerr=stdbin, fmt='s', ms=8, \
             mfc='r', mec='k', ecolor='r', elinewidth=2., ls='none', zorder=2000)

# regress data
slope, intercept, r_value, p_value, std_err  = linregress(medx, medamp)

# fit data
xplt = array([1., 5.5])
yplt = slope * xplt + intercept
plt.plot(xplt, yplt, 'k-', lw=2.5)

plt.savefig('ml_sensitivity_diff.png',format='png', dpi=300, bbox_inches='tight')
plt.show()

# write regression coeffs
outtxt = '# coeffs to correct ML from W-A sensitivity from 2080 to 2800, e.g.: ML2800 = c0 * ML2080 + c1\n'
outtxt += ','.join((str(slope), str(intercept)))
f = open('wa_sensitivity_coeffs.csv', 'w')
f.write(outtxt)
f.close()

################################################################################r
# plot MW vs MW
################################################################################r
fig = plt.figure(2, figsize=(8,8))
ax = plt.subplot(111)

plt.plot([1., 6.], [1., 6.], 'k--', label='1:1')
plt.plot(mlm92_2800_array, mw_array, 'o', ms=7, c='dodgerblue', label='W-A Sensitivity = 2800')
#plt.plot(mlm92_2080_array, mw_array, 'o', ms=7, c='orange', label='W-A Sensitivity = 2080') 

###############################################################################
# regress 2800
###############################################################################

def ortho_quadratic_reg(c, x):
    return c[0] * x**2 + c[1] * x + c[2]


idx = where((isnan(mlm92_2800_array) == False) & (isnan(mw_array) == False))[0] 
data = odrpack.RealData(mlm92_2800_array[idx], mw_array[idx])
'''
sx = np.interp(mw, [min(mw), max(mw)],[0.3, 0.1])
sy = sx
data = odrpack.RealData(mb[mb>=3.5], mw[mb>=3.5], sx=sx[mb>=3.5], sy=sy[mb>=3.5])
'''
quad_reg = odrpack.Model(ortho_quadratic_reg)
odr = odrpack.ODR(data, quad_reg, beta0=[0.05, 0.5, 1.5])

odr.set_job(fit_type=0) #if set fit_type=2, returns the same as least squares
out = odr.run()
c0 = out.beta[0]
c1 = out.beta[1]
c2 = out.beta[2]

xplt = arange(1, 6., 0.01)
yplt = c0 * xplt**2 + c1 * xplt + c2
plt.plot(xplt, yplt, 'k-', lw=2)

print('\nWeight by probability - but use sample of 10000!!!!\n')
################################################################################
# make pretty
plt.xlabel('ML', fontsize=20)
plt.ylabel('MW', fontsize=20)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.grid(which='both')
plt.legend(loc=2, numpoints=3)
plt.xlim([1.75, 7.0])
plt.ylim([1.75, 7.0])
fmt_axes_tick_decimals(ax)

plt.savefig('ml2mw_reg_var_wa_resp.png',format='png', dpi=300, bbox_inches='tight')
plt.show()
