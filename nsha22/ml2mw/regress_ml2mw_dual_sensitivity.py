from numpy import arange, around, array, random, delete, where, reshape, log10, sqrt, \
                  mean, floor, isnan, polyfit, poly1d, zeros_like
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

'''
loads pkl made by smsim/calculate_ml2mw_dual_sensitivity.py
'''

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

def fit_parabola(c, x):
    return c[0] * (x+c[2])**3 + c[1]


mlm92_diff = mlm92_2800_array - mlm92_2080_array
fig = plt.figure(1, figsize=(6,6))
#bins = arange(0., 0.2, 0.01)
#plt.hist(mlm92_diff, bins)
plt.plot(mlm92_2080_array, mlm92_diff, 'o', ms=6, c='0.5')
plt.ylabel('Wood-Anderson Sensitivity ML Difference', fontsize=14)
plt.xlabel('ML (Sensitivity = 2080)', fontsize=14)
plt.grid(which='both')

# get binned data
bins = arange(1.5, 6.6, 0.2)
medamp, stdbin, medx, binstrp, nperbin = get_binned_stats(bins, mlm92_2080_array, mlm92_diff)
plt.errorbar(medx, medamp, yerr=stdbin, fmt='s', ms=8, \
             mfc='r', mec='k', ecolor='r', elinewidth=2., ls='none', zorder=2000)

# regress linear data
slope, intercept, r_value, p_value, std_err  = linregress(medx, medamp)

# plot data
xplt = array([1., 6.5])
yplt = slope * xplt + intercept
#plt.plot(xplt, yplt, 'k-', lw=2.5)

# try polyfit porabola
x_norm = (medx - medx.mean())/medx.std()
fit_normalized = polyfit(x_norm, medamp, 4)
f = poly1d(fit_normalized)
plt.plot(medx, f(x_norm), 'k-', lw=2.5, zorder=10000)

'''
fn = fit_normalized
yplt = fn[0] * x_norm**4 + fn[1] * x_norm**3 + fn[2] * x_norm**2 + fn[3] * x_norm + fn[4]
plt.plot(medx, yplt,'m-', lw=10)
'''
"""
# now do parabola
idx = where((isnan(medx) == False) & (isnan(medamp) == False))[0] 
data = odrpack.RealData(x_norm[idx], medamp[idx])

test_txt = ''
for x,y in zip(medx, medamp):
	test_txt += ','.join((str(x), str(y))) + '\n'
	
f = open('par_fit.csv', 'w')
f.write(test_txt)
f.close()

par_reg = odrpack.Model(fit_parabola)
odr = odrpack.ODR(data, par_reg, beta0=[.002, 0.87, 0.0])

odr.set_job(fit_type=2) #if set fit_type=2, returns the same as least squares
out = odr.run()
p0 = out.beta[0]
p1 = out.beta[1]
p2 = out.beta[2]

yplt = p0 * (x_norm+p2)**3 + p1
plt.plot(x_norm, yplt, 'b--', lw=2.5)
"""

plt.savefig('ml_sensitivity_diff.png',format='png', dpi=300, bbox_inches='tight')	
plt.show()

# write regression coeffs
outtxt = '# coeffs to correct ML from W-A sensitivity from 2080 to 2800, e.g.: ML2800 = c0**4 * ML2080 ... + c4\n'
outtxt += ','.join([str(x) for x in fit_normalized]) + ','
outtxt += ','.join((str(medx.mean()), str(medx.std())))
f = open('wa_sensitivity_coeffs.csv', 'w')
f.write(outtxt)
f.close()

################################################################################r
# plot MW vs ML - 2800
################################################################################r
fig = plt.figure(2, figsize=(8,8))
ax = plt.subplot(111)

plt.plot([1., 7], [1., 7], 'k--', label='1:1')
plt.plot(mlm92_2800_array, mw_array, 'o', ms=7, c='dodgerblue', label='Simulated Data (Sensitivity = 2800)')
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

xplt = arange(1, 6.61, 0.01)
yplt = c0 * xplt**2 + c1 * xplt + c2
plt.plot(xplt, yplt, '-', c='orangered', lw=2.5, label='Qudratic Fit')

# write regression coeffs
outtxt = '# coeffs for MW-ML for W-A sensitivity = 2800, e.g.: MW = c0 * ML2800**2 + c1 * ML2800 + c2\n'
outtxt += ','.join((str(c0), str(c1), str(c2)))
f = open('mw-ml_coeffs_2800.csv', 'w')
f.write(outtxt)
f.close()

print('\nWeight by probability - but use sample of 10000!!!!\n')
################################################################################
# make pretty
plt.xlabel('ML', fontsize=20)
plt.ylabel('MW', fontsize=20)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.grid(which='both')
plt.legend(loc=2, numpoints=1)
plt.xlim([1.75, 7.0])
plt.ylim([1.75, 7.0])
fmt_axes_tick_decimals(ax)

plt.savefig('ml2mw_reg_2800.png',format='png', dpi=300, bbox_inches='tight')
plt.show()

################################################################################r
# plot MW vs ML - 2080
################################################################################r
fig = plt.figure(2, figsize=(8,8))
ax = plt.subplot(111)

plt.plot([1., 7], [1., 7], 'k--', label='1:1')
plt.plot(mlm92_2080_array, mw_array, 'o', ms=7, c='dodgerblue', label='Simulated Data (Sensitivity = 2080)')
#plt.plot(mlm92_2080_array, mw_array, 'o', ms=7, c='orange', label='W-A Sensitivity = 2080') 

###############################################################################
# regress 2080
###############################################################################

def ortho_quadratic_reg(c, x):
    return c[0] * x**2 + c[1] * x + c[2]
    
idx = where((isnan(mlm92_2080_array) == False) & (isnan(mw_array) == False))[0] 
data = odrpack.RealData(mlm92_2080_array[idx], mw_array[idx])
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

xplt = arange(1, 6.61, 0.01)
yplt = c0 * xplt**2 + c1 * xplt + c2
plt.plot(xplt, yplt, '-', c='orangered', lw=2.5, label='Qudratic Fit')

# write regression coeffs
outtxt = '# coeffs for MW-ML for W-A sensitivity = 2080, e.g.: MW = c0 * ML2080**2 + c1 * ML2080 + c2\n'
outtxt += ','.join((str(c0), str(c1), str(c2)))
f = open('mw-ml_coeffs_2080.csv', 'w')
f.write(outtxt)
f.close()

################################################################################
# get lin-quad
################################################################################
def lowside(x, hx):
    from numpy import zeros_like
    xmod = zeros_like(x)

    idx = x <= hx
    xmod[idx] = 1
    return xmod

def highside(x, hx):
    from numpy import zeros_like
    xmod = zeros_like(x)

    idx = x > hx
    xmod[idx] = 1
    return xmod
    
def lin2quad_reg(c, x):
   from numpy import zeros_like
   hx = c[4] # hinge magnitude
   ans2 = zeros_like(x)
   ans1 = zeros_like(x)

   modx_lo = lowside(x, hx)
   modx_hi = highside(x, hx)
   
   ans1 = modx_lo * (c[0] * x + c[1])
   yhinge = c[0] * hx + c[1]
   ans2 = modx_hi * (c[2]*(x-hx) + c[3]*(x-hx)**2 + yhinge)

   return ans1 + ans2

l2q_reg = odrpack.Model(lin2quad_reg)
odr = odrpack.ODR(data, l2q_reg, beta0=[0.6, 1.0, 0.7, 0.1, 4.5])

odr.set_job(fit_type=0) #if set fit_type=2, returns the same as least squares
out = odr.run()
lqc0 = out.beta[0]
lqc1 = out.beta[1]
lqc2 = out.beta[2]
lqc3 = out.beta[3]
lqc4 = out.beta[4]
print('\nLin to Quad')
out.pprint()

# write regression coeffs
header = 'For ML <= c4: MW = c0 * ML + c1\n'
header += 'For ML > c4: MW = c0 * c4 + c1 + c2 * (ML-c4) + c3 * (ML-c4)**2\n'
outtxt = header + ','.join((str(lqc0), str(lqc1), str(lqc2), str(lqc3), str(lqc4)))
f = open('mw-ml_lin2quad_coeffs_2080.csv', 'w')
f.write(outtxt)
f.close()

yplt = zeros_like(xplt)
idx = xplt <= lqc4
yplt[idx] = lqc0 * xplt[idx] + lqc1
idx = xplt > lqc4
yplt[idx] = lqc0 * lqc4 + lqc1 \
            + lqc2 * (xplt[idx]-lqc4) + lqc3 * (xplt[idx]-lqc4)**2
plt.plot(xplt, yplt, '-', lw=2.5, c='limegreen', label='Linear-Quadratic')


print('\nWeight by probability - but use sample of 10000!!!!\n')
################################################################################
# make pretty
plt.xlabel('ML', fontsize=20)
plt.ylabel('MW', fontsize=20)
plt.xticks(fontsize=17)
plt.yticks(fontsize=17)
plt.grid(which='both')
plt.legend(loc=2, numpoints=1)
plt.xlim([1.5, 7.0])
plt.ylim([1.5, 7.0])
fmt_axes_tick_decimals(ax)

plt.savefig('ml2mw_reg_2080.png',format='png', dpi=300, bbox_inches='tight')
plt.show()
