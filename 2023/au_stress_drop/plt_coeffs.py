import pickle
import matplotlib.pyplot as plt
from misc_tools import dictlist2array, savitzky_golay
from numpy import mean, median, where, array, log10, nanmean

# load atten coeffs
coeffs = pickle.load(open('atten_coeffs_1.3_5km.pkl', 'rb' ))

freqs = dictlist2array(coeffs, 'freq')

# plt nc0
nc0 = dictlist2array(coeffs, 'nc0')
nc0s = dictlist2array(coeffs, 'nc0s')

fig = plt.figure(1, figsize=(14, 10))

plt.subplot(421)
plt.semilogx(freqs, nc0, 'ro')
plt.semilogx(freqs, nc0s, 'bo')
plt.ylim([-2,0])
plt.ylabel('nc0')

#################################################################################
plt.subplot(422)
nc1 = dictlist2array(coeffs, 'nc1')
nc1s = dictlist2array(coeffs, 'nc1s')
nc1f = dictlist2array(coeffs, 'nc1f')

plt.semilogx(freqs, nc1, 'ro')
plt.semilogx(freqs, nc1s, 'bo')
plt.semilogx(freqs, nc1f, 'co')
plt.ylabel('nc1')

# fit quadratic
from numpy import polyfit
qc = polyfit(log10(freqs), nc1s, 2)
yplt = qc[0]*log10(freqs)**2 + qc[1]*log10(freqs) + qc[2]
plt.semilogx(freqs, yplt, 'g-')

# fit cubic
qc = polyfit(log10(freqs[15:]), nc1s[15:], 3)
yplt = qc[0]*log10(freqs)**3 + qc[1]*log10(freqs)**2 + qc[2]*log10(freqs) + qc[3]
plt.semilogx(freqs, yplt, 'k--')

# fit exponential
# fit mid
def fit_exponential(c, x):
    from numpy import sqrt, log10, exp

    ans = c[0]**log10(x) + c[1]

    return ans

from scipy.odr import Data, Model, ODR, models
import scipy.odr.odrpack as odrpack

data = odrpack.RealData(freqs, nc1s)
    
afit = odrpack.Model(fit_exponential)
odr = odrpack.ODR(data, afit, beta0=[1, 2])

odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
out = odr.run()
fn1e = out.beta

yplt = fn1e[0]**log10(freqs) + fn1e[1]
plt.semilogx(freqs, yplt, 'm-', lw=2)	

#################################################################################

# plt mc0
mc0 = dictlist2array(coeffs, 'mc0')
mc0s = dictlist2array(coeffs, 'mc0s')
mc0t = dictlist2array(coeffs, 'mc0f')
plt.subplot(423)
plt.semilogx(freqs, mc0, 'ro')
plt.semilogx(freqs, mc0s, 'bo')
plt.semilogx(freqs, mc0t, 'co')
print(median(mc0))
plt.ylabel('mc0')
#nc0s = dictlist2array(coeffs, 'nc0s')

# test mc0 fitting 
from scipy.stats import linregress

# fit to 3.5 hz
f0 = 0.03
f1 = 30
f2 = 8.
f3 = 20
idx = where((freqs > f0) & (freqs <= f1))[0]

reg1 = linregress(log10(freqs[idx]), mc0s[idx])
xplt = array([0.03162278, f1])
yplt = reg1.intercept + reg1.slope*log10(xplt)
plt.plot(xplt, yplt, 'g', label='fitted line')

# fit mid
def fit_mid_mc0(c, x):
    from numpy import sqrt, log10

    # set low-f
    ans = reg1.intercept + reg1.slope*log10(x)

    # set mid-f
    idx = where((x > f1) & (x <= f2))[0]  
    ans[idx] = reg1.intercept + reg1.slope*log10(f1) + c[0] * (log10(x[idx]) - log10(f1)) \

    return ans

'''
data = odrpack.RealData(freqs, mc0s)
    
afit = odrpack.Model(fit_mid_mc0)
odr = odrpack.ODR(data, afit, beta0=[-2])

odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
out = odr.run()
fcm = out.beta

xplt = array([f1, f2])
yplt = reg1.intercept + reg1.slope*log10(f1) + fcm[0] * (log10(xplt) - log10(f1))
plt.plot(xplt, yplt, 'g')
'''
# fit high
def fit_high_mc0(c, x):
    from numpy import sqrt, log10

    # set far-f
    ans = reg1.intercept + reg1.slope*log10(f1) + fcm[0] * (log10(f2) - log10(f1)) \
            + c[0] * (log10(x) - log10(f2))

    return ans

'''
idx = where((freqs > f2) & (freqs <= f3))[0]  
data = odrpack.RealData(freqs[idx], mc0s[idx])
    
afit = odrpack.Model(fit_far_mc0)
odr = odrpack.ODR(data, afit, beta0=[0.5])

odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
out = odr.run()
fcf = out.beta

xplt = array([f2, max(freqs)])
yplt = reg1.intercept + reg1.slope*log10(f1) + fcm[0] * (log10(f2) - log10(f1)) \
       + fcf[0] * (log10(xplt) - log10(f2))
plt.plot(xplt, yplt, 'g')
'''
##################################################################
# plt mc1
mc1 = dictlist2array(coeffs, 'mc1')
mc1s = dictlist2array(coeffs, 'mc1s')
mc1t = dictlist2array(coeffs, 'mc1f')
mc1h = dictlist2array(coeffs, 'mc1h')

plt.subplot(424)
plt.semilogx(freqs, mc1, 'ro')
plt.semilogx(freqs, mc1s, 'bo')
plt.semilogx(freqs, mc1t, 'co')
#plt.semilogx(freqs, mc1ts, 'mo')
plt.ylabel('mc1')

f0 = 0.03
f1 = 0.2
f2 = 0.5
f3 = 30
idx = where((freqs > f0) & (freqs <= f1))[0]

reg1 = linregress(log10(freqs[idx]), mc1s[idx])
xplt = array([0.03162278, f1])
yplt = reg1.intercept + reg1.slope*log10(xplt)
plt.plot(xplt, yplt, 'g', label='fitted line')

# fit mid
def fit_mid_mc0(c, x):
    from numpy import sqrt, log10

    # set mid-f
    ans = reg1.intercept + reg1.slope*log10(f1) + c[0] * (log10(x) - log10(f1))

    return ans

idx = where((freqs > f1) & (freqs <= f2))[0]
data = odrpack.RealData(freqs[idx], mc1s[idx])
    
afit = odrpack.Model(fit_mid_mc0)
odr = odrpack.ODR(data, afit, beta0=[-0.])

odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
out = odr.run()
fcm = out.beta

xplt = array([f1, f2])
yplt = reg1.intercept + reg1.slope*log10(f1) + fcm[0] * (log10(xplt) - log10(f1))
plt.plot(xplt, yplt, 'g')

# take mean < 2 Hz then smooth raw
idx = where(freqs < 2.0)[0]
mean_m1 = nanmean(mc1[idx])
hybrid_mc1 = mc1
hybrid_mc1[idx] = mean_m1

# smooth
sg_window = 41
sg_poly = 3
#smooth_hybrid_mc1 = savitzky_golay(hybrid_mc1, sg_window, sg_poly) # mid-slope
plt.plot(freqs, mc1h, 'mo')

# fit high
def fit_high_mc0(c, x):
    from numpy import sqrt, log10

    # set far-f
    ans = reg1.intercept + reg1.slope*log10(f1) + fcm[0] * (log10(f2) - log10(f1)) \
          + c[0] * (log10(x) - log10(f2))

    return ans

idx = where((freqs > f2) & (freqs <= f3))[0]  
data = odrpack.RealData(freqs[idx], mc1s[idx])
    
afit = odrpack.Model(fit_high_mc0)
odr = odrpack.ODR(data, afit, beta0=[0.0])

odr.set_job(fit_type=2) #if set fit_type=2, returns the same as leastsq, 0=ODR
out = odr.run()
fcf = out.beta

xplt = array([f2, max(freqs)])
yplt = reg1.intercept + reg1.slope*log10(f1) + fcm[0] * (log10(f2) - log10(f1)) \
       + fcf[0] * (log10(xplt) - log10(f2))
plt.plot(xplt, yplt, 'g')

############################################################

plt.subplot(425)
fc0 = dictlist2array(coeffs, 'fc0')
fc0s = dictlist2array(coeffs, 'fc0s')
plt.semilogx(freqs, fc0, 'ro')
plt.semilogx(freqs, fc0s, 'bo')
#plt.ylim([-2,0])
plt.ylabel('fc0')

plt.subplot(426)
fc1 = dictlist2array(coeffs, 'fc1')
fc1s = dictlist2array(coeffs, 'fc1s')
plt.semilogx(freqs, fc1, 'ro')
plt.semilogx(freqs, fc1s, 'bo')
#plt.ylim([-2,0])
plt.ylabel('fc1')
'''
plt.subplot(427)
fc2 = dictlist2array(coeffs, 'fc2')
#fc2s = dictlist2array(coeffs, 'fc2s')
plt.semilogx(freqs, fc2, 'ro')
#plt.semilogx(freqs, fc2s, 'bo')
#plt.ylim([-2,0])
plt.ylabel('fc2')

#fig = plt.figure(1, figsize=(12, 10))
'''
'''
#plt.semilogx(freqs, nc0s, 'bo')

# plt mc0
mc1 = dictlist2array(coeffs, 'mc1')
#nc0s = dictlist2array(coeffs, 'nc0s')


fig = plt.figure(1, figsize=(12, 9))

plt.subplot(313)
plt.semilogx(freqs, mc1, 'ro')
idx = where((freqs > 0.5) & (freqs < 10))[0]
print(median(mc1[idx]))
'''

plt.subplot(427)
fc1 = dictlist2array(coeffs, 'magc0')
fc1s = dictlist2array(coeffs, 'magc0s')
plt.semilogx(freqs, fc1, 'ro')
plt.semilogx(freqs, fc1s, 'bo')
plt.ylabel('mag0')

plt.subplot(428)
fc1 = dictlist2array(coeffs, 'magc1')
fc1s = dictlist2array(coeffs, 'magc1s')
plt.semilogx(freqs, fc1, 'ro')
plt.semilogx(freqs, fc1s, 'bo')
plt.ylabel('mag1')


plt.show()

######################################
# Plot regional FF corrections

fig = plt.figure(1, figsize=(14, 10))

NCCZ = dictlist2array(coeffs, 'NCCZ_rc')
EBGZ = dictlist2array(coeffs, 'EBGZ_rc')
CBGZ = dictlist2array(coeffs, 'CBGZ_rc')

plt.subplot(421)
plt.semilogx(freqs, NCCZ, 'ro')
#plt.ylim([-2,0])
plt.ylabel('NCCZ')

plt.subplot(422)
plt.semilogx(freqs, EBGZ, 'ro')
#plt.ylim([-2,0])
plt.ylabel('EBGZ')

plt.subplot(423)
plt.semilogx(freqs, CBGZ, 'ro')
#plt.ylim([-2,0])
plt.ylabel('CBGZ')

plt.show()
