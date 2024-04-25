import pickle
import matplotlib.pyplot as plt
from misc_tools import dictlist2array
from numpy import mean, median, where, array, log10


# load atten coeffs
coeffs = pickle.load(open('atten_coeffs.pkl', 'rb' ))

freqs = dictlist2array(coeffs, 'freq')

# plt nc0
nc0 = dictlist2array(coeffs, 'nc0')
nc0s = dictlist2array(coeffs, 'nc0s')
nc1 = dictlist2array(coeffs, 'nc1s')

fig = plt.figure(1, figsize=(14, 10))

plt.subplot(421)
plt.semilogx(freqs, nc0, 'ro')
plt.semilogx(freqs, nc0s, 'bo')
plt.ylim([-2,0])
plt.ylabel('nc0')

plt.subplot(422)
plt.semilogx(freqs, nc1, 'ro')
plt.ylabel('nc1')

# plt mc0
mc0 = dictlist2array(coeffs, 'mc0')
mc0s = dictlist2array(coeffs, 'mc0s')
mc0t = dictlist2array(coeffs, 'mc0t')
plt.subplot(423)
plt.semilogx(freqs, mc0, 'ro')
plt.semilogx(freqs, mc0s, 'bo')
plt.semilogx(freqs, mc0t, 'co')
print(median(mc0))
plt.ylabel('mc0')
#nc0s = dictlist2array(coeffs, 'nc0s')

# test mc0 fitting 
from scipy.stats import linregress
from scipy.odr import Data, Model, ODR, models
import scipy.odr.odrpack as odrpack

# fit to 3.5 hz
f0 = 0.03
f1 = 15
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

# plt mc1
mc1 = dictlist2array(coeffs, 'mc1')
mc1s = dictlist2array(coeffs, 'mc1s')
plt.subplot(424)
plt.semilogx(freqs, mc1, 'ro')
plt.semilogx(freqs, mc1s, 'bo')
plt.ylabel('mc1')

plt.subplot(425)
fc0 = dictlist2array(coeffs, 'fc0')
plt.semilogx(freqs, fc0, 'ro')
#plt.semilogx(freqs, nc0s, 'bo')
#plt.ylim([-2,0])
plt.ylabel('fc0')

plt.subplot(426)
fc1 = dictlist2array(coeffs, 'fc1')
plt.semilogx(freqs, fc1, 'ro')
#plt.semilogx(freqs, nc0s, 'bo')
#plt.ylim([-2,0])
plt.ylabel('fc1')

plt.subplot(427)
fc2 = dictlist2array(coeffs, 'fc2')
plt.semilogx(freqs, fc2, 'ro')
#plt.semilogx(freqs, nc0s, 'bo')
#plt.ylim([-2,0])
plt.ylabel('fc2')

#fig = plt.figure(1, figsize=(12, 10))

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
plt.show()
