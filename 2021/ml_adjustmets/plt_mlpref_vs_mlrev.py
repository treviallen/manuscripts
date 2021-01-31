import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('classic')
from numpy import loadtxt, arange, array
from scipy.stats import linregress
from misc_tools import get_binned_stats
from scipy.odr import RealData, Model, ODR, models
from os import getcwd

if getcwd().startswith('/nas'):
    magdiffcsv = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/catalogue/matlab/rev_mag_diff.csv'
else:
    magdiffcsv = '/Users/trev/Documents/Geoscience_Australia/NSHA2018/catalogue/matlab/rev_mag_diff.csv'

data = loadtxt(magdiffcsv, delimiter=',')

fig = plt.figure(1, figsize=(7,7))

plt.plot([2,7], [2,7], '--', c='0.25', lw=2.0, label='1:1')
plt.plot(data[:,2], data[:,3], '+', c='r', label='Data')
	
# bin data and regress
mbins = arange(1, 6.2, 0.2)
medres, stdres, medx, mbins, nperbin = get_binned_stats(mbins, data[:,2], data[:,3])
#plt.errorbar(mbins, medres, stdres, marker='s', c='orange', mec='orangered', ls='none', ms=8, label='Binned Data')
reg_pars = linregress(mbins, medres)
xplt = array([2.4, 6.5])
yplt = reg_pars[0] * xplt + reg_pars[1]
#plt.plot(xplt, yplt, 'r-', lw=2., zorder=100, label='Linear Fit')

#reg_pars2 = linregress(data[:,2], data[:,3])

# now do ODR
# Define a function (quadratic in our case) to fit the data with.
def linear_func(p, x):
   m, c = p
   return m*x + c

# Create a model for fitting.
linear_model = Model(linear_func)
# Create a RealData object using our initiated data from above.
data = RealData(data[:,2], data[:,3])
# Set up ODR with the model and data.
odr = ODR(data, linear_model, beta0=[1., 0.])
odr.set_job(fit_type=0) #if set fit_type=2, returns the same as leastsq
# Run the regression.
out = odr.run()
# Use the in-built pprint method to give us results.
out.pprint()
modr = out.beta[0]
codr = out.beta[1]
yplt_odr = modr * xplt + codr
plt.plot(xplt, yplt_odr, 'b-', lw=2, label='ODR')

plt.xlim([1.8,6.5])
plt.ylim([1.8,6.5])

plt.legend(numpoints=1, loc=2)
plt.grid(which='both')
plt.xlabel('Original $\mathregular{M_L}$ $\mathregular{(M_{LH})}$', fontsize=16)
plt.ylabel('Revised $\mathregular{M_L}$ $\mathregular{(M_{LR})}$', fontsize=16)

plt.savefig('hist_mag_fit.png', fmt='png', bbox_inches='tight')
plt.savefig('hist_mag_fit.eps', fmt='eps', bbox_inches='tight')

plt.show()