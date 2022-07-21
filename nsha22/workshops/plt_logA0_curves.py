from calculate_magnitudes import *
from numpy import logspace, sqrt, arange
import matplotlib.pyplot as plt
import matplotlib as mpl
from misc_tools import get_mpl2_colourlist
mpl.style.use('classic')

plt.figure(1, figsize=(8, 5))
cols = get_mpl2_colourlist()

dep = 10.
repi = arange(1, 1000, 10)
rhyp = sqrt(repi**2 + dep**2)

r35 = []
bj84 = []
gs86 = []
gg91 = []
mlm92 = []
w96 = []

for re, rh in zip(repi, rhyp):
	h = 1
	v = 0
	logA = 0.
	# do not worry about V/H correction
	r35.append(calc_R35(h, logA, re))
	bj84.append(calc_BJ84(h, logA, re))
	gs86.append(calc_GS86(v, logA, re)) # native V - but no correction
	gg91.append(calc_GG91(h, logA, re))
	mlm92.append(calc_MLM92(h, logA, re))
	w96.append(calc_WGW96(v, logA, re))

# now plot
plt.plot(repi, r35, 'o', ms=3, c=cols[0], label='R35')
plt.plot(repi, bj84, '-', c=cols[1], lw=2, label='BJ84')
plt.plot(repi, gs86, '-', c=cols[2], lw=2, label='GS86')
plt.plot(repi, gg91, '-', c=cols[3], lw=2, label='GG91')
plt.plot(repi, mlm92, '-', c=cols[4], lw=2, label='MLM92')
plt.plot(repi, w96, '-', c=cols[5], lw=2, label='W96')

plt.grid(which='both')
plt.xlabel('Epicentral Distance (km)')
plt.ylabel('-log A0')
plt.legend(loc=4)

plt.savefig('logA0_curves.png', fmt='png', bbox_inches='tight')
plt.show()
