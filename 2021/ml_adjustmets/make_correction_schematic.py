from calculate_magnitudes import calc_R35, calc_MLM92, calc_HB87
import matplotlib.pyplot as plt
from numpy import array, arange, sqrt

import matplotlib as mpl
mpl.style.use('classic')
plt.rc('xtick',labelsize=13)
plt.rc('ytick',labelsize=13)

# 'Original $\mathregular{M_{LH}}$'
dep = 10.
repi = arange(0, 601, 1)
rhyp = sqrt(repi**2 + dep**2)

# get log A0 corrections
r35 = []
mlm92 = []
hb87 = []
for i in range(0, len(repi)):
   r35.append(calc_R35(1, 0.0, repi[i]))
   mlm92.append(calc_MLM92(1, 0.0, rhyp[i]))
   hb87.append(calc_HB87(1, 0.0, rhyp[i]))
mldiff = array(mlm92) - array(r35)
mldiff_modern = array(mlm92) - array(hb87)

# now plot
fig = plt.figure(1, figsize=(12, 4))

plt.plot([0, 600], [0, 0], 'k--')
plt.plot(repi, mldiff, '-', c='orangered', lw=2, label='Richter (1935)')
plt.plot(repi, mldiff_modern, '--', c='seagreen', lw=2, label='Hutton & Boore (1987)')

plt.xlabel('Epicentral Distance (km)', fontsize=16)
plt.ylabel('-log $\mathregular{A_0}$ (MLM92 - CA)', fontsize=16)

plt.xlim([0, 600])
plt.ylim([-0.8, 0.8])

# now plot stns
plt.plot(100, 0.2, '^', ms=18, c='dodgerblue')
txt = 'Site 1\nClosest Site = 100 km\nCorrection = '+str('%0.2f' % mldiff[100])+' mu'
plt.text(100, 0.3, txt, ha='center', va='bottom', fontsize=14)

plt.plot(510, 0.2, '^', ms=18, c='dodgerblue')
plt.plot([510, 510], [0, mldiff[510]], ':', c='dodgerblue', lw=2)
txt = 'Site 2\nClosest Site = 510 km\nCorrection = '+str('%0.2f' % mldiff[510])+' mu'
plt.text(510, 0.3, txt, ha='center', va='bottom', fontsize=14)

plt.legend(loc=3, fontsize=14)

plt.savefig('ml_correction_schematic.png', fmt='png', dpi=300, bbox_inches='tight')
plt.savefig('ml_correction_schematic.eps', fmt='eps', dpi=300, bbox_inches='tight')
plt.show()
