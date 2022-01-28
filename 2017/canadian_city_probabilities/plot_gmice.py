#plot_gmice.py

from mmi_tools import mmi2pgm_wald99, mmi2pgm_worden12, mmi2pgm_atkinson07, \
                      mmi2pgm_dangkua11_ena, mmi2pgm_caprio15, mmi2pgv_newmark_rosenblueth, \
                      mmi2pga_gaull
                      
from numpy import arange, nan
import matplotlib.pyplot as plt

# set up figure
fig = plt.figure(1, figsize=(9, 8.4))
plt.tick_params(axis='both', labelsize=15)
plt.rcParams['pdf.fonttype'] = 42

# make mmi range
mmis = arange(2, 9.001, 0.001)

# get pga conversions
w99pga = mmi2pgm_wald99(mmis, 'pga')[0]
wea12pga = mmi2pgm_worden12(mmis, 'pga', nan, nan)[0]
ak07pga = mmi2pgm_atkinson07(mmis, 'pga')[0]
dc11pga = mmi2pgm_dangkua11_ena(mmis, 'pga', nan, nan)[0]
cea15pga = mmi2pgm_caprio15(mmis, 'pga')[0]
g79pga = mmi2pga_gaull(mmis)

ax = plt.subplot(2, 2, 1)
plt.semilogx(mmis*100000, mmis, '-', c='m', lw=2.) # dummmy plot
plt.semilogx(g79pga, mmis, '-', c='olive', lw=2.) # dummmy plot
plt.semilogx(w99pga, mmis, '-', color='c', lw=2.)
plt.semilogx(ak07pga, mmis, '-', color='limegreen', lw=2.)
plt.semilogx(dc11pga, mmis, 'b-', lw=2.)
plt.semilogx(wea12pga, mmis, '-', color='orange', lw=2.)
plt.semilogx(cea15pga, mmis, 'r-', lw=2.)
plt.xlabel('PGA $\mathregular{(cm/s^{2}}$)', fontsize=16)
plt.ylabel('MMI', fontsize=16)
plt.grid(which='both', color='0.5')
plt.legend(['NR71', 'G79', 'Wea99', 'AK07', 'DC11 (ENA)', 'Wea12', 'Cea15 (Global)'], loc=2, fontsize=10)
ticks = range(2,10)
rom_num = ['II', 'III', 'VI', 'V', 'VI','VII','VIII','IX']
plt.yticks(ticks, rom_num)
plt.ylim([4, 9])
plt.xlim([1, 1000])

# get pgv conversions
w99pgv = mmi2pgm_wald99(mmis, 'pgv')[0]
wea12pgv = mmi2pgm_worden12(mmis, 'pgv', nan, nan)[0]
ak07pgv = mmi2pgm_atkinson07(mmis, 'pgv')[0]
dc11pgv = mmi2pgm_dangkua11_ena(mmis, 'pgv', nan, nan)[0]
cea15pgv = mmi2pgm_caprio15(mmis, 'pgv')[0]
nr71pgv = mmi2pgv_newmark_rosenblueth(mmis)

ax = plt.subplot(2, 2, 2)
plt.semilogx(nr71pgv, mmis, '-', c='m', lw=2.)
plt.semilogx(w99pgv, mmis, '-', color='c', lw=2.)
plt.semilogx(ak07pgv, mmis, '-', color='limegreen',lw=2.)
plt.semilogx(dc11pgv, mmis, 'b-', lw=2.)
plt.semilogx(wea12pgv, mmis, '-', color='orange', lw=2.)
plt.semilogx(cea15pgv, mmis, 'r-', lw=2.)
plt.xlabel('PGV (cm/s)', fontsize=16)
#plt.ylabel('MMI', fontsize=16)
plt.grid(which='both', color='0.5')
plt.yticks(ticks, rom_num)
plt.ylim([4, 9])
plt.xlim([.1, 300])

# get Sa0.3 conversions
wea12pgv = mmi2pgm_worden12(mmis, 'sa03', nan, nan)[0]
ak07pgv = mmi2pgm_atkinson07(mmis, 'sa03')[0]
dc11pgv = mmi2pgm_dangkua11_ena(mmis, 'sa03', nan, nan)[0]
#cea15pgv = mmi2pgm_caprio15(mmis, 'pgv')[0]

ax = plt.subplot(2, 2, 3)
plt.semilogx(ak07pgv, mmis, '-', color='limegreen',lw=2.)
plt.semilogx(dc11pgv, mmis, 'b-', lw=2.)
plt.semilogx(wea12pgv, mmis, '-', color='orange', lw=2.)
#plt.semilogx(cea15pgv, mmis, 'r-', lw=2.)
plt.xlabel(r'SA 0.3 $\mathregular{(cm/s^{2}}$)', fontsize=16)
plt.ylabel('MMI', fontsize=16)
plt.grid(which='both', color='0.5')
plt.yticks(ticks, rom_num)
plt.ylim([4, 9])
plt.xlim([10, 10000])

# get Sa1.0 conversions
wea12pgv = mmi2pgm_worden12(mmis, 'sa10', nan, nan)[0]
ak07pgv = mmi2pgm_atkinson07(mmis, 'sa10')[0]
dc11pgv = mmi2pgm_dangkua11_ena(mmis, 'sa10', nan, nan)[0]
#cea15pgv = mmi2pgm_caprio15(mmis, 'pgv')[0]

ax = plt.subplot(2, 2, 4)
plt.semilogx(ak07pgv, mmis, '-', color='limegreen', lw=2.)
plt.semilogx(dc11pgv, mmis, 'b-', lw=2.)
plt.semilogx(wea12pgv, mmis, '-', color='orange', lw=2.)
#plt.semilogx(cea15pgv, mmis, 'r-', lw=2.)
plt.xlabel(r'SA 1.0 $\mathregular{(cm/s^{2}}$)', fontsize=16)
#plt.ylabel('MMI', fontsize=16)
plt.grid(which='both', color='0.5')
plt.yticks(ticks, rom_num)
plt.ylim([4, 9])
plt.xlim([.1, 1000])

plt.savefig('gmice_comparison.png', format='png', dpi=150, bbox_inches='tight')
plt.show()
