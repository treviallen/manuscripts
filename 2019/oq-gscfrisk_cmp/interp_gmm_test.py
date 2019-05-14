from numpy import exp
from misc_tools import get_log_xy_locs
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('classic')

# set consts
mag = 7.15
dep = 10.
ztor = 7.
dip = 30.
rake = 90.
hdf5file = '/Users/trev/Documents/NRCan/2015_National_Hazard/2015_gsc_nshm/gm_tables/Wcrust_med_clC.hdf5'
#hdf5file = '/Users/trev/Documents/NRCan/2020_National_Hazard/2020_gsc_nshm/gmm_tables/gmm_hdf5_tables/AbrahamsonEtAl2015SInter.vs450.h15.hdf5'

# parse frisk GMM
rhyp = []
gsc_sa02 = [] # in g
gsc_sa10 = []

lines = open('gscfrisk_gmm_interp.csv').readlines()
for line in lines[1:]:
	dat = line.split(',')
	rhyp.append(float(dat[1]))
	gsc_sa02.append(float(dat[-1]))
	gsc_sa10.append(float(dat[5]))
	

# now get gmm interpolator
from calc_oq_gmpes import hdf5_gsim

oq_sa02 = [] # in g
oq_sa10 = []

# loop thru dists
for r in rhyp:

	imt = hdf5_gsim(mag, dep, ztor, dip, rake, r, r, r, 450., hdf5file)
	oq_sa02.append(exp(imt['sa'][2]))
	oq_sa10.append(exp(imt['sa'][5]))
	
# now plot
fig = plt.figure(1, figsize=(12, 6))

ax = plt.subplot(121)
plt.loglog(rhyp, gsc_sa02, '-', c='dodgerblue', lw=2.0, label='GSCFRISK')
plt.loglog(rhyp, oq_sa02, '--', c='darkred', lw=3.0, label='OpenQuake-engine')
plt.xlim([8, 500])
plt.ylim([0.001, 10])
plt.xlabel('Hypocentral Distance (km)', fontsize=16)
plt.ylabel('Spectral Aceleration (g)', fontsize=16)
plt.title ('Sa[0.2 s]', fontsize=20)
plt.grid(which='both')
plt.legend(loc=1, fontsize=14)

ylims = ax.get_ylim()
ypos = get_log_xy_locs(ylims, 0.98)
plt.text(9, ypos, '(a)', fontsize=20, ha='left', va='top')


ax = plt.subplot(122)
plt.loglog(rhyp, gsc_sa10, '-', c='dodgerblue', lw=2.0, label='GSCFRISK')
plt.loglog(rhyp, oq_sa10, '--', c='darkred', lw=3.0, label='OpenQuake-engine')
plt.xlim([8, 500])
#plt.ylim([0.001, 10])
plt.xlabel('Hypocentral Distance (km)', fontsize=16)
#plt.ylabel('Spectral Aceleration (g)', fontsize=16)
plt.title ('Sa[1.0 s]', fontsize=20)
plt.grid(which='both')

ylims = ax.get_ylim()
ypos = get_log_xy_locs(ylims, 0.98)
plt.text(9, ypos, '(b)', fontsize=20, ha='left', va='top')

plt.savefig('oq-interp_test.png', fmt='png', bbox_inches='tight')
plt.show()
