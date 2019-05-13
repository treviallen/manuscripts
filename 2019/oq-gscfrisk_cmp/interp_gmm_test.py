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
hdf5file = '/Users/trev/Documents/NRCan/2020_National_Hazard/2020_gsc_nshm/gmm_tables/gmm_hdf5_tables/AbrahamsonEtAl2015SInter.vs450.h15.hdf5'

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