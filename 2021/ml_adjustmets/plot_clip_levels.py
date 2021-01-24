from calc_oq_gmpes import allen2012_gsim
from numpy import array, arange, log, exp, pi, interp, where
from scipy.constants import g
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('classic')


# parse simulated pgv
pgv_csv = 'smsim/simulated_pgv.csv' # output from running "smsim/simulate_acceleration_data.py"

lines = open(pgv_csv).readlines()

mags = []
dists = []
pgvs = []

for line in lines[1:]:
    dat = line.strip().split(',')
    mags.append(float(dat[0]))
    dists.append(float(dat[1]))
    pgvs.append(float(dat[2]))

mags = array(mags)
dists = array(dists)
pgvs = array(pgvs)

# plot pgv - dist
fig = plt.figure(1, figsize=(9,5))
idx = where(mags == 4.0)[0]
plt.loglog(dists[idx], pgvs[idx], '^', mec='r', mfc='none', label='ML 4.0')

idx = where(mags == 4.5)[0]
plt.loglog(dists[idx], pgvs[idx], 's', mec='g', mfc='none', label='ML 4.5')

idx = where(mags == 5.0)[0]
plt.loglog(dists[idx], pgvs[idx], 'o', mec='b', mfc='none', label='ML 5.0')

# plt approx wwssn SP level
plt.plot([10,500],[0.001, 0.001], 'k--', lw=2, label='WWSSN')
plt.plot([10,500],[1.3, 1.3], 'r--', lw=2, label='STS-2')

plt.xlabel('Hypocentral Distance (km)')
plt.ylabel('PGV (cm/s)')
plt.grid(which='both')
plt.xlim([10,500])
plt.legend(loc=3, fontsize=10)

plt.show()




'''
rrups = arange(1,301,10)
mags = arange(2,6.1,0.1)
dep = 10.
wwssn_clip_vel = 10**(-100. / 20.) # m/s (or -100 dB)
wwssn_clip_vel = 10**(-50. / 20.) # m/s (or -100 dB)

f = 2 # Hz
wwssn_clip_acc = wwssn_clip_vel * 2 * pi * f / g # in g
ln_wwssn_clip_acc = log(wwssn_clip_acc) # in ln g

rrup_clip = []
for rrup in rrups:
    sa05 = []
    
    for mag in mags:
        a12 = allen2012_gsim(mag, dep, rrup)
        
        sa05.append(a12['sa'][11]) # in ln g
        
    rrup_clip.append(interp(ln_wwssn_clip_acc, sa05, mags))
'''