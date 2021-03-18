from numpy import loadtxt, array, sqrt
import matplotlib.pyplot as plt
import matplotlib as mpl
from os import getcwd
mpl.style.use('classic')

fig = plt.figure(1, figsize=(7,7))

# parse NGA-E
csvfile = 'nga-e_mag-dist.csv'

ngaDat = loadtxt(csvfile, skiprows=1, delimiter=',')

plt.semilogx(ngaDat[:,0], ngaDat[:,1], '+', c='0.5', ms=7, label='NGA-East')

# parse WCA data
if getcwd().startswith('/Users'):
    csvfile = '../Western_Australia/aecom_record_list.csv'
else:
    csvfile = '20201210_aecom_record_list_no_dup.csv'
lines = open(csvfile).readlines()[1:]

cwa_mag = []
cwa_epi = []
cwa_dep = []
for line in lines:
    dat = line.strip().split(',')
    
    if getcwd().startswith('/Users'):
        cwa_epi.append(float(dat[10]))
        cwa_mag.append(float(dat[5]))
    else:
        cwa_epi.append(float(dat[11]))
        cwa_mag.append(float(dat[6]))
    cwa_dep.append(float(dat[3]))

cwa_mag = array(cwa_mag)
cwa_epi = array(cwa_epi)
cwa_dep = array(cwa_dep)

cwa_hyp = sqrt(cwa_epi**2 + cwa_dep**2)

plt.semilogx(cwa_hyp, cwa_mag, 'o', mec='r', mfc='none', ms=7, label='CWA')

plt.xlabel('Hypocentral Distance (km)', fontsize=20)
plt.ylabel('Magnitude', fontsize=20)
plt.legend(loc=2, numpoints=1)
plt.grid(which='both')
plt.xlim([1, 2000])
plt.ylim([2.5, 6.5]) 

plt.savefig('wca_ngae_mag_dist.png', dpi=600, fmt='png', bbox_inches='tight')
plt.show()