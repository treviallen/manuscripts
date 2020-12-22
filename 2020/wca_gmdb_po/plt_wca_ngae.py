from numpy import loadtxt, array, sqrt
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('classic')

fig = plt.figure(1, figsize=(7,7))

# parse NGA-E
csvfile = '/Users/trev/Documents/Earthquake_Data/NGA-East/nga-e_mag-dist.csv'

ngaDat = loadtxt(csvfile, skiprows=1, delimiter=',')

plt.semilogx(ngaDat[:,0], ngaDat[:,1], '+', c='0.5', ms=7, label='NGA-East')

# parse WCA data
csvfile = '../Western_Australia/aecom_record_list.csv'
lines = open(csvfile).readlines()[1:]

cwa_mag = []
cwa_epi = []
cwa_dep = []
for line in lines:
    dat = line.strip().split(',')
    cwa_mag.append(float(dat[5]))
    cwa_epi.append(float(dat[10]))
    cwa_dep.append(float(dat[3]))

cwa_mag = array(cwa_mag)
cwa_epi = array(cwa_epi)
cwa_dep = array(cwa_dep)

cwa_hyp = sqrt(cwa_epi**2 + cwa_dep**2)

plt.semilogx(cwa_hyp, cwa_mag, 'o', mec='r', mfc='none', ms=7, label='CWA')

plt.xlabel('Hypocentral Distance (km)')
plt.ylabel('Magnitude')
plt.legend(loc=2, numpoints=1)
plt.grid(which='both')
plt.xlim([1, 2000])
plt.ylim([2.5, 6.5]) 

plt.savefig('wca_ngae_mag_dist.png', fmt='png', bbox_inches='tight')
plt.show()