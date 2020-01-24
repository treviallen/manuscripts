# -*- coding: utf-8 -*-
"""
Script to parse 'SRL_table.csv' and perform ordinary least squares and orthogonal
distance regression for SCR fault scaling coefficients.

Generates figure as used in:
    Clark, D., Lawrie, S., Brenn, G., Allen, T. I., Garthwaite, M. C., 
    and Standen, S. (2019). Surface deformation relating to the 2018 Lake 
    Muir earthquake sequence, south west Western Australia: insights into the 
    behaviour of Australian SCR earthquakes, Solid Earth.
    

Created on Wed May 01 14:09:02 2019

@author: u56903
"""
from numpy import array, arange, log10, argsort, copyto, zeros_like, where
from scipy import stats
from regressions import odr_lin_reg, lsq_lin_reg, lsq_reg_with_confidence
from misc_tools import get_log_xy_locs, get_mpl2_colourlist
from fault_tools import mag2ruplen_WC94, mag2srl_WC94, mag2srl_Cea14, mag2len_L14, \
                        mag2lsr_L14, len2mag_J94_SCR_quad, mag2len_J94_SCR_lin
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patheffects as PathEffects
mpl.style.use('classic')

cols = get_mpl2_colourlist()

# parse SRL table
eqn = []
year = []
mag = []
m0  = []
dsrl = []
vsrl = []
vdsp = []
minlen = []
maxlen = []
insar = []

lines = open('SRL_table.csv').readlines()

for line in lines[1:]:
    dat = line.strip().split(',')
    
    eqn.append(dat[0])
    year.append(int(dat[1]))
    mag.append(float(dat[2]))
    m0.append(float(dat[3]))
    dsrl.append(float(dat[4]))
    vsrl.append(float(dat[5]))
    vdsp.append(float(dat[6]))
    minlen.append(float(dat[7]))
    maxlen.append(float(dat[8]))
    insar.append(int(dat[-1]))
    
mag = array(mag)
dsrl = array(dsrl)
vsrl = array(vsrl)

yerr = []
for mil, mal in zip(minlen, maxlen):
    yerr.append(array([mil, mal]))
yerr=[dsrl-array(minlen), array(maxlen)-dsrl]
# sort
'''
idx = argsort(mag)
mag = mag[idx]
dsrl = dsrl[idx]
vsrl = vsrl[idx]
'''
##############################################################################
# plot DSRL data
##############################################################################
fig = plt.figure(1, figsize=(7,7))

# plot error first
plt.errorbar(mag, dsrl, yerr=yerr, marker='.', c='k', ls='none')
# loop through data and plot one-by-one
plt.semilogy(mag[0], dsrl[0], 's', c='k', ms=9, label='Original Data', zorder=100) #meckering
plt.semilogy(mag[1], dsrl[1], 's', c='k', ms=9, zorder=100) #calingiri
plt.semilogy(mag[2], dsrl[2], 's', c='k', ms=9, zorder=100) #cadoux
plt.semilogy(mag[3], dsrl[3], 's', c='k', ms=9, zorder=100) #MC
plt.semilogy(mag[4], dsrl[4], 's', c='k', ms=9, zorder=100) #TC
plt.semilogy(mag[5], dsrl[5], '^', c='k', ms=11, zorder=100) #katanning
plt.semilogy(mag[6], dsrl[6], 's', c='k', ms=9, zorder=100) #Ernabella
plt.semilogy(mag[7], dsrl[7], '^', c='0.75', ms=11, zorder=100) #PR
plt.semilogy(mag[8], dsrl[8], '^', c='0.75', ms=11, zorder=100) # LM
plt.semilogy(100, 10000, 's', c='0.75', ms=9, label='New Data', zorder=100) #random for legend

plt.grid(which='both')
plt.xlabel('Moment Magnitude', fontsize=14)
plt.ylabel('Detectable SRL (km)', fontsize=14)

# plt J94
loglrng = arange(0.1, 2, 0.02)
j94_mrng = len2mag_J94_SCR_quad(10**loglrng)
idx=where((j94_mrng >= 4.7) & (j94_mrng <= 6.8))[0]
#plt.semilogy(j94_mrng[idx], 10**loglrng[idx], '-', c=cols[0], lw=1.5, label='Johnston (1994; Quad)')

# plt WC94
mrng = arange(4.7, 6.81, 0.01)
j94_lrng = mag2len_J94_SCR_lin(mrng)
plt.semilogy(mrng, j94_lrng, '-', c=cols[0], lw=1.5, label='Johnston (1994; Lin)')

# plt WC94
mrng = arange(4.7, 6.81, 0.01)
wc94_lrng = mag2ruplen_WC94(mrng, 'rs')
plt.semilogy(mrng, wc94_lrng, '-', c=cols[1], lw=1.5, label='Wells & Coppersmith (1994; RLD)')

# plt Cea14
c14_lrng, sig = mag2srl_Cea14(mrng, 'au')
plt.semilogy(mrng, c14_lrng, '-', c=cols[2], lw=1.5, label='Clark et al. (2014)')

# plt L14
l14_lrng = mag2len_L14(mrng, 'scrrs')
plt.semilogy(mrng, l14_lrng[0], '-', c=cols[3], lw=1.5, label='Leonard (2014; FRL)')

leg1 = plt.legend(loc=4, numpoints=1, fontsize=11)

# make second legend
d1 = plt.semilogy(100, 10000, 's', c='w', ms=9, zorder=100) #random for legend
d2 = plt.semilogy(100, 10000, '^', c='w', ms=11, zorder=100) #random for legend
plt.legend([d1[0], d2[0]], ['Non-InSAR-based', 'InSAR-based'], loc=2, numpoints=1, fontsize=11)

plt.ylim([0.7, 60])
plt.xlim([4.5, 7])

plt.gca().add_artist(leg1)
plt.savefig('2020_scr_scaling_relations.png', fmt='png', dpi=300, bbox_inches='tight')
plt.show()











