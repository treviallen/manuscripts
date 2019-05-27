# -*- coding: utf-8 -*-
"""
Created on Wed May 01 14:09:02 2019

@author: u56903
"""
from numpy import array, arange, log10, argsort
from scipy import stats
from regressions import odr_lin_reg, lsq_lin_reg, lsq_reg_with_confidence
from fault_tools import mag2ruplen_WC94, mag2srl_WC94, mag2srl_Cea14, mag2len_L14
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('classic')

# parse SRL table
eqn = []
year = []
mag = []
m0  = []
dsrl = []
vsrl = []
vdsp = []

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
    
mag = array(mag)
dsrl = array(dsrl)
vsrl = array(vsrl)

# sort
idx = argsort(mag)
mag = mag[idx]
dsrl = dsrl[idx]
vsrl = vsrl[idx]

##############################################################################
# plot VSRL data
##############################################################################

fig = plt.figure(1, figsize=(13,6))
ax = plt.subplot(121)

plt.semilogy(mag, vsrl, 's', c='seagreen', label='VSRL Data')
plt.grid(which='both')
plt.xlabel('Moment Magnitude', fontsize=14)
plt.ylabel('Visible SRL (km)', fontsize=14)

# plt WC94
mrng = arange(4.7, 6.81, 0.001)
wc94_lrng, sig = mag2srl_WC94(mrng, 'rs')
plt.semilogy(mrng, wc94_lrng, '-', c='dodgerblue', lw=1.5, label='Wells & Coppersmith (1994; SRL)')

# plt Cea14
c14_lrng, sig = mag2srl_Cea14(mrng, 'au')
plt.semilogy(mrng, c14_lrng, '-', c='orange', lw=1.5, label='Clark et al. (2014)')

# plt L14
l14_lrng = mag2len_L14(mrng, 'scrrs')
plt.semilogy(mrng, l14_lrng, '-', c='r', lw=1.5, label='Leonard (2014)')

# do regression
start = [1., -5.]
m, c = lsq_lin_reg(mag, log10(vsrl), start)
lower, upper = lsq_reg_with_confidence(mag, log10(vsrl))
plt.fill_between(mag, 10**lower, 10**upper, color='0.85', edgecolor='', zorder=0)

# plot new relationship
lrng = 10**(m * mrng + c)
plt.semilogy(mrng, lrng, '-', c='k', lw=1.5, label='Present Study')

plt.legend(loc=4, numpoints=1, fontsize=11)
plt.ylim([0.1, 100])

##############################################################################
# plot DSRL data
##############################################################################

ax = plt.subplot(122)

plt.semilogy(mag, dsrl, 's', c='seagreen', label='DSRL Data')
plt.grid(which='both')
plt.xlabel('Moment Magnitude', fontsize=14)
plt.ylabel('Detectable SRL (km)', fontsize=14)

# plt WC94
mrng = arange(4.7, 6.81, 0.01)
wc94_lrng = mag2ruplen_WC94(mrng, 'rs')
plt.semilogy(mrng, wc94_lrng, '-', c='dodgerblue', lw=1.5, label='Wells & Coppersmith (1994; RLD)')

# plt Cea14
c14_lrng, sig = mag2srl_Cea14(mrng, 'au')
plt.semilogy(mrng, c14_lrng, '-', c='orange', lw=1.5, label='Clark et al. (2014)')

# plt L14
l14_lrng = mag2len_L14(mrng, 'scrrs')
plt.semilogy(mrng, l14_lrng, '-', c='red', lw=1.5, label='Leonard (2014)')

# do regression
start = [1., -5.]
m, c = lsq_lin_reg(mag, log10(dsrl), start)
lower, upper = lsq_reg_with_confidence(mag, log10(dsrl))
plt.fill_between(mag, 10**lower, 10**upper, color='0.85', edgecolor='', zorder=0)

# plot new relationship
lrng = 10**(m * mrng + c)
plt.semilogy(mrng, lrng, '-', c='k', lw=1.5, label='Present Study')

plt.ylim([0.1, 100])
plt.legend(loc=4, numpoints=1, fontsize=11)

##############################################################################
# export & show
##############################################################################

plt.savefig('2019_srl_regression.png', fmt='png', dpi=300, bbox_inches='tight')
plt.show()

##############################################################################
# get ratio of DSRL vs VSRL
##############################################################################


fig = plt.figure(1, figsize=(13,6))
ax = plt.subplot(121)

rat_vdsrl = vsrl/dsrl
idx = rat_vdsrl != 1.
plt.semilogy(mag[idx], rat_vdsrl[idx], 's', c='seagreen', label='DSRL Data')

plt.grid(which='both')
plt.xlabel('Moment Magnitude', fontsize=14)
plt.ylabel('VSRL/DSRL', fontsize=14)

plt.savefig('2019_vsrl-dsrl_ratio.png', fmt='png', dpi=300, bbox_inches='tight')
plt.show()
