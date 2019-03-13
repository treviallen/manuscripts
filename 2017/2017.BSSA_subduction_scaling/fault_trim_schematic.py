# -*- coding: utf-8 -*-
"""
Created on Fri Oct 09 11:13:40 2015

@author: tallen
"""

# http://stackoverflow.com/questions/23477775/matplotlib-mixed-subplots-with-plot-and-pcolormesh-to-have-the-same-x-axis-len

import matplotlib.pyplot as plt
from numpy import array, arange, random, linspace, meshgrid, dstack, amax, \
                  percentile, shape, where, mean
from scipy.stats import multivariate_normal
from gmt_tools import cpt2colormap
import matplotlib 
matplotlib.rc('xtick', labelsize=16) 
matplotlib.rc('ytick', labelsize=16) 

cptfile = 'U:\\DATA\\GMT\\cpt\\grays.cpt'
cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//grays.cpt'

ncolours = 11
cmap, zvals = cpt2colormap(cptfile, ncolours, rev=False)
#cmap = plt.cm.get_cmap('gray_r', ncolours)
flen = 400 # km
fwid = 200 # km
subfsize = 10 # km
#cmap = cm.rainbow
vmin = 0.
vmax = 33.

random.seed(1000)
# Generate some interesting-looking random data...
xnum = flen/subfsize + 1
ynum = fwid/subfsize + 1
rndgrid = random.rand(ynum, xnum)
x = linspace(0, flen, xnum)
y = linspace(0, fwid, ynum)

xv, yv = meshgrid(x, y)

xx = linspace(0, 1, xnum)
yy = linspace(0, 1, ynum)
xm, ym = meshgrid(xx, yy)
pos = dstack((xm, ym))

##############################################################################
# make random asperities
##############################################################################

# asperity 1
rv = multivariate_normal([0.5, .55], [[.05, 0.05], [.02, 0.04]])
asp1 = rv.pdf(pos)

# asperity 2
rv = multivariate_normal([0.7, .5], [[.02, 0.02], [.02, 0.05]])
asp2 = rv.pdf(pos)

asperities = 21*(asp1/amax(asp1) + asp2/amax(asp2)) - rndgrid
idx = asperities < 0.
asperities[idx] = 0.

'''
# test for cells LT threshold
idx = asperities / amax(asperities) < 0.1
asperities[idx] = amax(asperities)
'''
'''
# test to see what gets dropped
for i in range(0, shape(asperities)[1]):
    idx = xx <= 0.
    asperities[:, idx] = amax(asperities)/2.
    idx = xx >= 1.
    asperities[:, idx] = amax(asperities)
'''    
# !!! pcolormesh drops last row and column !!!    

fig = plt.figure(1, figsize=(16, 6))
im = plt.pcolormesh(xv, yv, asperities, edgecolors='0.3', lw=0.1, cmap=cmap)
plt.clim([vmin,vmax])

##############################################################################
# do trimming
##############################################################################

# set vectors
widvect = []
lenvect = []
slipthresh = 0.15
tcol = '0.4'

# get downdip width
for i in range(0, shape(asperities)[1]-1):
    normslip = asperities[:-1,i] / amax(asperities) # drop last row
    idx = where(normslip > slipthresh)[0]
    
    if len(idx) > 0: 
        # get max/min indices
        maxidx = max(idx)
        minidx = min(idx)    
        # plot bounds
        plt.plot([x[i], x[i+1]], [y[maxidx+1], y[maxidx+1]], tcol, lw=4)
        plt.plot([x[i], x[i+1]], [y[minidx], y[minidx]], tcol, lw=4)
        widvect.append(subfsize*(maxidx-minidx+1))
    else:
        widvect.append(0)

widpctle = percentile(array(widvect), 75.)

# get along-stk length
for j in range(0, shape(asperities)[0]-1):
    normslip = asperities[j,:-1] / amax(asperities) # drop last row
    idx = where(normslip > slipthresh)[0]
    
    if len(idx) > 0: 
        # get max/min indices
        maxidx = max(idx)
        minidx = min(idx)    
        # plot bounds
        plt.plot([x[maxidx+1], x[maxidx+1]], [y[j], y[j+1]], tcol, lw=4)
        plt.plot([x[minidx], x[minidx]], [y[j], y[j+1]], tcol, lw=4)
        lenvect.append(subfsize*(maxidx-minidx+1))
    else:
        lenvect.append(0)

lenpctle = percentile(array(lenvect), 75.)

print lenpctle, widpctle

##############################################################################
# plot trimmed fault
##############################################################################

# get centroid
normslip = asperities[:-1,:-1] / amax(asperities)
idx = normslip >= 0.25
centx = mean(xv[:-1,:-1][idx]) + subfsize/2
centy = mean(yv[:-1,:-1][idx]) + subfsize/2

#plt.plot(centx, centy, 'go', ms=10)

# get fault bounds
lx = centx - lenpctle/2.
rx = centx + lenpctle/2.
ty = centy - widpctle/2.
by = centy + widpctle/2.
fx = [lx, rx, rx, lx, lx]
fy = [ty, ty, by, by, ty]

# now plt
plt.plot(fx, fy, 'k-', lw=4)


##############################################################################
# plt colorbar
##############################################################################

#plt.axis('equal')
plt.xlim([0, flen])
plt.ylim([fwid, 0])
plt.ylabel('FFRM Width (km)', fontsize=20)
plt.xlabel('FFRM Length (km)', fontsize=20)

cb = fig.colorbar(im, ticks=arange(0,vmax+1,3))
cb.set_label('Slip (m)', fontsize=20, rotation=270, labelpad=30)
#

# save fig
plt.savefig('fault_trim_scematic.pdf', dpi=150, format='pdf', bbox_inches='tight')

plt.show()
