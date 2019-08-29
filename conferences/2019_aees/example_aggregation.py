import matplotlib.pyplot as plt
from numpy import random, arange
from gmt_tools import cpt2colormap, remove_last_cmap_colour
from os import getcwd
import matplotlib as mpl
from matplotlib import colors, colorbar
mpl.style.use('classic')


##########################
# get colours
##########################

if getcwd().startswith('/nas'):
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/mi_pop.cpt'
else:
    cptfile = '//Users//trev//Documents//DATA//GMT//cpt//mi_pop.cpt'
ncols = 10
cmap, zvals = cpt2colormap(cptfile, ncols+1, rev=False)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncols)))

##########################
# set mmi
##########################

rand_mmi = 10*random.rand(25) / 2 + 1
rand_mmi[7] = 7.

idx = rand_mmi < 2.5
rand_mmi[idx] = 3.
rand_mmi[21] = 2.
rand_mmi[4] = 2.

locs = arange(0,51,10)

##########################
# set figs
##########################
fig = plt.figure(1, figsize=(12, 6))

ax = plt.subplot(121)

k = 0
for i in arange(0, 5):
    for j in arange(0, 5):

        # now plot
        pltx = [locs[i], locs[i+1], locs[i+1], locs[i], locs[i]]
        plty = [locs[j], locs[j], locs[j+1], locs[j+1], locs[j]]
        
        colidx = int(round(rand_mmi[k]))-1
        c= tuple(cs[colidx][:-1])
        plt.fill(pltx, plty, fc=c, ec='k', lw=0.75)
        
        k += 1

plt.xlabel('Easting (km)', fontsize=16)
plt.ylabel('Northing (km)', fontsize=16)
plt.text(1, 49, '(a)', va='top', ha='left', fontsize=20)

##########################
# set figs 2
##########################
ax = plt.subplot(122)

# now plot
i = 0
j = 0
pltx = [locs[i], locs[-1], locs[-1], locs[i], locs[i]]
plty = [locs[j], locs[j], locs[-1], locs[-1], locs[j]]

colidx = 7-1
c= tuple(cs[colidx][:-1])
plt.fill(pltx, plty, fc=c, ec='k', lw=0.75)

plt.xlabel('Easting (km)', fontsize=16)
plt.ylabel('Northing (km)', fontsize=16)
plt.text(1, 49, '(b)', va='top', ha='left', fontsize=20, zorder=100)

##########################
# make colorbar
##########################

# set colourbar
plt.gcf().subplots_adjust(bottom=0.22)
cax = fig.add_axes([0.25,0.05,0.5,0.06]) # setup colorbar axes.

norm = mpl.colors.Normalize(vmin=0.5, vmax=10.5)#myb
cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')

# set cb labels
ticks = range(1,11)
rom_num = ['I', 'II', 'III', 'IV', 'V', 'VI','VII','VIII','IX','X']
cb.set_ticks(ticks)
cb.set_ticklabels(rom_num)

titlestr = 'Macroseismic Intensity'
cb.set_label(titlestr, fontsize=18)



plt.savefig('example_aggregation.png', fmt='png', bbox_inches='tight')
plt.show()