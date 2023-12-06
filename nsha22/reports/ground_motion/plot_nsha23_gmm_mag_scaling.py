import warnings
warnings.filterwarnings('ignore')


'''
start main
'''
# start of plotting routine
from numpy import array, arange, sqrt, exp, log, log10, logspace, argwhere, interp, unique, vstack, nan
from os import path, walk, system
from sys import argv
import matplotlib.pyplot as plt
from gmt_tools import cpt2colormap
from calc_oq_gmpes import nsha23_gsims
from misc_tools import get_log_xy_locs, remove_last_cmap_colour, get_mpl2_colourlist
from mapping_tools import distance
import matplotlib as mpl
mpl.style.use('classic')

'''
usage: run plot_scr_attenuation.py Moe_5.4/psa True
'''

plt.rcParams['pdf.fonttype'] = 42
import matplotlib 
matplotlib.rc('xtick', labelsize=13) 
matplotlib.rc('ytick', labelsize=13) 

prefix = argv[1]
colTrue = 'True'

# set event details
if prefix.startswith('cratonic'):
    
    dep = 10.0


ztor = 8. # guess
rake = 0. # USGS CMT
dip  = 90.

letters = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)', '(g)', '(h)', '(i)']

vs30 = 760.

rjb = 20 # km - hardwired
rrup = sqrt(rjb**2 + dep**2) # assume point source; i.e. repi = rjb
mags = arange(4.5, 7.6, 0.5)

fig = plt.figure(1, figsize=(11, 13))

ncols = 7
cptfile = '/Users/trev/Documents/DATA/GMT/cpt/gay-flag-1978.cpt'
cptfile = '/Users/trev/Documents/DATA/GMT/cpt/keshet.cpt'
cmap, zvals = cpt2colormap(cptfile, ncols+1)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncols)))

cs = get_mpl2_colourlist()

titles = ['PGA', 'Sa(0.2)','Sa(0.5)','Sa(2.0)']
Tplot = [0.0, 0.2, 0.5, 2.0]
props = dict(boxstyle='round', facecolor='w', alpha=1)

'''
titles = ['Sa(0.01)','Sa(0.2)','Sa(1.0)','Sa(2.0)']
Tplot = [0.01, 0.2, 1.0, 2.0]

titles = ['PGA','Sa(0.2)','Sa(1.0)','Sa(2.0)']
Tplot = [0.0, 0.2, 1.0, 2.0]
'''
#Tplot = [0.0]
# loop thru periods
for j, t in enumerate(Tplot):
    ax = plt.subplot(3, 2, j+1)
    AB06r = []
    Sea09_NCr = []
    Sea09_YCr = []
    A12r = []
    D15r = []
    ESHM20r = []
    NGAEr = []
    for i,mag in enumerate(mags):
        
        r = rrup
        # get ground motion estimates from GMPEs
        vs30 = 760.
        AB06imt, Sea09imt, Sea09YCimt, A12imt, D15imt, NGAEimt, ESHM20Cimt \
                 = nsha23_gsims(mag, dep, ztor, dip, rake, rrup, rjb, vs30)
        
        syms = ['s', 'x', '+', '1', '2', 'o']
        
                    
        if t == 0.0:
            AB06r.append(AB06imt['pga'][0])
            Sea09_NCr.append(Sea09imt['pga'][0])
            Sea09_YCr.append(Sea09YCimt['pga'][0])
            A12r.append(A12imt['sa'][0])
            #AA13r.append(AA13imt['pga'][0])
            D15r.append(D15imt['pga'][0])
            NGAEr.append(NGAEimt['sa'][0])
            ESHM20r.append(ESHM20Cimt['pga'][0])
            
        else:
            # interpolate log values to correct period
            AB06r.append(interp(t, AB06imt['per'], AB06imt['sa']))
            Sea09_NCr.append(interp(t, Sea09imt['per'], Sea09imt['sa']))
            Sea09_YCr.append(interp(t, Sea09YCimt['per'], Sea09YCimt['sa']))
            A12r.append(interp(t, A12imt['per'], A12imt['sa']))
            D15r.append(interp(t, D15imt['per'], D15imt['sa']))
            NGAEr.append(interp(t, NGAEimt['per'], NGAEimt['sa']))
            ESHM20r.append(interp(t, ESHM20Cimt['per'], ESHM20Cimt['sa']))
            
    if colTrue == 'True':
        h1 = plt.loglog(rjb, exp(AB06r),  ls='-', lw=1., color=cs[0])
        h2 = plt.loglog(rjb, exp(Sea09_NCr), ls='-', lw=1., color=cs[1])
        h3 = plt.loglog(rjb, exp(Sea09_YCr), ls='-', lw=1., color=cs[2])
        h4 = plt.loglog(rjb, exp(A12r),   ls='-', lw=1., color=cs[3])
        h5 = plt.loglog(rjb, exp(D15r), ls='-', lw=1., color=cs[4])
        h6 = plt.loglog(rjb, exp(NGAEr),  ls='-', lw=1., color=cs[5])
        h7 = plt.loglog(rjb, exp(ESHM20r), ls='-', lw=1., color=cs[6])
        
        
    if j >= 2:
        plt.xlabel('$\mathregular{R_{JB}}$ (km)', fontsize=16)
    
    if j == 0 or j == 2 or j == 4:
        plt.ylabel('Spectral Acceleration (g)', fontsize=16)
        
    plt.xlim([8, 600])
    plt.ylim([5E-4, 2])
       
    #plt.title(titles[j])
    xtxt = get_log_xy_locs(ax.get_xlim(), 0.95)
    ytxt = get_log_xy_locs(ax.get_ylim(), 0.95)
    plt.text(xtxt, ytxt, titles[j], size=17, horizontalalignment='right', verticalalignment='top', weight='normal', bbox=props)
    plt.grid(which='both', color='0.5')
    ylims = ax.get_ylim()
    plt.text(7., ylims[1]*1.25, letters[j], va='bottom', ha ='right', fontsize=16)

    if j == 0:
        '''
        plt.legend((h1[0], h2[0], h3[0], h4[0], h5[0], h6[0], h7[0]) \
                   ['AB06','Sea09(NC)','Sea09(YC)*', 'A12', 'D15', 'NGA-E', 'ESHM20*'],loc=3,numpoints=1,fontsize=11)
        '''           
        plt.legend((h1[0], h2[0], h3[0], h4[0], h5[0], h6[0], h7[0]), \
                   ['AB06','Sea09(NC)','Sea09(YC)', 'A12', 'D15', 'NGA-E', 'ESHM20'],loc=3,numpoints=1,fontsize=11)
                   
#                   'AB06','Sea09(NC)','Sea09(YC)', 'A12', 'D15', 'NGA-E', 'ESHM20'

'''
# export table for russ
f = open('D15_gmm_table.csv', 'w')
f.write(tabtxt)
f.close()
'''      
plt.savefig(prefix+'_mag_scaling.png', format='png', dpi=300, bbox_inches='tight')
plt.show()