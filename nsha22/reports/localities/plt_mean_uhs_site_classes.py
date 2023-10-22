import matplotlib.pyplot as plt
from numpy import arange, array, interp, unique, vstack, log, exp
from os import path, getcwd, sep
from sys import argv
from gmt_tools import cpt2colormap 
from misc_tools import remove_last_cmap_colour, get_log_xy_locs, get_ga_master_colours_2022, get_ga_master_colours, get_mpl2_colourlist
#from mpl_toolkits.basemap import Basemap
#from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes
from hazard_tools import return_AS1170_4_shape

import matplotlib as mpl
mpl.style.use('classic')

'''
example run:
    run plt_uhs_curves.py results_fractilesUHS/hazard_uhs-mean_1.csv False

'''    
###############################################################################
# set params
###############################################################################

plt1170 = True


def split_config_lines(line):
    return '='.join(line.split('=')[1:]).strip().strip('\n').split(';')
    #return line.split('=')[1].split('#')[0].strip('\n').split(';')

uhsfile_b = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/complete_model/2023_final/results_uhs_fractiles_SC_B/hazard_uhs-mean_1.csv'
uhsfile_c = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/complete_model/2023_final/results_uhs_fractiles_SC_C/hazard_uhs-mean_1.csv'
uhsfile_d = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/complete_model/2023_final/results_uhs_fractiles_SC_D/hazard_uhs-mean_1.csv'
uhsfile_e = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/complete_model/2023_final/results_uhs_fractiles_SC_E/hazard_uhs-mean_1.csv'

uhsfile_b_84 = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/complete_model/2023_final/results_uhs_fractiles_SC_B/quantile_uhs-0.84_1.csv'
uhsfile_c_84 = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/complete_model/2023_final/results_uhs_fractiles_SC_C/quantile_uhs-0.84_1.csv'
uhsfile_d_84 = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/complete_model/2023_final/results_uhs_fractiles_SC_D/quantile_uhs-0.84_1.csv'
uhsfile_e_84 = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/source_models/complete_model/2023_final/results_uhs_fractiles_SC_E/quantile_uhs-0.84_1.csv'

'''
if pltProb == '10':
    probidx = 0
elif pltProb == '9':
    probidx = 1
elif pltProb == '2':
    probidx = 2
'''

pltLog = False

# get colours
if getcwd().startswith('/nas'):
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/Paired_12.cpt'
else:
    cptfile = '//Users//trev//Documents//DATA//GMT//cpt//Paired_12.cpt'

ncolours = 13
cmap, zvals = cpt2colormap(cptfile, ncolours)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncolours-1)))

#cs = get_mpl2_colourlist()

#get_mpl2_colourlist


###############################################################################
# parse uhs file
###############################################################################
def parse_uhs_file(uhsfile):
    lines = open(uhsfile).readlines()
    headers = [x for x in lines[1].strip().split(',')]
    
    # get keys from uhs file
    keys = lines[1].strip().split(',')[2:]
    
    # get peridos in keys
    periods = []
    tmpProb = []
    for key in keys:
        tmpProb.append(key.split('~')[0])
        
        if key.startswith('0.1'):
            if key.endswith('PGA'):
                periods.append(0.0)
            else:
                periods.append(float(key.split('(')[-1][:-1]))
    
    # get unique probabilities
    probabilities = unique(tmpProb)[::-1] # reorder
    
    
    # site site data
    uhsDict = []
    for line in lines[2:]:
        dat = [float(x) for x in line.strip().split(',')]
        tmpdict = {'lon':float(dat[0]), 'lat':float(dat[1])}
        
        for i, prob in enumerate(probabilities):
            startidx = i * len(periods) + 2
            stopidx = startidx + len(periods)
            siteUHS = [float(x) for x in dat[startidx:stopidx]]
            
            tmpdict[prob] = array(siteUHS)
    
        uhsDict.append(tmpdict)
        
    return probabilities, uhsDict, periods
    
probabilities, uhsDict_b, periods = parse_uhs_file(uhsfile_b)
probabilities, uhsDict_c, periods = parse_uhs_file(uhsfile_c)
probabilities, uhsDict_d, periods = parse_uhs_file(uhsfile_d)
probabilities, uhsDict_e, periods = parse_uhs_file(uhsfile_e)

probabilities, uhsDict_b_84, periods = parse_uhs_file(uhsfile_b_84)
probabilities, uhsDict_c_84, periods = parse_uhs_file(uhsfile_c_84)
probabilities, uhsDict_d_84, periods = parse_uhs_file(uhsfile_d_84)
probabilities, uhsDict_e_84, periods = parse_uhs_file(uhsfile_e_84)

probkey3pc = 1 # 3.3% in 50-yr
probkey3pc = 0

###################################################################################
# get mean UHS for all sites
###################################################################################
def get_mean_uhs(uhsDict, probkey):
    
    for i, uhs in enumerate(uhsDict):
        if i == 0:
            stackedUHS = log(uhs[probabilities[probkey]])
        else:
            stackedUHS = vstack((stackedUHS, log(uhs[probabilities[probkey]])))
            
    meanUHS = exp(stackedUHS.mean(axis=0))
    #plt_au_haz_map_ratios
    #meanUHS[0] = 2 * meanUHS[0] / 3 # applying 2/3 of 1/1500 AEP
    norm_meanUHS = meanUHS / meanUHS[0] 
    
    return meanUHS, norm_meanUHS, meanUHS[0]

meanUHS_b, norm_meanUHS_b, mean_pga_b = get_mean_uhs(uhsDict_b, probkey3pc)
meanUHS_c, norm_meanUHS_c, mean_pga_c = get_mean_uhs(uhsDict_c, probkey3pc)
meanUHS_d, norm_meanUHS_d, mean_pga_d = get_mean_uhs(uhsDict_d, probkey3pc)
meanUHS_e, norm_meanUHS_e, mean_pga_e = get_mean_uhs(uhsDict_e, probkey3pc)

meanUHS_b_84, norm_meanUHS_b_84, mean_pga_b_84 = get_mean_uhs(uhsDict_b_84, probkey3pc)
meanUHS_c_84, norm_meanUHS_c_84, mean_pga_c_84 = get_mean_uhs(uhsDict_c_84, probkey3pc)
meanUHS_d_84, norm_meanUHS_d_84, mean_pga_d_84 = get_mean_uhs(uhsDict_d_84, probkey3pc)
meanUHS_e_84, norm_meanUHS_e_84, mean_pga_e_84 = get_mean_uhs(uhsDict_e_84, probkey3pc)

# recalc norm 84th normalising by actual mean
norm_meanUHS_b_84 =  meanUHS_b_84 / mean_pga_b
norm_meanUHS_c_84 =  meanUHS_c_84 / mean_pga_c
norm_meanUHS_d_84 =  meanUHS_d_84 / mean_pga_d
norm_meanUHS_e_84 =  meanUHS_e_84 / mean_pga_e

###################################################################################
# set AS1170.4 shape for site class Be
###################################################################################

per1170 = arange(0, 4.1, 0.1)
        
shp1170_b = return_AS1170_4_shape(per1170, 'B')
shp1170_c = return_AS1170_4_shape(per1170, 'C')
shp1170_d = return_AS1170_4_shape(per1170, 'D')
shp1170_e = return_AS1170_4_shape(per1170, 'E')
           
###################################################################################
# plt mean UHS for different site classes
###################################################################################

fig = plt.figure(1, figsize=(11, 6))

plt.plot(periods, norm_meanUHS_b, lw=2.0, c=cs[1], label='NSHA23 SC B')
plt.plot(periods, norm_meanUHS_b_84, ls=':', lw=2.0, c=cs[1], label='NSHA23 SC B - 84th')
plt.plot(per1170, shp1170_b, ls='--', lw=2.0, c=cs[0], label='AS1170.4 SC B')
plt.plot(periods, norm_meanUHS_c, lw=2.0, c=cs[3], label='NSHA23 SC C')
plt.plot(periods, norm_meanUHS_c_84, ls=':', lw=2.0, c=cs[3], label='NSHA23 SC C - 84th')
plt.plot(per1170, shp1170_c, ls='--', lw=2.0, c=cs[2], label='AS1170.4 SC C')
plt.plot(periods, norm_meanUHS_d, lw=2.0, c=cs[5], label='NSHA23 SC D')
plt.plot(periods, norm_meanUHS_d_84, ls=':', lw=2.0, c=cs[5], label='NSHA23 SC D - 84th')
plt.plot(per1170, shp1170_d, ls='--', lw=2.0, c=cs[4], label='AS1170.4 SC D')
plt.plot(periods, norm_meanUHS_e, lw=2.0, c=cs[7], label='NSHA23 SC E')
plt.plot(periods, norm_meanUHS_e_84, ls=':', lw=2.0, c=cs[7], label='NSHA23 SC E - 84th')
plt.plot(per1170, shp1170_e, ls='--', lw=2.0, c=cs[6], label='AS1170.4 SC E')

plt.xlabel('Period (s)', fontsize=14)
plt.ylabel('Normalised '+str('%0.0f' % (100*float(probabilities[probkey3pc])))+'% in 50-year AEP Sa', fontsize=14)
plt.grid(which='both')
plt.xlim([0.0, 3.0])
plt.legend(loc=1, fontsize=11)

###############################################################################
# save figs
###############################################################################
pngname = 'norm_mean_uhs_site_classes.png'
plt.savefig(pngname, fmt='png', dpi=300, bbox_inches='tight')

plt.show()
