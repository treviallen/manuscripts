# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 14:50:20 2015

@author: tallen

Code generates a python dictionary of amplification factors
relative to B/C site conditions

To access amp factors for a given PGAref, site class, type:

	> C_ampfactors[PGAref][Vs30]

	e.g. 
	
	IN > C_ampfactors['0.1g']['vs1100']
	
	OUT > array([ 0.84782554,  0.75796444,  0.69797956,  0.64357284,  0.62994067,
                0.62824066,  0.6411572 ,  0.69124828,  0.85792907,  0.67049578]) 
	
Outputs CSV file of amp factors relative to B/C site conditions

"""

#from boore_atkinson_site_2008 import boore_atkinson_siteamp
#import matplotlib as mpl
from atkinson_boore_site_2006 import atkinson_boore_siteamp
from seyhan_stewart_2014 import seyhan_stewart_siteamp
from hazard_tools import return_AS1170_4_shape
from numpy import array, arange, where, log, exp, interp, hstack
from gmt_tools import cpt2colormap
from misc_tools import remove_last_cmap_colour, get_log_xy_locs
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('classic')

# periods for NBCC2015 (PGA = 0; PGV = -1)     
T = array([0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.25, \
           1.5, 2., 2.25, 2.5, 2.75, 3., 3.25, 3.5, 4., 4.5, 5., 0., -1.]) # in s
asT = arange(0.05, 5., 0.01)
#T = [0.1, 0.2, 0.3, 0.5, 1., 2., 5., 10., 0., -1.] # in s - NBCC


# Set Vs30 array - equivalent to 
vs30 = [150, 270, 464, 760, 1100] # in m/s 
rev_vs30 = vs30[::-1] # reverse vs30 order for tables
as1170_site_class = ['E', 'D', 'C', 'B', 'A']

# frist get B/C pga4nl
target_pga4nl_BC = array([0.1, 0.2, 0.3, 0.4, 0.5]) # in g

# get B/C linear ground-motions that provide target_pga4nl_C values - AW rounding here led to minor differences in final factors
refVs30 = 760.  # reset to B/C
refT    = 0.0 # i.e. PGA
refPGA  = 0.1 # in g
pga4nl_BC = target_pga4nl_BC / seyhan_stewart_siteamp(refVs30, refT, refPGA) # B/C to C conversion factor for PGA of 1.208 based on BA08

#pga4nl_BC = target_pga4nl_C / 1.  # for testing purposes - needed to match Gail's BA08 C->B/C factors

##########################################################################
'''
Get BA08 amp factors for site class BC - this loop is not required for 
NBCC2015 site class tables, but kept for testing purposes
'''

# set empty dictionary
BC_ampfactors = {}

# loop thru pga4nl
for i, pga in enumerate(pga4nl_BC):
    # set dictionary and loop thru velocity
    vdict = {}
    # loop thru Vs30
    for v, asc in zip(vs30, as1170_site_class):
        # loop thru periods
        amp = []
        for t in T:
            if v < 1000:
                # get amp factor for each period and append to amp
                amp.append(seyhan_stewart_siteamp(v, t, pga))
                
            else:
                amp.append(atkinson_boore_siteamp(v, t, pga))
        
        # add amp factors to velocity dictionary
        vdict['vs'+str(v)] = array(amp)
    
    # add velocity dictionary to pga4nl dictionary
    BC_ampfactors[str(target_pga4nl_BC[i])+'g'] = vdict
    
BC_ampfactors['periods'] = array(T)

##########################################################################
'''
Re-cast amp factors in terms of site class B/C 
'''
# set empty dictionary
C_ampfactors  = {}

# loop thru pga4nl
for i, pga in enumerate(pga4nl_BC):
    # set dictionary and loop thru velocity
    vdict = {}
    asvdict = {}
    for v, asc in zip(vs30, as1170_site_class):
        # loop thru periods
        amp = []
        for t in T:
            #if v < 1000:
            # get amp factor for each period relative to C and append to amp
            amp.append(seyhan_stewart_siteamp(v, t, pga) / 
                       seyhan_stewart_siteamp(760., t, pga))
                       
            '''
            else:
                amp.append(atkinson_boore_siteamp(v, t, pga))
            '''
        
        # add amp factors to velocity dictionary
        vdict['vs'+str(v)] = array(amp)
        
        # add AS1170 factors
        asv = return_AS1170_4_shape(asT, asc) / return_AS1170_4_shape(asT, 'B')
        asvdict['as_vs'+str(v)] = asv
        
    
    # add velocity dictionary to amp factor dictionary
    C_ampfactors[str(target_pga4nl_BC[i])+'g'] = vdict
    C_ampfactors['as'+str(target_pga4nl_BC[i])+'g'] = asvdict

# add periods to amp factor dictionary 
C_ampfactors['periods'] = array(T)

##########################################################################

##########################################################################
'''
plt amp factors
'''

#import matplotlib as mpl

pltpga = [0.1, 0.4]

mpl.rcParams['pdf.fonttype'] = 42
plt.rcParams['xtick.labelsize'] = 14
plt.rcParams['ytick.labelsize'] = 14 

numcols = 6
cptfile = '/Users/trev/Documents/DATA/GMT/cpt/temperature.cpt'
cptfile = '//Users//trev//Documents//DATA//GMT//cpt//qual-dark-06.cpt' # try me!
cptfile = '//Users//trev//Documents//DATA//GMT//cpt//humidity.cpt'
cmap, zvals = cpt2colormap(cptfile, numcols)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(numcols)))

# remove yellow
#idx=array([0, 1, 3, 4, 5])
#cs = cs[idx]

legtxt = ['Site Class B (Ae)', 'Site Class B/C (Be)', 'Site Class C (Ce)', 'Site Class D (De)', \
          'Site Class E (Ee)'] #, 'F(T)', 'Fa & Fv']

figure = plt.figure(1,figsize=(9,13))

# Boarcherdt 1994 for 0.1 and 0.4 g
Fa01 = [0.6, 0.7, 1.0, 1.5, 1.5]
Fa04 = [1.1, 1.0, 1.0, 0.9, 0.9]
Fv01 = [0.4, 0.6, 1.0, 2.0, 2.0]
Fv04 = [0.6, 0.7, 1.0, 1.6, 1.6]
FaT  = [0.1, 0.5]
FvT  = [0.4, 2.0]


# Finn & Whiteman for 0.1 and 0.4 g - Relative to site Sa(0.2 s)
Fa01 = array([0.7, 0.8, 1.0, 1.3, 2.1])  
Fa04 = array([1.1, 1.0, 1.0, 0.9, 0.9])  
Fv01 = array([0.4, 0.6, 1.0, 2.0, 2.0])  
Fv04 = array([0.6, 0.7, 1.0, 1.6, 1.6])  
FaT = [0.1, 0.5]
FvT = [0.4, 2.0] 


for i, pga in enumerate(pltpga):
    ax = plt.subplot(2, 1, i+1)
    plt.grid(which='major', color='gray', linestyle='--', linewidth=0.5)
    plt.grid(which='minor', color='gray', linestyle='--', linewidth=0.5)
    # loop thru vs
    for j, v in enumerate(rev_vs30):
        pltamps = C_ampfactors[str(pga)+'g']['vs'+str(v)][:-2]
        plt.semilogx(T[:-2], pltamps, '-', lw=2.5,color=cs[j], label=legtxt[j]) 
        	
        aspltamps = C_ampfactors['as'+str(pga)+'g']['as_vs'+str(v)]
        plt.semilogx(asT, aspltamps, '--', lw=2.5,color=cs[j])       
        	
    plt.ylim([0.5, 6.])
    plt.title('PGAref = ' + str(pga) + ' g', fontsize=20)
    plt.xlabel('Period (s)', fontsize=16)
    plt.ylabel('Amplification Factor (Relative to B/C)', fontsize=16)
    
    '''
    if i == 0:
        plt.legend(legtxt ,loc=2, fontsize=14)
    else:
        h1 = plt.semilogx([0.1, 1.],[9999., 9999.], 'k-', lw=3.)
        h2 = plt.semilogx([0.1, 1.],[9999., 9999.], 'k--', lw=3.)
        plt.legend((h1[0], h2[0]),('F(T)', 'Fa & Fv') ,loc=2, fontsize=14)
    '''
    # at John's request
    h1 = plt.semilogx([0.1, 1.],[9999., 9999.], 'k-', lw=2.5)
    h2 = plt.semilogx([0.1, 1.],[9999., 9999.], 'k--', lw=2.5)
    
    if i == 0:
        plt.legend(loc=2, fontsize=16)
        plt.title('Linear Amplification Only', fontsize=20)
        letter = '(a)'
    
    elif i == 1:
        plt.legend((h1[0], h2[0]), ('Seyhan & Stewart (2014)', 'AS1170.4-2007'), loc=2, fontsize=16)
        plt.title('Non-Linear Amplification', fontsize=20)
        letter = '(b)'
    
    '''    
    for j, v in enumerate(rev_vs30):	
        # plt Fa, Fv
        if i == 0:
            plt.semilogx(FaT, [Fa01[j], Fa01[j]], '--', lw=4.,color=[cs[j][0],cs[j][1],cs[j][2]])
            plt.semilogx(FvT, [Fv01[j], Fv01[j]], '--', lw=4.,color=[cs[j][0],cs[j][1],cs[j][2]])
        else:                                              
            plt.semilogx(FaT, [Fa04[j], Fa04[j]], '--', lw=4.,color=[cs[j][0],cs[j][1],cs[j][2]])
            plt.semilogx(FvT, [Fv04[j], Fv04[j]], '--', lw=4.,color=[cs[j][0],cs[j][1],cs[j][2]])
    '''
            
    # set new xticls
    ax.set_xticks([0.1, 0.3, 1., 3., 10])
    xlabels = ['0.1', '0.3', '1.0', '3.0', '10']
    ax.set_xticklabels(xlabels)
    plt.xlim([0.05, 5.])
    
    xpos = get_log_xy_locs((0.05, 5.), .98)
    plt.text(xpos, 6*.98, letter, va='top', ha='right', fontsize=20)
        
plt.savefig('cmp_AS1170_SS14_factors_BC.png', format='png', bbox_inches='tight')

plt.show()