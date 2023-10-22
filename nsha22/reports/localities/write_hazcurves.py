from __future__ import unicode_literals
from numpy import array, exp, log, interp, ceil, around, arange
from tools.oq_tools import return_annualised_haz_curves
import matplotlib.pylab as plt
from os import path, mkdir, getcwd
import warnings, sys
from gmt_tools import cpt2colormap 
from misc_tools import remove_last_cmap_colour, get_log_xy_locs, get_mpl2_colourlist, checkfloat
import matplotlib as mpl
mpl.style.use('classic')

###################################################################################
# load AS1170.4 table
###################################################################################

as1170table = '/Users/trev/Documents/Geoscience_Australia/NSHA2023/shared/nsha18_AS1170.4_cmp.csv'

as1170dict = []
lines = open(as1170table).readlines()
for line in lines[1:]:
    dat = line.strip().split(',')
    tmpdict = {'city':dat[0].strip(), 'lon':float(dat[1]), 'lat':checkfloat(dat[2]), 'as1170':checkfloat(dat[3]), \
    	         'nshm12_500':checkfloat(dat[4]), 'nsha18_475':float(dat[5]), 'nsha18_500':float(dat[6])} 
    as1170dict.append(tmpdict)

###################################################################################
# set files
###################################################################################
warnings.filterwarnings("ignore")

"""
#reload(sys) # for unicode chars
#sys.setdefaultencoding("latin-1")
plt_places = ['Adelaide', 'Canberra', 'Melbourne', 'Perth']
#plt_places = ['York', 'Meckering', 'Adelaide', 'Morwell']
if getcwd().startswith('/nas'):
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/Paired_08.cpt'
else:
    cptfile = '//Users//trev//Documents//DATA//GMT//cpt//Paired_08.cpt'
"""
'''
ncolours = 9
cmap, zvals = cpt2colormap(cptfile, ncolours)
cmap = remove_last_cmap_colour(cmap)
p = -1
'''
cols = get_mpl2_colourlist()
###############################################################################
# read config file
###############################################################################

def split_config_lines(line):
    return '='.join(line.split('=')[1:]).strip().strip('\n').split(';')
    #return line.split('=')[1].split('#')[0].strip('\n').split(';')
    
if __name__ == "__main__":
       
    # parse config file
    conf_file = sys.argv[1]
    
    # get paths for input files
    lines = open(conf_file).readlines()
    period         = split_config_lines(lines[0])
    jobs           = split_config_lines(lines[1])
    relativepaths  = split_config_lines(lines[2])
    hazcurvelabels = split_config_lines(lines[3])
    outputdir      = split_config_lines(lines[4])[0]
    sitelistfile   = split_config_lines(lines[5])
    plt_places     = split_config_lines(lines[6])
    
    #period = str(float(period[0])) # fix general weirdness
    period = period[0]
    jobsstr = '_'.join([x.strip() for x in jobs])
    hazcurvelabels = [x.strip().replace('\\n', '\n') for x in hazcurvelabels]
    
    # check to see if exists
    if path.isdir(outputdir) == False:
        mkdir(outputdir)

    ###############################################################################
    # parse site file
    ###############################################################################
    
    lines = open(sitelistfile[0]).readlines()
    places = []
    place_lat = []
    place_lon = []
    
    for line in lines:
        dat = line.strip().split(',')
        place_lon.append(float(dat[0]))
        place_lat.append(float(dat[1]))
        places.append(dat[2])
        
    """
    ###############################################################################
    # make colormap
    ###############################################################################
    
    #ncolours = len(jobs)+1
    ##cmap, zvals = cpt2colormap(cptfile, ncolours, rev=False)
    #cmap = plt.cm.get_cmap('hsv', ncolours)
    cs = (cmap(arange(ncolours-1)))
    """
    ###############################################################################
    # def to get hazard curve data
    ###############################################################################
    
    i = 1
    ii = 0
    fig = plt.figure(i, figsize=(16, 10))
    yhaz2 = 1./2475.
    yhaz10 = 1./475.
    
    ###############################################################################
    # parse first job file to define plotting order
    ###############################################################################
    '''
    # make path to hazcurvefile
    #hazcurvefile1 = path.join(relativepaths[0], ''.join(('hazard_curve-mean_',jobs[0],'-SA(',period,').xml'))) # will need to be changed with OQ V2.2
    rootpath = path.split(conf_file)[0]
    hazcurvefile1 = path.join(rootpath, relativepaths[0], 'hazard_curve-mean-PGA_1.csv')
    
    # get data from first job
    #curves1, curlon1, curlat1, metadata1, imls1 = get_oq_haz_curves(hazcurvefile1)
    siteDict1, imls, investigation_time = return_annualised_haz_curves(hazcurvefile1)
    '''
    rootpath = path.split(conf_file)[0]
    # loop thru sites in first job file and plot
    for p, pp in enumerate(plt_places):
        pltTrue = False
        #ii += 1
        ax = plt.subplot(2,3,p+1)
        
        # get data from jobs
        for j, relativepath in enumerate(relativepaths):
            '''
            # for testing 
            try:
            	hazcurvefile = path.join(relativepath, 'quantile_curve-0.5-PGA_1.csv')
            	siteDict, imls, investigation_time = return_annualised_haz_curves(hazcurvefile)
            	
            except:
            '''
            hazcurvefile = path.join(relativepath, 'hazard_curve-mean-PGA_1.csv')
            siteDict, imls, investigation_time = return_annualised_haz_curves(hazcurvefile)
            
                    
            # loop thru siteDict
            for sd in siteDict:
            
                ###############################################################################
                # loops thru places to get title - check if want to plot
                ###############################################################################
                for place, plon, plat in zip(places, place_lon, place_lat):
                    if around(plon, decimals=2) == around(sd['lon'], decimals=2) \
                       and around(plat, decimals=2) == around(sd['lat'], decimals=2):
                        
                        # now loop through places we want to plot
                        if pp == place:
                            
                            label_place = place
                            
                            # plot first curve
                            h = plt.loglog(imls, sd['poe_probs_annual'], color=cols[j], lw=1.5, label=hazcurvelabels[j])
                            
                            pltTrue = True
                            
                            plt.title(pp, fontsize=16)
                            
            ###############################################################################
            # plot fractiles
            ###############################################################################
            if relativepath.startswith('2023'):
                fractfile16 = path.join(relativepath, 'quantile_curve-0.05-PGA_1.csv')
                siteDict16, imls, investigation_time = return_annualised_haz_curves(fractfile16)
                
                fractfile84 = path.join(relativepath, 'quantile_curve-0.95-PGA_1.csv')
                siteDict84, imls, investigation_time = return_annualised_haz_curves(fractfile84)
                        
                # loop thru siteDict
                for sd16, sd84 in zip(siteDict16, siteDict84):
                
                    ###############################################################################
                    # loops thru places to get title - check if want to plot
                    ###############################################################################
                    for place, plon, plat in zip(places, place_lon, place_lat):
                        if around(plon, decimals=2) == around(sd16['lon'], decimals=2) \
                           and around(plat, decimals=2) == around(sd16['lat'], decimals=2):
                            
                            # now loop through places we want to plot
                            if pp == place:
                                
                                h = plt.loglog(imls, sd16['poe_probs_annual'], '--', color=cols[j], lw=1.5, label='5-95th Percentile')
                                h = plt.loglog(imls, sd84['poe_probs_annual'], '--', color=cols[j], lw=1.5)
                                
        ###############################################################################
        # make plot pretty
        ###############################################################################
        plt.loglog([0, 2.5], [yhaz2, yhaz2], 'k--')
        plt.loglog([0, 2.5], [yhaz10, yhaz10], 'k--')
        
        # plt haz floor
        plt.loglog(0.08, 1/500., 's', ms=7, c=cols[j+1], mec=cols[j+1], lw=1.5, label='AS1170.4 Floor')
        
        # add AS1170.4-2018
        for as1170 in as1170dict:
            if as1170['city'] == pp:
                plt.loglog(max([as1170['as1170'], 0.08]), 1/500., 'o', ms=7, c=cols[j+2], mec=cols[j+2], label='AS1170.4-2018')
        
        yoff = get_log_xy_locs([1e-4, .1], .015)
        #plt.text(0.295, yhaz10+yoff, '1/475 AEP', va='bottom',ha='right',fontsize=14)
        plt.text(0.0012, yhaz10+yoff, '1/475 AEP', va='bottom',ha='left',fontsize=12)
        yoff = get_log_xy_locs([1e-4, .1], .005)
        #plt.text(0.295, yhaz2+yoff/5., '1/2475 AEP', va='bottom',ha='right',fontsize=14)
        plt.text(0.0012, yhaz2+yoff/5., '1/2475 AEP', va='bottom',ha='left',fontsize=12)
        
        if p == 0:
            plt.legend(loc=1, fontsize=10, numpoints=1)
        
        plt.grid(which='both')
         
        # get x lims from haz curve 1
        thaz = exp(interp(log(1e-4), log(sd['poe_probs_annual'][::-1]), log(imls[::-1])))
        
        # round to neares t 0.1
        xmax = ceil(thaz / 0.1)  * 0.1
        plt.xlim([0, xmax])
        plt.ylim([1e-4, .1])
        plt.xlim([1e-3, .5])
        ax.tick_params(labelsize=10)
        
        if p == 0 or p == 3:
            plt.ylabel('Annual Probabability of Exceedance', fontsize=14)
        
        if p >= 3 or len(plt_places) <= 3:
            plt.xlabel(' '.join(('Mean', period, 'Hazard (g)')), fontsize=14)
                       
# adjust x axes
#fig.subplots_adjust(hspace=0.2)

# save
if period == 'PGA' or period == 'PGV':
    plt.savefig(path.join(outputdir, '_'.join((period, jobsstr,str(i)+'.png'))), format='png',bbox_inches='tight', dpi=300)
else:
    plt.savefig(path.join(outputdir, '_'.join(('SA('+period+')',jobsstr,str(i)+'.png'))), format='png',bbox_inches='tight', dpi=300)
  
plt.show()
    
    
