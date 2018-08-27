from __future__ import unicode_literals
from numpy import array, exp, log, interp, ceil, around, arange
from tools.oq_tools import return_annualised_haz_curves
import matplotlib.pylab as plt
from os import path, mkdir
import warnings, sys
from gmt_tools import cpt2colormap 
from misc_tools import remove_last_cmap_colour, get_log_xy_locs


reload(sys) # for unicode chars
sys.setdefaultencoding("latin-1")
warnings.filterwarnings("ignore")

plt_places = ['Adelaide', 'Canberra', 'Melbourne', 'Perth']
cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//Paired_08.cpt'
#cptfile = '//Users//tallen//Documents//DATA//GMT//cpt//qual-dark-06.cpt'
ncolours = 9
cmap, zvals = cpt2colormap(cptfile, ncolours)
cmap = remove_last_cmap_colour(cmap)
p = -1

###############################################################################
# read config file
###############################################################################

def split_config_lines(line):
    return line.split('=')[1].split('#')[0].strip().split(';')

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
        
    ###############################################################################
    # make colormap
    ###############################################################################
    
    #ncolours = len(jobs)+1
    ##cmap, zvals = cpt2colormap(cptfile, ncolours, rev=False)
    #cmap = plt.cm.get_cmap('hsv', ncolours)
    cs = (cmap(arange(ncolours-1)))
    
    ###############################################################################
    # def to get hazard curve data
    ###############################################################################
    
    def get_oq_haz_curves(hazcurvefile):
        if hazcurvefile.endswith('xml'):
            # Change the number 0.5 to 0.4 in hazard_curve-mean.xml so that it will run with the built-in parser.
            lines = open(hazcurvefile, 'r').readlines()
            lines[2] = 'xmlns="http://openquake.org/xmlns/nrml/0.4"\n'
            out = open(hazcurvefile, 'w')
            out.writelines(lines)
            out.close()
        
        # get annualize the curves
        curves, curlon, curlat, metadata = return_annualised_haz_curves(hazcurvefile)
        imls = array(metadata['imls'])
        
        return curves, curlon, curlat, metadata, imls
    
    i = 1
    ii = 0
    fig = plt.figure(i, figsize=(9, 10))
    yhaz2 = 1./2475.
    yhaz10 = 1./475.
    
    ###############################################################################
    # parse first job file to define plotting order
    ###############################################################################
    
    # make path to hazcurvefile
    #hazcurvefile1 = path.join(relativepaths[0], ''.join(('hazard_curve-mean_',jobs[0],'-SA(',period,').xml'))) # will need to be changed with OQ V2.2
    rootpath = path.split(conf_file)[0]
    hazcurvefile1 = path.join(rootpath, relativepaths[0], 'hazard_curve-mean_1.csv')
    
    # get data from first job
    curves1, curlon1, curlat1, metadata1, imls1 = get_oq_haz_curves(hazcurvefile1)
    
    
    
    # loop thru sites in first job file and plot
    for lon1, lat1, curve1 in zip(curlon1, curlat1, curves1):
        pltTrue = False
        #ii += 1
        ax = plt.subplot(1,1,1)
                    
                
        ###############################################################################
        # loops thru places to get title - check if want to plot
        ###############################################################################
        for place, plon, plat in zip(places, place_lon, place_lat):
            if around(plon, decimals=2) == around(lon1, decimals=2) \
               and around(plat, decimals=2) == around(lat1, decimals=2):
                
                # now loop through places we want to plot
                for pp in plt_places:
                    if pp == place:
                        
                        label_place = place
                        
                        # plot first curve
                        p += 1
                        h1 = plt.semilogy(imls1, curve1, color=cs[p*2], lw=2.0, label=label_place+' (F)')
                        
                        pltTrue = True
        
    
        ###############################################################################
        # loops thru subsequent job files
        ###############################################################################
    
        for j, hazlab in enumerate(hazcurvelabels[1:]):
            
            # make path for subsequent jobs
            '''
            root = path.split(conf_file)[0]
            filestr = ''.join(('hazard_curve-mean-SA(',period,')_',jobs[jj].strip(),'.csv'))
            hazcurvefilex = path.join(root, relativepaths[jj].strip(), filestr)
            '''
            hazcurvefilex = path.join(rootpath, relativepaths[1].strip(), 'hazard_curve-mean_1.csv')
            
            # get data from subsequent jobs
            curvepath = path.join(conf_file.split(path.sep)[0:-1])
            	
            curvesx, curlonx, curlatx, metadatax, imlsx = get_oq_haz_curves(hazcurvefilex)
    
            # loop thru Xnd OQ curves in job
            for lonx, latx, curvex in zip(curlonx, curlatx, curvesx):
                
                # if matches lon1 & lat1
                if around(lonx, decimals=2) == around(lon1, decimals=2) \
                   and around(latx, decimals=2) == around(lat1, decimals=2):
                    
                    # plt haz curves
                    if pltTrue == True:
                        hx = plt.semilogy(imlsx, curvex, '--', color=cs[p*2+1], lw=2.0, label=label_place+' (NF)')
        
###############################################################################
# make plot pretty
###############################################################################
plt.semilogy([0, 2.5], [yhaz2, yhaz2], 'k--')
plt.semilogy([0, 2.5], [yhaz10, yhaz10], 'k--')
yoff = get_log_xy_locs([1e-4, .1], .015)
plt.text(0.245, yhaz10+yoff, '1/475-year AEP', va='bottom',ha='right',fontsize=16)
yoff = get_log_xy_locs([1e-4, .1], .005)
plt.text(0.245, yhaz2+yoff/5., '1/2475-year AEP', va='bottom',ha='right',fontsize=16)
plt.legend()

plt.grid(which='both')
 
# get x lims from haz curve 1
thaz = exp(interp(log(1e-4), log(curve1[::-1]), log(imls1[::-1])))

# round to neares t 0.1
xmax = ceil(thaz / 0.1)  * 0.1
plt.xlim([0, xmax])
plt.ylim([1e-4, .1])
plt.xlim([0, .25])
ax.tick_params(labelsize=14)

plt.ylabel('Annual Probabability of Exceedance', fontsize=16)

plt.xlabel(' '.join(('Mean', period, 'Hazard (g)')), fontsize=16)
                       
# adjust x axes
#fig.subplots_adjust(hspace=0.2)

# save
if period == 'PGA' or period == 'PGV':
    plt.savefig(path.join(outputdir, '_'.join(('fault_nofault', period, jobsstr,str(i)+'.png'))), format='png',bbox_inches='tight')
else:
    plt.savefig(path.join(outputdir, '_'.join(('fault_nofault','SA('+period+')',jobsstr,str(i)+'.png'))), format='png',bbox_inches='tight')

  
plt.show()
    
    
