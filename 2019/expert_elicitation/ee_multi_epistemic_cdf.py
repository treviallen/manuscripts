import matplotlib.pyplot as plt
from numpy import arange, array, interp, mean
from os import path, getcwd, sep
from sys import argv
import matplotlib as mpl
from gmt_tools import cpt2colormap 
from misc_tools import remove_last_cmap_colour, get_log_xy_locs
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as PathEffects

mpl.style.use('classic')
mpl.rcParams['pdf.fonttype'] = 42
plt.rc('xtick',labelsize=13)
plt.rc('ytick',labelsize=13)


paramfile = argv[1] # param file with locs folder where fractile files sit 

'''
run cmp_multi_epistemic_cdf.py fractiles_NSHA12_regional_background.param
'''

##############################################################################
# parse param file
##############################################################################

lines = open(paramfile).readlines()
fracpaths = []
modnames = []
outfile = lines[0].strip()
for line in lines[1:]:
    modnames.append(line.strip().split(';')[0])
    fracpaths.append(line.strip().split(';')[1])
    
fractiles = arange(0., 1.01, 0.01)

##############################################################################
# parse cpt
##############################################################################

# get colours
if getcwd().startswith('/nas'):
    cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/gay-flag-1978.cpt'
    #cptfile = '/nas/active/ops/community_safety/ehp/georisk_earthquake/hazard/DATA/cpt/sst.cpt'
else:
    cptfile = '/Users/trev/Documents/DATA/GMT/cpt/gay-flag-1978.cpt'
ncolours = len(modnames)+1
#ncolours = 8
cmap, zvals = cpt2colormap(cptfile, ncolours)
cmap = remove_last_cmap_colour(cmap)
cs = (cmap(arange(ncolours-1)))
ls = '-'
lw = 1.75

if ncolours < 4:
    cs = ['orangered', 'dodgerblue']
    ls = ['-', '--']

##############################################################################
# call parsing & plotting function
##############################################################################

def parse_plot_fractiles(fracfolder):
    # set params 
    fractiles = arange(0., 1.01, 0.01)
    
    # get keys to fill
    fracFile = path.join(fracfolder, 'quantile_map-0.5_1.csv' )
    
    # parse file
    lines = open(fracFile).readlines()
    
    # get keys from medianfile
    keys = lines[1].strip().split(',')[2:]
    	
    # load mean lines to add
    meanFile = path.join(fracfolder, 'hazard_map-mean_1.csv' )
    meanLines = open(meanFile).readlines()[2:]
    
    # loop thru fractiles and fill city Dict
    fracDict = []
    for i, fractile in enumerate(fractiles):
        fracStr = str('%0.2f' % fractile)
        if fracStr.endswith('0'):
            fracStr = fracStr[0:-1]
        fracFile = path.join(fracfolder, 'quantile_map-'+fracStr+'_1.csv' )
        
        # parse file
        lines = open(fracFile).readlines()[2:]
        
        # loop through lines
        j = 0
        for line, mLine in zip(lines, meanLines):
            dat = line.strip().split(',')
            mdat = mLine.strip().split(',')
            
            # if first fractile
            if i == 0:
                tmpdict = {'lon':float(dat[0]), 'lat':float(dat[1])}
                for k, key in enumerate(keys):
                    tmpdict['quant_'+key] = [float(dat[k+2])]
                    
                    # add mean values too
                    tmpdict['mean_'+key]  = float(mdat[k+2])
    
                fracDict.append(tmpdict)
                
            # now append to existing dict
            else:
                for k, key in enumerate(keys):
                    fracDict[j]['quant_'+key].append(float(dat[k+2]))
                    
            j+=1
                    
    # now turn lists to arrays
    for j in range(0, len(fracDict)):
        for k, key in enumerate(keys):
            fracDict[j]['quant_'+key] = array(fracDict[j]['quant_'+key])
            
    ###################################################################################
    # match city name to fracDict
    ###################################################################################
    
    # first parse city file
    cwd = getcwd()
    if cwd.startswith('/Users'): #mac
        citycsv = '/Users/trev/Documents/Geoscience_Australia/NSHA2018/shared/nsha_cities.csv'
    else:
        citycsv = '/nas/active/ops/community_safety/ehp/georisk_earthquake/modelling/sandpits/tallen/NSHA2018/shared/nsha_cities.csv'
    lines = open(citycsv).readlines()
        
        
        
    # make city dict
    cityDict = []
    for line in lines:
        dat = line.strip().split(',')
        tmpdict = {'city':dat[2], 'lon':float(dat[0]), 'lat':float(dat[1])} 
        cityDict.append(tmpdict)
    
    
    for j in range(0, len(fracDict)):
        for city in cityDict:
            if city['lon'] == fracDict[j]['lon'] \
               and city['lat'] == fracDict[j]['lat']:
               
               # add place
               fracDict[j]['place'] = city['city']
    
    return fracDict, keys
    
###################################################################################
# start main code here
###################################################################################
fracDict, keys = parse_plot_fractiles(fracpaths[0])
    
###################################################################################
# let's make the plots
###################################################################################

altPlaces = False

if altPlaces == False:
    places = ['Perth', 'Darwin', 'Adelaide', 'Melbourne', 'Hobart', 'Canberra', 'Sydney', 'Brisbane']
else:
    places = ['Wongan Hills', 'Kalgoorlie', 'Karratha', 'Yulara', 'Port Pirie', 'Hawker', 'Cooma', 'Morwell']

letters = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)']

gs1 = gridspec.GridSpec(4, 2)
gs1.update(wspace=0.095, hspace=0.2) # negative looks bad in "show", but ok in pngs

# just do PGA
for k, key in enumerate(keys[:1]):
    # set up figure
    fig = plt.figure(k+1, figsize=(15, 17))

    # loop thru places to plot
    for i, place in enumerate(places):
        #ax = plt.subplot(4, 2, i+1)
        ax = fig.add_subplot(gs1[i])
        
        # loop through models
        j = 0
        for modname, fracpath in zip(modnames, fracpaths):
         
            # loop thru fracDict
            fracDict, keys = parse_plot_fractiles(fracpath)
         
            for frac in fracDict:
                if place == frac['place']:
                    print(frac['place'], mean(frac['quant_'+key]))
                    if j == 4:
                        col = 'gold'
                    else:
                        col = cs[j]
                    # plot fig
                    if modname.endswith('(B)'):
                        plt.semilogx(frac['quant_'+key], fractiles, '-', c=col, lw=lw, label=modname)
                    else:
                        plt.semilogx(frac['quant_'+key], fractiles, '-', c=col, lw=lw, label=modname)
                    
                    j += 1        
        
        # make pretty
        plt.title(place, fontsize=15)
        
        if i == 0 or i == 2 or i == 4 or i == 6:
            plt.ylabel('Fractile', fontsize=16)
                  
        if i >= 6:
            plt.xlabel(key.replace('(','').replace(')','').split('-')[0] + ' (g)', fontsize=16)
            
                    
        plt.grid(which='both')
        
        ticks = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2]
        
        # reset ticks
        if i == 7 and i == 0:
            if altPlaces == True:
                plt.xlim([0.01, 0.3])        
                ticks = [0.01, 0.02, 0.05, 0.1, 0.2]
            else:
                plt.xlim([0.003, 0.2])
        else:
            plt.xlim([0.003, 0.2])
        
        # set letters
        if altPlaces == False:
            xplt = get_log_xy_locs(ax.get_xlim(), 0.02)
            txt = plt.text(xplt, 0.98, letters[i], fontsize=18, va='top', ha='left')
        else:
            xplt = get_log_xy_locs(ax.get_xlim(), 0.98)
            txt = plt.text(xplt, 0.98, letters[i], fontsize=18, va='top', ha='right')
        txt.set_path_effects([PathEffects.withStroke(linewidth=5, foreground='w')])

        ax.set_xticks(ticks)
        ax.set_xticklabels([str(x) for x in ticks])
        
        if altPlaces == False:
            if i == 7:
                plt.legend(loc=4, fontsize=14)
        else:
            if i == 0:
                plt.legend(loc=2, fontsize=13)
            
     # set fig file
    if altPlaces == True:
        figFile = '_'.join(('regional_'+outfile,key,'CDF.png'))
    else:
        figFile = '_'.join(('capital_'+outfile,key,'CDF.png'))
    plt.savefig(figFile, fmt='png', bbox_inches='tight')
        
    
    
plt.show()