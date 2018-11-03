import matplotlib.pyplot as plt
from numpy import arange, array, interp
from os import path, getcwd, sep
from sys import argv
import matplotlib as mpl

mpl.style.use('classic')
mpl.rcParams['pdf.fonttype'] = 42

paramfile = argv[1] # param file with locs folder where fractile files sit

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
        citycsv = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/shared/nsha_cities.csv'
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

fracDict1, keys = parse_plot_fractiles(fracpaths[0])
fracDict2, keys = parse_plot_fractiles(fracpaths[1])
    
###################################################################################
# let's make the plots
###################################################################################

places = ['Perth', 'Darwin', 'Adelaide', 'Melbourne', 'Hobart', 'Canberra', 'Sydney', 'Brisbane']
letters = ['(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)','(i)']


# loop through keys

# just do PGA
for k, key in enumerate(keys[:2]):
    # set up figure
    fig = plt.figure(k+1, figsize=(15, 15))
    
    # loop thru places to plot
    for i, place in enumerate(places):
        # loop thru fracDict
        for frac1, frac2 in zip(fracDict1, fracDict2):
            if place == frac1['place']:
                print frac1['place'], frac2['place']
            
                # plot fig
                ax = plt.subplot(4, 2, i+1)
                plt.semilogx(frac1['quant_'+key], fractiles, '-', c='g', lw=1.5, label=modnames[0])
                plt.semilogx(frac2['quant_'+key], fractiles, '-', c='b', lw=1.5, label=modnames[1])
                
                
                # make pretty
                plt.title(place, fontsize=15)
                #plt.text(0.095, 0.02, place, fontsize=18, va='bottom', ha='right')
                plt.text(0.095, 0.98, letters[i], fontsize=18, va='top', ha='right')
                if i == 0 or i == 2 or i == 4 or i == 6:
                    plt.ylabel('Fractile', fontsize=16)
                    
                
                if i >= 6:
                    plt.xlabel(key.replace('(','').replace(')','').split('-')[0] + ' (g)', fontsize=16)
                    
                
                plt.grid(which='both')
                plt.xlim([0.003, 0.1])
                
                
                # plt mean
                plt.semilogx([frac1['mean_'+key],frac1['mean_'+key]], [0,1], '-', c='orange', lw=1.5,   label=modnames[0]+' Mean')
                plt.semilogx([frac2['mean_'+key],frac2['mean_'+key]], [0,1], '--', c='r', lw=2.5, label=modnames[1]+' Mean')
                
                if i == 0:
                    plt.legend(loc=2, fontsize=11)
                
                '''
                # plt median
                plt.semilogx([frac['quant_'+key][50],frac['quant_'+key][50]], [0,1], '--', c='dodgerblue', lw=1.5, label='50th Percentile')
                # plt 84th
                plt.semilogx([frac['quant_'+key][84],frac['quant_'+key][84]], [0,1], '--', c='orange', lw=1.5, label='84th Percentile')
                # plt 95th
                plt.semilogx([frac['quant_'+key][95],frac['quant_'+key][95]], [0,1], '--', c='r', lw=1.5, label='95th Percentile')
                '''
                
                
    #plt.suptitle(fracFolder.split(sep)[1] + ' ' + key, fontsize=20)

    # set fig file
    figFile = '_'.join((path.join('cdf',outfile),key,'CDF.png'))
    plt.savefig(figFile, fmt='png', bbox_inches='tight')
    
    
    
plt.show()