# -*- coding: utf-8 -*-
"""
Created on Thu Apr 30 11:47:02 2015

convert gridded hazard to map

Usage:
    python map_nsha18.py <path to csv file>
    

@author: tallen
"""
from sys import argv
from matplotlib.mlab import griddata
from matplotlib import colors, colorbar #, cm
from os import path, mkdir, getcwd
#import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from numpy import arange, array, log10, mean, mgrid, ogrid, percentile, ma, isnan, nan, where, delete, floor, ceil
from mapping_tools import distance # drawshapepoly, labelpolygon, 
import shapefile
from scipy.constants import g
import matplotlib.gridspec as gridspec

#from gmt_tools import cpt2colormap
#from shapely.geometry import Point, Polygon

##############################################################################
# set some default values here
##############################################################################
mpl.use('Agg')
mpl.style.use('classic')
mpl.rcParams['pdf.fonttype'] = 42

drawshape = False # decides whether to overlay seismic sources

# set map resolution
res = 'i' 

# get current working directory (and computer!)
cwd = getcwd()

##############################################################################
# define inputs
##############################################################################

# set with paths to hazard grids - format = plotting name; file path
paramfile = argv[1]

# which probability - acceptable values are: 2 (2%), 9 (9.5%) or 10 (10%)
pltProbability = argv[2]

##############################################################################
# parse param file
##############################################################################

lines = open(paramfile).readlines()
hazfiles = []
modnames = []
outfile = lines[0].strip()
for line in lines[1:]:
    modnames.append(line.strip().split(';')[0])
    hazfiles.append(line.strip().split(';')[1])

# set figure size based on number of models
rows = 1
figure = plt.figure(1,figsize=(10,5*rows))

trenchlat = -9.
trenchlon = 129.5

pltlett = ['a)', 'b)', 'c)', 'd)', 'e)']
collist = ['orange', 'dodgerblue']
##############################################################################
# loop thru files and parse hazard grids
##############################################################################
ax = figure.add_subplot(111)
for ii, hazfile in enumerate(hazfiles):

    # parse hazard grid file 
    lines = open(hazfile).readlines()
    
    # get keys for model
    if lines[0].startswith('#'):
        line = lines[1]
    else:
        line = lines[0]
    
    # get dictionary keys
    keys = line.strip().split(',')[2:]
    
    # make grid dictionary
    profdict = []

    print '\nReading', modnames[ii]
    for line in lines[2:]:
        dat = line.strip().split(',')
        
        tmpdict = {'lon':float(dat[0]), 'lat':float(dat[1])}
        
        # fill keys
        idx = 2
        for key in keys:
            tmpdict[key] = float(dat[idx])
            idx += 1
        
        # add to grid list
        profdict.append(tmpdict)
        
    ##############################################################################    
    # get key index for plotting
    ##############################################################################
    
    for i, key in enumerate(keys):
        keyProb = str(int(floor(100*float(key.split('-')[-1]))))
        if keyProb == pltProbability:
            mapidx = i
    
    ##############################################################################    
    # now make profile
    ##############################################################################
    
    #keys = ['PGA_10', 'PGA_02', 'SA02_10', 'SA02_02', 'SA10_10', 'SA10_02']
    for i, key in enumerate([keys[mapidx]]): # just plot 1 for now!
        print key
    
        # get IM period
        period = key.split('-')[0]
        period = period.replace('(','')
        period = period.replace(')','')
        period = period.replace('.','')
        
        # get map probability of exceedance
        probFraction = str(float(key.split('-')[-1]))
        probability = str(100*float(key.split('-')[-1])).split('.')[0]+'%'
        #probability = str(100*float(key.split('-')[-1]))+'%'
        if probability == '9%':
            probability = '9.5%'
        print 'Probability', probability
        
        # add buffer to data
        # build data to plot
        hazvals = []
        latlist = []
        lonlist = []
        distlist = []
        for profval in profdict:
            lonlist.append(profval['lon'])
            latlist.append(profval['lat'])
            rngkm, az, baz = distance(trenchlat, trenchlon, profval['lat'], profval['lon'])
            distlist.append(rngkm)
            
            if profval[key] == 0.0:
                hazvals.append(0.0)
            else:
                hazvals.append(profval[key])
        print hazvals
        # now plot
        plt.plot(distlist, hazvals, 'o-', c=collist[ii], ms=8, lw=2.0, label=modnames[ii])

plt.legend(loc=1, numpoints=1, fontsize=13)
plt.ylabel(' '.join(('Sa(1.0)', probability, 'in 50-Year Mean Hazard (g)')), fontsize=16)
plt.xlabel('Approximate Distance from Timor Trough (km)', fontsize=16)
plt.grid()

# now save png file
plt.savefig(path.join('multi_hazard_profile_'+outfile.replace(' ','_')+'.'+probFraction+'.png'), \
            dpi=300, format='png', bbox_inches='tight')

plt.show()