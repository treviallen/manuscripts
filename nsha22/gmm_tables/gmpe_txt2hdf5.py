# -*- coding: utf-8 -*-
"""
Created on Thu Dec 24 11:40:40 2015

@author: tallen
"""

from misc_tools import listdir_extension
from os import system, path, remove
from sys import argv

# get txt gmm files
infolder = argv[1]
outfolder = argv[2]

#folder = 'gmm_txt_tables'
extension = 'txt'
gmpetables = listdir_extension(infolder, extension)

'''
run gmpe_txt2hdf5.py weighted_geometric_ngae_tables weighted_geometric_ngae_tables

NOTE: need to delete line from
'''

for tab in gmpetables:

    '''
    # get params for converting to hdf5
    if tab.startswith('ENA') or tab.startswith('Wcrust_') \
       or tab.startswith('Winslab') or tab.startswith('Woffshore'):
        dist = 'Rhypo'
    elif tab.startswith('WcrustFRjb'):
        dist = 'Rjb'
    elif tab.startswith('WinterfaceCombo'):
        dist = 'Rrup'
    else:
        dist = 'Rrup' # for NGA-E
    '''
        
    tabin  = path.join(infolder, tab)
    
    hdf5out = path.join(outfolder, path.split(tab)[-1].strip('txt')+'hdf5')
    
    try:
        remove(hdf5out)
    except:
        print('\nContinue: hdf5 file does not exist')
    
    # read in table and get diatance metric
    lines = open(tabin).readlines()
    distType = lines[0].split('. Log10')[0].split()[-1].capitalize()
    print('\n' + distType)
    
    system(''.join(('python /Users/trev/Documents/Code/my_codes/gmpe_table_builder.py --input-file=', \
                    tabin, ' --output-file=', hdf5out, ' --distance-key=', distType, ' --is-log True')))
    
    

