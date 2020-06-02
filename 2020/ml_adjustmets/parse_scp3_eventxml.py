from misc_tools import dictlist2array, listdir_extension
from os import path
from numpy import nan, nanmedian, arange
import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.style.use('classic')

def strip_val(line):
    value = line.strip().split('>')[1].split('</')[0]
    return float(value)

def strip_val_str(line):
    return line.strip().split('>')[1].split('</')[0]

#scpxmlfile = 'bulk_sc3xml/ga2019aaxinx_sc3_event_parameters.xml'


evla = []
evlo = []
evdep = []
evmag = []
evmag_uncert = []
nstas = []
evmag_type = []

evdict = []

# get scp3 xml files
scpxmlfiles = listdir_extension('bulk_sc3xml', 'xml')

# loop thru files
for scpxmlfile in scpxmlfiles:
    print(scpxmlfile)
    lines = open(path.join('bulk_sc3xml', scpxmlfile)).readlines()
    fillParams = True
    
    # loop thru lines in xml file
    for i, line in enumerate(lines):
        # get location
        if line.strip().startswith('<latitude>'):
            ev = {'lat':strip_val(lines[i+1])}
            	
        if line.strip().startswith('<longitude>'):
            ev['lon'] = strip_val(lines[i+1])
            
        if line.strip().startswith('<depth>'):
            ev['dep'] = strip_val(lines[i+1])
        
        # get magnitude
        if line.strip().startswith('<magnitude publicID') and fillParams == True:
            ev['mag'] = strip_val(lines[i+2])
            if lines[i+3].strip().startswith('<uncertainty>'):
                ev['nstas'] = strip_val(lines[i+7])
                ev['mag_uncert'] = strip_val(lines[i+3])
                if ev['mag_uncert'] == 0.:
                    ev['mag_uncert'] = nan
                ev['mag_type'] = strip_val_str(lines[i+5])
            else:
                ev['nstas'] = strip_val(lines[i+6])
                ev['mag_uncert'] = nan
                ev['mag_type'] = strip_val_str(lines[i+4])
            
            fillParams = False
    
    # only add MLa with uncert
    if 'mag_uncert' in ev.keys():
        if ev['mag_type'] == 'MLa':
            evdict.append(ev)

###############################################################################
# get array and plot hist
###############################################################################
evmag_uncert = dictlist2array(evdict, 'mag_uncert')
print('median =', nanmedian(evmag_uncert))

plt.hist(evmag_uncert, bins=arange(0.0, 1.0, 0.05), facecolor='0.7')

plt.xlabel('MLa Uncertainty', fontsize=16)
plt.ylabel('Count', fontsize=16)
plt.xlim([0, 0.9])
plt.show()

