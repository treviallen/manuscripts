from misc_tools import dictlist2array, listdir_extension, get_mpl2_colourlist, \
                       get_binned_stats
from os import path
from numpy import nan, nanmedian, arange
import matplotlib.pyplot as plt
from obspy import UTCDateTime

import matplotlib as mpl
mpl.style.use('classic')

def strip_val(line):
    value = line.strip().split('>')[1].split('</')[0]
    return float(value)

def strip_val_str(line):
    return line.strip().split('>')[1].split('</')[0]

#scpxmlfile = 'bulk_sc3xml/ga2019aaxinx_sc3_event_parameters.xml'

# get scp3 xml files
scpxmlfiles = listdir_extension('bulk_sc3xml', 'xml')

evdict = []
# loop thru files
for scpxmlfile in scpxmlfiles:
    print(scpxmlfile)
    lines = open(path.join('bulk_sc3xml', scpxmlfile)).readlines()
    fillParams = True
    ev = {}
    
    # loop thru lines in xml file
    for i, line in enumerate(lines):
        # get location
        if line.strip().startswith('<origin publicID'):
            ev['origin'] = UTCDateTime(strip_val_str(lines[i+2]))
            	
        if line.strip().startswith('<latitude>'):
            ev['lat'] = strip_val(lines[i+1])
            	
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
fig = plt.figure(1, figsize=(11,11))
evmag_uncert = dictlist2array(evdict, 'mag_uncert')
evmag = dictlist2array(evdict, 'mag')
origin = dictlist2array(evdict, 'origin')
allmags = evmag
alluncert = evmag_uncert
print(min(origin))
print(max(origin))
print('median =', nanmedian(evmag_uncert))

ax = plt.subplot(221)
plt.hist(evmag_uncert, bins=arange(0.0, 1.0, 0.05), facecolor='0.7')

#plt.xlabel('MLa Uncertainty', fontsize=16)
plt.ylabel('Count', fontsize=16)
plt.xlim([0, 1.])

# add label
props = dict(boxstyle='round', facecolor='w', alpha=1)
text = ' '.join((r'$\tilde{x}_{All}$', '=', str('%0.2f' % nanmedian(evmag_uncert))+'\n',\
                 r'$\it{n}$', '=', str(len(evmag_uncert)),''))
ylims = ax.get_ylim()
plt.text(0.95, ylims[1]*0.95, text, ha='right', va='top', fontsize=16, bbox=props)
plt.text(0.02*0.9, ylims[1]*0.98, '(a)', ha='left', va='top', fontsize=16)

###############################################################################
# get array and plot hist
###############################################################################

import shapefile
from shapely.geometry import Point, Polygon
from mapping_tools import get_field_data

shpfile = 'shapefile/australia_ml_regions.shp'
sf = shapefile.Reader(shpfile)
shapes = sf.shapes()
polygons = []
for poly in shapes:
    polygons.append(Polygon(poly.points))
    
ml_reg = get_field_data(sf, 'ML_REGION', 'str')
mmin = 2.45

i = 2
letter = ['(b)', '(c)', '(d)']
reg_dict = []
for poly, reg in zip(polygons, ml_reg):
    reg_evdict = []
    for ev in evdict:
        pt = Point(ev['lon'], ev['lat'])
        if pt.within(poly) and ev['mag'] >= mmin:
            reg_evdict.append(ev)
            
    # now plot regions
    evmag_uncert = dictlist2array(reg_evdict, 'mag_uncert')
    evmag = dictlist2array(reg_evdict, 'mag')
    
    print('median =', nanmedian(evmag_uncert))
    
    ax = plt.subplot(2,2,i)
    plt.hist(evmag_uncert, bins=arange(0.0, 1.0, 0.05), facecolor='0.7')
    
    if i > 2:
        plt.xlabel('$\mathregular{M_L}$ Uncertainty', fontsize=16)
    if i == 3:
        plt.ylabel('Count', fontsize=16)
    plt.xlim([0, 1.])
    
    text = ' '.join((r'$\tilde{x}_{'+reg+'}$', '=', str('%0.2f' % nanmedian(evmag_uncert))+'\n',\
                     r'$\it{n}$', '=', str(len(evmag_uncert)),''))
    ylims = ax.get_ylim()
    plt.text(0.95, ylims[1]*0.95, text, ha='right', va='top', fontsize=16, bbox=props)
    plt.text(0.02*0.9, ylims[1]*0.98, letter[i-2], ha='left', va='top', fontsize=20)
   
    i += 1
    
    reg_tmp = {'mag_uncert': evmag_uncert, 'mag': evmag, 'reg': reg}
    print (len(evmag_uncert))
    print (len(evmag))
    reg_dict.append(reg_tmp)

plt.savefig('ml_uncertainty_histograms.png', fmt='png', bbox_inches='tight')
plt.savefig('ml_uncertainty_histograms.eps', fmt='eps', bbox_inches='tight')
plt.show()

###############################################################################
# plot mag uncert with mag
###############################################################################

fig = plt.figure(2, figsize=(10,5))
sym = ['^', 'D', 'o']
cols = get_mpl2_colourlist()
for i, reg in enumerate(reg_dict):
    plt.plot(reg['mag'], reg['mag_uncert'], sym[i], mec=cols[i], mfc='none', mew=1.5, label = reg['reg'])

# get binned stats
bins = arange(2.5, 6.1, 0.1)
med_uncert, std_uncert, medx, bins, nperbin = get_binned_stats(bins, allmags, alluncert)
plt.errorbar(bins, med_uncert, yerr=std_uncert, fmt='rs', label='Binned Median')  
plt.legend(loc=1, numpoints=1)
plt.ylim([0, 1])
plt.xlim([2.4, 5])

plt.savefig('ml_uncertainty_vs_ml.png', fmt='png', bbox_inches='tight')    
plt.savefig('ml_uncertainty_vs_ml.eps', fmt='eps', bbox_inches='tight') 
plt.show()