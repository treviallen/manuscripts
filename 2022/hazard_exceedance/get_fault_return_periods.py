from openquake.hmtk.parsers.source_model.base import BaseSourceModelParser
from openquake.hmtk.parsers.source_model.nrml04_parser import nrmlSourceModelParser
from openquake.hazardlib.nrml import to_python 
from numpy import array
from mapping_tools import distance

# middle coords for Willunga Fault
will_lat = -35.3453885334
will_lon = 138.459539572 
target_mag = 6.9

fault_nrml = '/Users/trev/Documents/Geoscience_Australia/NSHA2018/source_models/faults/National_Fault_Source_Model_2018_GR/National_Fault_Source_Model_2018_GR_source_model.xml'

#src_parser = BaseSourceModelParser(fault_nrml)
src_parser = nrmlSourceModelParser(fault_nrml)
fsm = src_parser.read_file("AU FSM")

fault_dict = []
for fs in fsm.sources:
    addFault = False
    
    # get GR-MFD
    rates = array(fs.mfd.get_annual_occurrence_rates())
    
    # get trace
    flons = []
    flats = []
    tp = str(fs.fault_trace.points).split(', ')
    for i, p in enumerate(tp):
       if p.replace('[','').startswith('<Latitude='):
           flats.append(float(tp[i].split('=')[-1]))
           flons.append(float(tp[i+1].split('=')[-1]))
    
    # extract faults within 100 km
    minDist = 9999.
    for fla, flo in zip(flats, flons):
        rngkm = distance(fla, flo, will_lat, will_lon)[0]
        if rngkm <= 100.:
            addFault = True
    
    
    if addFault == True:
        tmp = {'name':fs.name, 'mags':rates[:,0], 'rates':rates[:,1], 'lats':array(flats), 'lons':array(flons)}
        
        fault_dict.append(tmp)
        
	
###############################################################################
# get rates above target MW

sum_rates = 0

for fd in fault_dict:
    idx = fd['mags'] > target_mag
    sum_rates += sum(fd['rates'][idx])

        
print(sum_rates)
print('\n')
print(1./sum_rates)