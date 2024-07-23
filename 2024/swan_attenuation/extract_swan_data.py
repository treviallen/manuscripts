
from data_fmt_tools import get_swan_data
from obspy import UTCDateTime


################################################################################
# parse events
################################################################################
# id,isotime,latitude,longitude,depth,magnitude

events = []
lines = open('swan.dd.dec2023.csv').readlines()[1:]

for line in lines:
    dat = line.strip().split(',')
    tmpdict = {'id': dat[0], 'datetime': UTCDateTime(dat[1]),
               'lat': float(dat[2]), 'lon': float(dat[3]), 
               'dep': float(dat[4]), 'mag': float(dat[5])}
    events.append(tmpdict)
    
################################################################################
# loop thru events and get data
################################################################################

for ev in events: #[40:]:
    dt = ev['datetime']
    
    if dt.year > 2020 and ev['mag'] >= 2.95:
        dateTuple = (dt.year, dt.month, dt.day, dt.hour, dt.minute)
        
        get_swan_data(dateTuple, durn=600, network='2P', station='*', channel='*')
        