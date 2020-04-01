from pykml import parser
from mapping_tools import distance
from numpy import arange, array, where, mean

csvfile = 'usb000ajek.csv'

lines = open(csvfile).readlines()[7:]
#root = parser.fromstring(kmllines)

mmi = []
lat = []
lon = []
i = 0
for line in lines:
    dat = line.strip().split(',')
    mmi.append(float(dat[2]))
    lat.append(float(dat[3]))
    lon.append(float(dat[4]))
	
	
mmi = array(mmi)
lat = array(lat)
lon = array(lon)

# grid data
eqlat = -38.279
eqlon = 146.269
ztor = 17.24

dd = 0.05
d2 = dd / 2.

lonrng = arange(eqlon-4., eqlon+4., dd)
latrng = arange(eqlat-4., eqlat+4., dd)

cnt = 0
mmitxt = ''
obstxt = 'Source: Geoscience Australia\n'
for lo in lonrng:
    for la in latrng:
        idx = where((lon >= lo-d2) & (lon < lo+d2) & \
                    (lat >= la-d2) & (lat < la+d2))[0]
        if len(idx) > 1:
        	cnt += 1
        	meanlo = mean(lon[idx])
        	meanla = mean(lat[idx])
        	mmitxt += ','.join((str(meanlo), str(meanla), str(mean(mmi[idx])))) + '\n'
        	obstxt += '\t'.join((str(mean(mmi[idx])), str(meanla), str(meanlo), 'Unknown')) + '\n'
        	
outfile = 'USGS_MMI_0.05d.csv'
f = open(outfile, 'wb')
f.write(mmitxt)
f.close()

# write SM raw
outfile = 'USGS_mmi_obs.dat'
f = open(outfile, 'wb')
f.write(obstxt)
f.close()

        	

