from mapping_tools import reckon, distance

lon1 = 131.0
lat1 = -11.5
lon2 = 134.5
lat2 = -17.0

rngkm, az, baz = distance(lat1, lon1, lat2, lon2)

# get inc distance
npts = 16
inckm = rngkm / (npts-1)

plats = []
plons = []
csvtxt = ''

for i in range(0, npts):

    lon2d, lat2d = reckon(lat1, lon1, i*inckm, az)
    plats.append(lat2d)
    plons.append(lon2d)
    
    csvtxt += ','.join((str('%0.4f' % lon2d), str('%0.4f' % lat2d))) + '\n'
    
f = open('north_aus_profile.csv', 'wb')
f.write(csvtxt)
f.close()
