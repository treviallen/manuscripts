from netCDF4 import Dataset as NetCDFFile
from os import system

nc = NetCDFFile('avg_0.0s_500yr_180_30km_60km.grd')

pgas = nc.variables['z'][:]
lons = nc.variables['x'][:]
lats = nc.variables['y'][:]

xyztxt = ''
for i, lat in enumerate(lats):
    for j, lon in enumerate(lons): 
        xyztxt += ' '.join((str(lon), str(lat), str(pgas[i, j]))) + '\n'
    
print('Writing txt file...')
f = open('nshm12.xyz', 'w')
f.write(xyztxt)
f.close()

# make grd file
print('Making surface...')
system('gmt5 surface nshm12.xyz -Gnshm12_resampled.0.05.grd -R110/156/-46/-9 -I0.05')