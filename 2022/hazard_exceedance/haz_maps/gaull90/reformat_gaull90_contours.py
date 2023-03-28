from gmt_tools import shp2xyz
from os import system

shpfile = 'Gaull_1990_PGA_WGS84_appended.shp'

shp2xyz(shpfile, 'gaull90_contours.xyz', 'PGA_value2')

# make grd file
print('Making surface...')
system('gmt5 surface gaull90_contours.xyz -Ggaull90_interp.0.05.grd -R110/156/-46/-9 -I0.05 -T0.15')

# convert from m/s**2 to g
system('gmt5 grdmath gaull90_interp.0.05.grd 9.80665 DIV = gaull90_interp.0.05.grd')