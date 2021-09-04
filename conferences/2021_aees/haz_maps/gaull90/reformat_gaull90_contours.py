from gmt_tools import shp2xyz
from os import system

shpfile = 'Gaull_1990_PGA_WGS84.shp'

shp2xyz(shpfile, 'gaull90_contours.xyz', 'PGA_value')

# make grd file
print('Making surface...')
system('gmt5 surface gaull90_contours.xyz -Ggaull90_interp.0.05.grd -R110/156/-46/-9 -I0.05')