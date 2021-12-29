import shapefile
from mapping_tools import get_field_data, distance, reckon, get_line_parallels
from fault_tools import mag2wid_L10, mag2rupwid_WC94
from numpy import radians, cos, sin

shpfile = 'simple_aus_ruptures.shp'

print 'Reading source shapefile...'
sf = shapefile.Reader(shpfile)

ids = get_field_data(sf, 'ID', 'str')
mw = get_field_data(sf, 'MW', 'float')
dip = get_field_data(sf, 'DIP', 'float')


shapes = sf.shapes()

for i, shape in enumerate(shapes):
	
	#rupwid = mag2wid_L10(mw[i], 'scr')
	rupwid = mag2rupwid_WC94(mw[i], 'rs')
	
	rngkm = rupwid * cos(radians(dip[i]))
	rupdep = rupwid * sin(radians(dip[i]))
	
	
	posazpts, negazpts = get_line_parallels(shape.points, rngkm)
	
	ftxt = '\t'.join((str('%0.3f' % shape.points[0][0]),str('%0.3f' % shape.points[0][1]), '0.0')) + '\n'
	ftxt += '\t'.join((str('%0.3f' % shape.points[1][0]),str('%0.3f' % shape.points[1][1]), '0.0')) + '\n\n'
	
	ftxt += '\t'.join((str('%0.3f' % posazpts[1][0]),str('%0.3f' % posazpts[1][1]), str('%0.1f' % rupdep))) + '\n'
	ftxt += '\t'.join((str('%0.3f' % posazpts[1][0]),str('%0.3f' % posazpts[0][1]), str('%0.1f' % rupdep))) + '\n\n'
	
	ftxt += '\t'.join((str('%0.3f' % negazpts[1][0]),str('%0.3f' % negazpts[1][1]), str('%0.1f' % rupdep))) + '\n'   
	ftxt += '\t'.join((str('%0.3f' % negazpts[0][0]),str('%0.3f' % negazpts[0][1]), str('%0.1f' % rupdep))) + '\n\n'
	
	ftxt += '\t'.join((str('%0.3f' % shape.points[0][0]),str('%0.3f' % shape.points[0][1]), '0.0'))
		
	f = open(str(ids[i])+'_fault.txt', 'wb')
	f.write(ftxt)
	f.close()
