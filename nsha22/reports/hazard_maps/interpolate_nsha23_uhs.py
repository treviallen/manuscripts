#def get_nsha12_hazard_spectra(lon, lat, probability_50, place):
'''
lon = longitude (decimal degrees)

lat = latitude (decimal degrees)

probability_50: list of strings, i.e. 10, 3.3, 2

place = free text

alias gmt='/opt/local/bin/gmt5'

'''
# test data
place = 'test site'
lon = 134
lat = -30
probability_50 = 10

from os import path, system
from numpy import array

# provide warning if lon is negative
if lon < 0:
    print('!!! CHECK LAT LON ORDER !!!')
    
# provide warning if lat is positive
if lat > 0:
    print('!!! CHECK THAT LAT IS NEGATIVE VALUE !!!')

# set grid spectral periods for hazard curve
spectral_periods_str = ['PGA', '01', '015', '02', '03', '04', '05', '07', \
                        '10', '15', '20', '30']

spectral_periods = [0.0, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0]

# rename place
place = place.replace(' ','_')
probability_50_decimal = probability_50/100.

# write gmt lon/lat file
f = open('lola.txt', 'w')
f.write('\t'.join((str(lon), str(lat))))
f.close()

# set base path for grd location on the NAS
basepath = 'netcdf'

uhs = []
#periods = []

uhstxt = "Uniform Hazard Spectra based on Geoscience Australia's 2023 National Seismic Hazard Assessment\n"
uhstxt += 'Spectral acceleration for an exceedance probability of '+str(probability_50)+'% in 50 years in units of g (where PGA=0.0s)\n'
uhstxt += 'PERIOD,ACCEL\n'
for spectral_period_str, spectral_period in zip(spectral_periods_str, spectral_periods):
    # NSHA23.PGA.0.1.B.grd
    if spectral_period_str.startswith('PGA'):
        grdfile = '.'.join(('NSHA23', spectral_period_str, str(probability_50_decimal), 'B', 'grd'))
    else:
        grdfile = '.'.join(('NSHA23', 'SA'+spectral_period_str, str(probability_50_decimal),'B', 'grd'))
    grdpath = path.join(basepath, grdfile)

    try:
        # do grdtrack to extract hazard value
        print(''.join(('gmt5 grdtrack lola.txt -G', grdpath, ' > lolahaz.txt')))
        system(''.join(('gmt5 grdtrack lola.txt -G', grdpath, ' > lolahaz.txt')))
    
    
        # parse in hazard value
        hazval = float(open('lolahaz.txt').read().strip().split('\t')[-1])
        
        # append to hazArray
        uhs.append(hazval)
        #periods.append(float(spectral_period))
        
        # add text to output
        uhstxt += ','.join((spectral_period_str, str('%0.4f' % hazval))) + '\n'
        
    except:
        print('File not found:', grdpath)
        
# write to file

uhsFile = '_'.join(('NSHA23_UHS',place,str(probability_50_decimal),str(lon),str(lat)))+'.csv'
f = open(uhsFile, 'w')
f.write(uhstxt)
f.close()    
    
#return array(periods), array(uhs)