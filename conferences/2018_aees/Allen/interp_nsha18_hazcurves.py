from tools.oq_tools import return_annualised_haz_curves
from numpy import interp, exp, log, array

probs = array([0.02,0.01375,0.01,0.00445,0.0021,0.001,0.0005,0.000404,0.0002,0.0001])

hazcurvefile = '/Users/tallen/Documents/Geoscience_Australia/NSHA2018/source_models/complete_model/final/results_fractiles/hazard_curve-mean-SA(0.2)_1.csv'

# get data from first job
siteDict, imls, investigation_time = return_annualised_haz_curves(hazcurvefile)

interpTXT = 'LON,LAT,' + ','.join(('P'+str(x) for x in probs)) + '\n'

# loop through site dict
for site in siteDict:
    interpHaz = exp(interp(log(probs[::-1]), log(site['poe_probs_annual'][::-1]), log(imls[::-1])))[::-1]
    	
    # make text
    interpTXT += ','.join((str('%0.2f' % site['lon']), str('%0.2f' % site['lat']))) + ',' \
                 + ','.join((str(x) for x in interpHaz)) + '\n'
    
f = open('interp_haz_curves.csv', 'wb')
f.write(interpTXT)
f.close()