import pickle
import matplotlib.pyplot as plt
from misc_tools import dictlist2array
from numpy import mean, median, where, array, log10


# load atten coeffs
coeffs = pickle.load(open('atten_coeffs.pkl', 'rb' ))

freqs = dictlist2array(coeffs, 'freq')

# plt nc0
nc0 = dictlist2array(coeffs, 'nc0')
nc0s = dictlist2array(coeffs, 'nc0s')

mc0f = dictlist2array(coeffs, 'mc0f')
mc1 = dictlist2array(coeffs, 'mc1')
mc1s = dictlist2array(coeffs, 'mc1s')
mc1f = dictlist2array(coeffs, 'mc1f')

magc0s = dictlist2array(coeffs, 'magc0')
magc1s = dictlist2array(coeffs, 'magc1')

#mc1ts = dictlist2array(coeffs, 'mc1fs')

##################################################################
txt = 'freq,c0,c1,c2,m0,m1\n'

for i in range(0, len(freqs)):
	txt += ','.join((str('%0.4f' % freqs[i]), str('%0.4f' % nc0[i]), str('%0.4f' % mc0f[i]), str('%0.4e' % mc1f[i]), \
	                 str('%0.4f' % magc0s[i]), str('%0.4f' % magc1s[i]))) + '\n'
	                
f = open('swan_fas_coeffs.csv', 'w')
f.write(txt)
f.close()
