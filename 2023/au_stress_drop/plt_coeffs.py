import pickle
import matplotlib.pyplot as plt
from misc_tools import dictlist2array


# load atten coeffs
coeffs = pickle.load(open('atten_coeffs.pkl', 'rb' ))

freqs = dictlist2array(coeffs, 'freq')

# plt c1
nc0 = dictlist2array(coeffs, 'nc0')
nc0s = dictlist2array(coeffs, 'nc0s')

fig = plt.figure(1, figsize=(12, 9))

plt.subplot(311)
plt.semilogx(freqs, nc0, 'ro')
plt.semilogx(freqs, nc0s, 'bo')



plt.show()
