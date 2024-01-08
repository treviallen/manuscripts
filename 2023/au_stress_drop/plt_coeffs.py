import pickle
import matplotlib.pyplot as plt
from misc_tools import dictlist2array
from numpy import mean, median, where


# load atten coeffs
coeffs = pickle.load(open('atten_coeffs.pkl', 'rb' ))

freqs = dictlist2array(coeffs, 'freq')

# plt nc0
nc0 = dictlist2array(coeffs, 'nc0')
nc0s = dictlist2array(coeffs, 'nc0s')

fig = plt.figure(1, figsize=(12, 9))

plt.subplot(311)
plt.semilogx(freqs, nc0, 'ro')
plt.semilogx(freqs, nc0s, 'bo')
plt.ylim([-2,0])

# plt mc0
mc0 = dictlist2array(coeffs, 'mc0')
print(median(mc0))
#nc0s = dictlist2array(coeffs, 'nc0s')

fig = plt.figure(1, figsize=(12, 9))

plt.subplot(312)
plt.semilogx(freqs, mc0, 'ro')
#plt.semilogx(freqs, nc0s, 'bo')

# plt mc0
mc1 = dictlist2array(coeffs, 'mc1')
#nc0s = dictlist2array(coeffs, 'nc0s')


fig = plt.figure(1, figsize=(12, 9))

plt.subplot(313)
plt.semilogx(freqs, mc1, 'ro')
idx = where((freqs > 0.5) & (freqs < 10))[0]
print(median(mc1[idx]))

plt.show()
