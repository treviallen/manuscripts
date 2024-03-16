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
nc1 = dictlist2array(coeffs, 'nc1s')

fig = plt.figure(1, figsize=(12, 9))

plt.subplot(411)
plt.semilogx(freqs, nc0, 'ro')
plt.semilogx(freqs, nc0s, 'bo')
plt.ylim([-2,0])
plt.ylabel('nc0')

plt.subplot(412)
plt.semilogx(freqs, nc1, 'ro')
plt.ylabel('nc1')

# plt mc0
mc0 = dictlist2array(coeffs, 'mc0')
mc0s = dictlist2array(coeffs, 'mc0s')
plt.subplot(413)
plt.semilogx(freqs, mc0, 'ro')
plt.semilogx(freqs, mc0s, 'bo')
print(median(mc0))
plt.ylabel('mc0')
#nc0s = dictlist2array(coeffs, 'nc0s')

# plt mc1
mc1 = dictlist2array(coeffs, 'mc1')
plt.subplot(414)
plt.semilogx(freqs, mc1, 'ro')
plt.ylabel('mc1')

fig = plt.figure(1, figsize=(12, 10))

'''
#plt.semilogx(freqs, nc0s, 'bo')

# plt mc0
mc1 = dictlist2array(coeffs, 'mc1')
#nc0s = dictlist2array(coeffs, 'nc0s')


fig = plt.figure(1, figsize=(12, 9))

plt.subplot(313)
plt.semilogx(freqs, mc1, 'ro')
idx = where((freqs > 0.5) & (freqs < 10))[0]
print(median(mc1[idx]))
'''
plt.show()
