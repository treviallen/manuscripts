import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.style.use('classic')


print('Loading pkl file...')
stdict = pickle.load(open("stdict.pkl", "rb" ))
