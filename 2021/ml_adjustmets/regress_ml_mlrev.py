from scipy.io import loadmat


#load mat file

matfile = '/Users/trev/Documents/Geoscience_Australia/NSHA2018/catalogue/matlab/mdat_ml_rev.mat'

mat_contents = loadmat(matfile)

mdat = mat_contents['mdat']

-4