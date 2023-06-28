import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
#from psrc_lib import twoD_Gaussian
from astropy.io import fits
from skimage import data, feature
import scipy.optimize as opt

def twoD_Gaussian(xdata_tuple, amplitude, xo, yo, sigma):
    (x,y) = xdata_tuple
    xo = float(xo)
    yo = float(yo)
    g = amplitude*np.exp( - ((x-xo)**2 + (y-yo)**2)/2*sigma**2)
    return g.ravel()

def minimize(pars,xdata_tuple,map_ravel):
    res = map_ravel - twoD_Gaussian(xdata_tuple, *pars)
    return np.sqrt(np.mean(res**2))

dir_map = "/home/scratch/sdicker/AGBT21B_298"
reduc = "_2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1"
reduction = "injected.txt"
reduction_list = "/users/ksarmien/mmpsrc_project/reductions_lists/"+reduction

nsigma = 4.0


file_list = list(pd.read_csv(reduction_list))
file_name = dir_map+file_list[0][1:]

hdu = fits.open(file_name)
snr_map = hdu[0].data
th= np.std(snr_map)*nsigma

theta1 = 2
theta2 = 6

blobs = feature.blob_dog(snr_map,theta1,theta2,threshold=th)
x_coord = np.array(blobs[:,0],int)
y_coord = np.array(blobs[:,1],int)
sigma_dog = np.array(blobs[:,2],float)

import re

map_file_name = re.sub("snr","map",file_name)
hdu_map = fits.open(map_file_name)
signal_map = hdu_map[0].data

xlen,ylen=signal_map.shape
x = np.linspace(0,xlen-1,xlen)
y = np.linspace(0,ylen-1,ylen)
x,y=np.meshgrid(x,y)
signal_map_ravel = signal_map.ravel()

for src in blobs:
	p = opt.minimize(minimize,[signal_map[int(src[1]),int(src[0])],src[0],src[1],3],args=((x,y),signal_map_ravel)).x

	g_r = twoD_Gaussian((x,y),p[0],p[1],p[2],p[3])
	g_r = g_r.reshape(xlen,ylen)

subtracted = signal_map - g_r


plt.imshow(snr_map)
plt.savefig("test_snr.png")
