import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
from skimage import data, feature

def blob_dog_flatten(map_arr,theta1,theta2,th):
	b = feature.blob_dog(image=map_arr,min_sigma=theta1,max_sigma=theta2,threshold=th)
	return np.split(b,b.shape[0],axis=0)

def dog_finder(map_arr,theta1,theta2,th):
	blob_func_v = np.vectorize(blob_dog_flatten,excluded=["map_arr","th"],signature='(k,h),(),(),()->(m,n)')
	#blob_func_v = np.vectorize(blob_dog_flatten,excluded=["image","sigma_ratio","threshold","overlap","threshold_rel","exclude_border"],signature='(),()->(m,n)')
	blobs = blob_func_v(map_arr,theta1,theta2,th)
	return blobs

dir_map = "/home/scratch/sdicker/AGBT21B_298"
reduc = "_2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1"
reduction = "injected.txt"
reduction_list = "/users/ksarmien/mmpsrc_project/reductions_lists/"+reduction
nsigma = 5.

file_list = list(pd.read_csv(reduction_list))
file_name = dir_map+file_list[0][1:]

hdu = fits.open(file_name)
snr_map = hdu[0].data

print(snr_map.shape)
theta1 = [2]
theta2 = [6]

th = np.std(snr_map)*nsigma

res = dog_finder(snr_map,theta1,theta2,th)

print(blob_dog_flatten(snr_map,2,6,th))
print(res)
print(res.shape)

