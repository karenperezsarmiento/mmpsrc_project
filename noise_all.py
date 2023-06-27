import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import re

reduc = "_2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1_files.txt"
direc = "ACT_Sources_2023_0f09-to-35f5_PCA0/"
cluster_name = []
good_noise_area = []
b = np.linspace(0.0,0.01,101)
H_total = np.zeros(len(b)-1)
A_total = np.zeros(len(b)-1)
min_rad = [0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5]
max_rad = [0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0]
mask = np.zeros((len(min_rad),len(b)-1))


cluster_list = np.array(pd.read_csv("/users/ksarmien/mmpsrc_project/reductions_lists/"+reduc,header=None))
median_noise = np.zeros((len(min_rad),len(cluster_list)))
ii=0
for i in cluster_list:
	cluster = i[0][2:20]
	c_snr =  "/home/scratch/cromero/MUSTANG2/Reductions/"+direc+i[0][2:]
	hdu_noise = fits.open(c_snr)[1]
	img_noise = hdu_noise.data
	arcsec_per_pixel = (3600*np.abs(hdu_noise.header["CD1_1"]))
	pixel_per_arcsec = 1./arcsec_per_pixel
	pixel_per_arcmin = pixel_per_arcsec*60.
	arcsec_per_pixel_sq = arcsec_per_pixel**2
	x = img_noise.shape[1]
	y = img_noise.shape[0]
	x_arr = np.arange(np.floor(-1*x/2),np.floor(x/2))
	y_arr = np.arange(np.floor(-1*y/2),np.floor(y/2))
	x_mesh,y_mesh = np.meshgrid(x_arr,y_arr)
	r = np.sqrt(x_mesh**2+y_mesh**2)
	for j in range(len(min_rad)):
		n = img_noise[np.logical_and((r<=max_rad[j]*pixel_per_arcmin),(r>min_rad[j]*pixel_per_arcmin))]
		mask[j] = mask[j]+np.cumsum(np.histogram(n,bins=b)[0])
		med = np.median(n)
		median_noise[j,ii] = med
	H,x_H =np.histogram(img_noise.ravel(),bins=b)
	H_c = np.cumsum(H)
	H_total = H_total+H_c
	arcsec_per_pixel_sq = (3600*np.abs(hdu_noise.header["CD1_1"]))**2
	A_c = H_c*arcsec_per_pixel_sq
	A_total = A_total+A_c
	ii+=1

index_labels = list(np.array(list(zip(min_rad,max_rad))))
index_labels = map(str,index_labels)
df = pd.DataFrame(mask,columns = list(b[1:]),index=index_labels)
df.to_csv("/users/ksarmien/mmpsrc_project/map_quality_tables/noise_by_radii_hist.csv")

fig = plt.figure()
plt.plot(b[1:],H_total)
plt.xlabel("Noise level")
plt.title("Cumulative count of pixels above a noise level")
plt.xlim(0.0,0.003)
plt.savefig("noise_histogram/cumulative_noise_pixels.png")
plt.close(fig)


fig = plt.figure()
for i in range(len(min_rad)-4):
	plt.plot(b[1:],mask[i],label="r<="+str(max_rad[i])+"  r>"+str(min_rad[i]))
plt.xlabel("Noise level")
plt.legend()
plt.title("Cumulative histogram of pixels above a noise level in given annuli")
plt.xlim(0.0,0.003)
plt.savefig("noise_histogram/cumulative_noise_pixels_radii.png")
plt.close(fig)

fig = plt.figure()
for i in range(len(min_rad)-4):
	hist,bins = np.histogram(median_noise[i],bins=np.linspace(0,0.003,100))
	plt.step(bins[:-1],hist,where="post",label="r<="+str(max_rad[i])+"  r>"+str(min_rad[i]))
plt.title("Histogram of median noise per annuli across all cluster maps")
plt.legend()
plt.xlim(0.0,0.002)
plt.savefig("noise_histogram/hist_median_noise.png")
plt.close(fig)

for i in range(len(min_rad)-4):
	fig = plt.figure()
	hist,bins = np.histogram(median_noise[i],bins=np.linspace(0,0.003,100))
	plt.step(bins[:-1],hist,where="post",label="r<="+str(max_rad[i])+"  r>"+str(min_rad[i]))
	plt.title("Histogram of median noise in annuli across all cluster maps")
	plt.yscale("log")
	plt.legend()
	if i<5:
		plt.xlim(0.0,0.0005)
	else:
		plt.xlim(0.0,0.001)
	plt.savefig("noise_histogram/hist_median_noise"+str(max_rad[i])+".png")
	plt.close(fig)
