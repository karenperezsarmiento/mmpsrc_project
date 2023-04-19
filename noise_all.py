import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import re

reduc_list = ["_2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1_files.txt"]
direc_list= ["ACT_Sources_2023_0f09-to-35f5_PCA0/"]
cluster_name = []
good_noise_area = []
x = np.linspace(-0.1,0.1,2001)
H_total = np.zeros(len(x)-1)
A_total = np.zeros(len(x)-1)


for k in range(len(reduc_list)):
	cluster_list = np.array(pd.read_csv("/users/ksarmien/mmpsrc_project/reductions_lists/"+reduc_list[k],header=None))
	for i in cluster_list:
		cluster = i[0][2:20]
		c_snr =  "/home/scratch/cromero/MUSTANG2/Reductions/"+direc_list[k]+i[0][2:]
		c_noise = re.sub('snr','noise',c_snr)
		hdu_noise = fits.open(c_noise)[0]
		img_noise = hdu_noise.data
		H,x_H =np.histogram(img_noise.ravel(),bins=x)
		H_c = np.cumsum(H)
		H_total = H_total+H_c
		arcsec_per_pixel_sq = (3600*np.abs(hdu_noise.header["CD1_1"]))**2
		A_c = H_c*arcsec_per_pixel_sq
		A_total = A_total+A_c


df = np.array(list(zip(x[1:],H_total,A_total)))
df = pd.DataFrame(df,columns = ["noise","cumulative_pix_num","cumulative_area_sqarcsec"])
df.to_csv("cumulative_noise_hist.csv")

fig = plt.figure()
plt.plot(x[1:],H_total)
plt.xlabel("Noise level")
plt.title("Cumulative count of pixels above a noise level")
plt.xlim(-0.02,0.02)
plt.savefig("cumulative_noise_pixels.png")
plt.close(fig)

fig = plt.figure()
plt.plot(x[1:],A_total)
plt.xlabel("Noise level")
plt.title("Cumulative area (arcsec sq) above a noise level")
plt.xlim(-0.02,0.02)
plt.savefig("cumulative_noise_area.png")
plt.close(fig)
