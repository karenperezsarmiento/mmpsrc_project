import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import re

reduc_list = ["_2aspcmsubqm2_fitel_0f11-to-25f1Hz_qc_1p2rr_M_PdoCals_dt20_snr_iter1_files.txt","_2asp_pca3_qm2_fitel_0f11-to-25f5Hz_qc_1p2rr_M_PdoCals_dt20_snr_iter1_files.txt"]

cluster_name = []
good_noise_area = []

for k in range(len(reduc_list)):
	cluster_list = np.array(pd.read_csv("/users/ksarmien/mmpsrc_project/reductions_lists/"+reduc_list[k]))
	for i in cluster_list:
		cluster = i[0][2:20]
		c_snr =  "/home/scratch/cromero/MUSTANG2/Reductions/"+cluster+i[0][20:]
		c_map = re.sub('_snr_iter1.fits','_map_iter1.fits',c_snr)
		hdu_noise = fits.open(c_map)[1]
		img_noise = hdu_noise.data
		fig = plt.figure()
		plt.hist(img_noise)
		plt.savefig("noise_histogram/"+cluster+"reduc_"+str(k)+".png")
		plt.close(fig)

		 
