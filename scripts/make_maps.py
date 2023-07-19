import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import re
import pandas as pd

reduc_list = ["_2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1_files.txt"]
direc_list = ["ACT_Sources_2023_0f09-to-35f5_PCA0/"]
code = ["20"]

for k in range(len(reduc_list)):
	cluster_list = np.array(pd.read_csv("../reductions_lists/"+reduc_list[k],header=None))
	for i in cluster_list:
		cluster = i[0][2:20]
		c_snr =  "/home/scratch/cromero/MUSTANG2/Reductions/"+direc_list[k]+i[0][2:]
		print(c_snr)
		c_map = re.sub("snr","map",c_snr)
		hdu_snr = fits.open(c_snr)[0]
		hdu_map = fits.open(c_map)[0]
		img_snr = hdu_snr.data
		img_map = hdu_map.data
		fig = plt.figure()
		plt.imshow(img_snr,cmap="bwr")
		plt.title(cluster)
		plt.savefig("../snr_maps/"+cluster+"reduc_"+code[k]+".png")
		plt.close(fig)
		fig = plt.figure()
		plt.imshow(img_map,cmap="bwr")
		plt.title(cluster)
		plt.savefig("../maps/"+cluster+"reduc_"+code[k]+".png")
		plt.close(fig)
