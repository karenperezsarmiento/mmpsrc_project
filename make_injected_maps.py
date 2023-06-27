import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import re
import pandas as pd

reduc_list = ["injected.txt"]
direc_list = ["/home/scratch/sdicker/AGBT21B_298/"]
code = ["20"]

for k in range(len(reduc_list)):
	cluster_list = np.array(pd.read_csv("/users/ksarmien/mmpsrc_project/reductions_lists/"+reduc_list[k],header=None))
	for i in cluster_list:
		cluster = i[0][2:20]
		c_snr =  direc_list[k]+i[0][2:]
		print(c_snr)
		c_map = re.sub("snr","map",c_snr)
		hdu_snr = fits.open(c_snr)[0]
		hdu_map = fits.open(c_map)[0]
		img_snr = hdu_snr.data
		img_map = hdu_map.data
		fig = plt.figure()
		plt.imshow(img_snr,cmap="bwr")
		plt.title(cluster)
		plt.savefig("snr_maps_injected/"+cluster+"_inj_reduc_"+code[k]+".png")
		plt.close(fig)
		fig = plt.figure()
		plt.imshow(img_map,cmap="bwr")
		plt.title(cluster)
		plt.savefig("maps_injected/"+cluster+"_inj_reduc_"+code[k]+".png")
		plt.close(fig)
