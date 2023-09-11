import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import re
import pandas as pd

reduc_list = ["injected_filter.txt"]
direc_list = ["/home/scratch/sdicker/AGBT21B_298/"]

for k in range(len(reduc_list)):
	cluster_list = np.array(pd.read_csv("../reductions_lists/"+reduc_list[k],header=None))
	for i in cluster_list:
		cluster = i[0][4:22]
		session = i[0][97:99]
		c_snr =  direc_list[k]+i[0][1:]
		print(cluster)
		print(session)
		hdu_snr = fits.open(c_snr)[0]
		img_snr = hdu_snr.data
		fig = plt.figure()
		plt.imshow(img_snr,cmap="bwr")
		plt.title(cluster)
		plt.savefig("../snr_filter_inj/"+cluster+"_"+session+".png")
		plt.close(fig)
