import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import re

reduc_list = ["_2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1_files.txt"]
direc_list = ["ACT_Sources_2023_0f09-to-35f5_PCA0/"]
code = ["20"]
data_qual = pd.read_csv("data_quality_code_20.csv")

cluster_name = []
good_noise_area = []

def gaussian(xdata,amp,mean,std):
        return amp*np.exp(-1*(xdata-mean)**2/(2*std**2))

for k in range(len(reduc_list)):
	cluster_list = np.array(pd.read_csv("/users/ksarmien/mmpsrc_project/reductions_lists/"+reduc_list[k]))
	for i in cluster_list:
		cluster = i[0][2:20]
		filtered = data_qual.loc[data_qual["Source"]==cluster]
		amp = filtered["amp_fit_snr"].values[0]
		mean = filtered["mean_fit_snr"].values[0]
		std = filtered["std_fit_snr"].values[0]
		c_snr =  "/home/scratch/cromero/MUSTANG2/Reductions/"+direc_list[k]+i[0][2:]
		c_weight = re.sub("snr","map",c_snr)
		print(c_snr)
		hdu = fits.open(c_snr)[0]
		hdu_weight = fits.open(c_weight)[1]
		img = hdu.data
		img_weight = hdu_weight.data
		histo = np.histogram(img[img_weight>0.0].ravel(),bins = np.linspace(-7.0,7.0,100))[0]
		x_histo = np.linspace(-7.0,7.0,100)[:-1] + np.diff(np.linspace(-7.0,7.0,100))/2
		fit = gaussian(x_histo,amp,mean,std)
		textstr = '\n'.join((
    			r'$\mu=%.2f$' % (mean, ),
    			r'$\sigma=%.2f$' % (std, )))
		fig,ax = plt.subplots()
		ax.hist(img[img_weight>0.0].ravel(),bins=np.linspace(-7.0,7.0,100))
		ax.plot(x_histo,fit,c="r")
		ax.set_xlim((-7.0,7.0))
		props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
		ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,verticalalignment='top', bbox=props)
		ax.set_title(cluster)
		plt.savefig("snr_histogram/"+cluster+"reduc_"+code[k]+".png")
		plt.close(fig)

		 
