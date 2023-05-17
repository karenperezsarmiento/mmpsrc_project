import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import pandas as pd
import re
import scipy.optimize as opt

reduc_list = ["_2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1_files.txt"]
direc_list = ["ACT_Sources_2023_0f09-to-35f5_PCA0/"]
code = ["20"]
data_qual = pd.read_csv("data_quality_code_20.csv")

cluster_name = []
good_noise_area = []

b = np.linspace(-7.0,7.0,100)
min_rad = [0.0,1.0,2.0,3.0,4.0,5.0,6.0]
max_rad = [1.0,2.0,3.0,4.0,5.0,6.0,7.0]
H = np.zeros((len(min_rad),len(b)-1))
N_pix = np.zeros(len(min_rad))


def gaussian(xdata,amp,mean,std):
        return amp*np.exp(-1*(xdata-mean)**2/(2*std**2))

def minimize(pars,xdata,hist_data):
        res = hist_data - gaussian(xdata, *pars)
        return np.sqrt(np.mean(res**2))

for k in range(len(reduc_list)):
	cluster_list = np.array(pd.read_csv("/users/ksarmien/mmpsrc_project/reductions_lists/"+reduc_list[k],header=None))
	all_means = np.zeros((len(min_rad),len(cluster_list)))
	all_std = np.zeros((len(min_rad),len(cluster_list)))
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
		histo = np.histogram(img[img_weight>0.0].ravel(),bins = b)[0]
		x_histo = b[:-1] + np.diff(b)/2
		fit = gaussian(x_histo,amp,mean,std)
		arcsec_per_pixel = (3600*np.abs(hdu.header["CD1_1"]))
		pixel_per_arcsec = 1./arcsec_per_pixel
		pixel_per_arcmin = pixel_per_arcsec*60.
		arcsec_per_pixel_sq = arcsec_per_pixel**2
		x = img.shape[1]
		y = img.shape[0]
		x_arr = np.arange(np.floor(-1*x/2),np.floor(x/2))
		y_arr = np.arange(np.floor(-1*y/2),np.floor(y/2))
		x_mesh,y_mesh = np.meshgrid(x_arr,y_arr)
		r = np.sqrt(x_mesh**2+y_mesh**2)
		for j in range(len(min_rad)):
			n = img[np.logical_and((r<=max_rad[j]*pixel_per_arcmin),(r>min_rad[j]*pixel_per_arcmin))]
			H[j] = H[j] + np.histogram(n,bins=b)[0]
			N_pix[j] = N_pix[j] + len(n)
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


np.savetxt("histo_data.csv",H) 
fig = plt.figure()
for i in range(len(min_rad)-3):
	p = opt.minimize(minimize,[np.max(H[i]),0,1],args=(b[1:],H[i])).x
	plt.step(b[:-1],H[i]/N_pix[i],label="r<="+str(max_rad[i])+"  r>"+str(min_rad[i]))
plt.title("Histogram of SNR per annuli across all cluster maps")
plt.legend()
plt.xlim(-7,7)
plt.savefig("hist_snr_all.png")
plt.close(fig)
