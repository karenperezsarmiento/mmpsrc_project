import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u

reduc_list = "injected.txt"
direc_list = "/home/scratch/sdicker/AGBT21B_298/"
code = "20"
theta_1 = 2.0
theta_2 = 6.0
fake = pd.read_csv("/users/ksarmien/mmpsrc_project/psrc_lists/fake_pnt_amplitudes.txt",names=["Cluster","amp(jy)","X_loc(pixels)","Y_loc(pixels)"],header=None,delimiter=r"\s+")
injected_found = pd.read_csv("/users/ksarmien/mmpsrc_project/psrc_lists/injected_psrcs__2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1_4.0sigma.csv")

ff = np.unique(np.array(pd.read_csv("/users/ksarmien/mmpsrc_project/reductions_lists/injected.txt",header=None)))
for i in ff:
	c_snr = direc_list+i[2:]
	hdu_snr = fits.open(c_snr)[0]
	wcs = WCS(hdu_snr.header)
	c = i[2:20]
	fake_cluster = fake.loc[fake["Cluster"]==c]
	found_cluster = injected_found.loc[injected_found["cluster"]==c]
	fig = plt.figure()
	img_snr = hdu_snr.data
	plt.subplot(projection=wcs)
	plt.imshow(img_snr,origin="lower",cmap="bwr")
	plt.title(str(i))	
	#plt.scatter(x_first,y_first,color="none",marker="s",edgecolor="blue")
	plt.scatter(np.array(fake_cluster["X_loc(pixels)"]),np.array(fake_cluster["Y_loc(pixels)"]),color="none",edgecolor="green")	
	plt.scatter(np.array(found_cluster["x"]),np.array(found_cluster["y"]),color="none",edgecolor="black")    
	plt.savefig("psrc_img_injected/"+str(i)+"_injected_found.png")
	plt.close(fig)
