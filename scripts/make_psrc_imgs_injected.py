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
fake = pd.read_csv("../psrc_lists/fake_pnt_amplitudes.csv")
injected_found = pd.read_csv("../psrc_lists/injected_2_6_psrcs__2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1_4.0sigma.csv")

ff = np.unique(np.array(pd.read_csv("../reductions_lists/injected.txt",header=None)))
for i in ff:
	c_snr = direc_list+i[1:]
	hdu_snr = fits.open(c_snr)[0]
	wcs = WCS(hdu_snr.header)
	c = i[4:22]
	session = i[88:-15]
	print(c)
	print(session)
	fake_cluster = fake.loc[fake["cluster"]==c]
	found_cluster = injected_found.loc[injected_found["cluster"]==c]
	fig = plt.figure()
	img_snr = hdu_snr.data
	plt.subplot(projection=wcs)
	plt.imshow(img_snr,origin="lower",cmap="bwr")
	plt.title(str(c))	
	#plt.scatter(x_first,y_first,color="none",marker="s",edgecolor="blue")
	plt.scatter(np.array(fake_cluster["x_loc"]),np.array(fake_cluster["y_loc"]),color="none",edgecolor="green")	
	plt.scatter(np.array(found_cluster["x"]),np.array(found_cluster["y"]),color="none",edgecolor="black")    
	plt.savefig("../psrc_img_injected/"+str(c)+"_"+str(session)+"_injected_found.png")
	plt.close(fig)
