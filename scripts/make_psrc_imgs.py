import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u

reduc_list = "_2aspcmsubqm2_fitel_0f05-to-49f5Hz_qc_0p6rr_M_PdoCals_dt10_snr"
direc_list = "real_maps/"
code = "21"
theta_1 = 2.0
theta_2 = 6.0
all_vlass = pd.read_csv("../catalogs/all_vlass_sources_in_cluster_maps_rad_6arcmin.csv")
all_first = pd.read_csv("../catalogs/all_first_sources_in_cluster_maps_rad_6arcmin.csv")

psrcs_found = pd.read_csv("../psrc_lists/final_cat_snr_5_3arcmin_snr_7.csv")
matched = pd.read_csv("../psrc_lists/final_matched_snr_5_3arcmin_snr_7.csv")

clusters = np.array(list(psrcs_found["cluster"]))
for i in range(len(psrcs_found)):
	ci = np.array(psrcs_found["cluster"])[i]
	xi = int(np.array(psrcs_found["x"])[i])
	yi = int(np.array(psrcs_found["y"])[i])
	c_snr = "/home/scratch/sdicker/AGBT21B_298/"+direc_list+"Jy_"+ci+reduc_list+".fits"
	hdu_snr = fits.open(c_snr)[0]
	wcs = WCS(hdu_snr.header)
	fig = plt.figure()
	img_snr = hdu_snr.data
	plt.subplot(projection=wcs)
	plt.imshow(img_snr[yi-29:yi+29,xi-29:xi+29],origin="lower",cmap="bwr")
	plt.title(str(i))	
	#plt.scatter(x_first,y_first,color="none",marker="s",edgecolor="blue")
	plt.scatter(30,30,color="none",edgecolor="green",s=500)	
	plt.savefig("../psrc_img/"+str(i)+"_psrc_first.png")
	print(str(i))
	plt.close(fig)
