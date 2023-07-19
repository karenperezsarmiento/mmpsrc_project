import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u

reduc_list = "_2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1"
direc_list = "ACT_Sources_2023_0f09-to-35f5_PCA0/"
code = "20"
theta_1 = 2.0
theta_2 = 6.0
all_vlass = pd.read_csv("all_vlass_sources_in_cluster_maps_rad_6arcmin.csv")
all_first = pd.read_csv("all_first_sources_in_cluster_maps_rad_6arcmin.csv")

psrcs_found = pd.read_csv("psrc_lists/final_cat_snr_5_3arcmin_snr_7.csv")
matched = pd.read_csv("psrc_lists/final_matched_snr_5_3arcmin_snr_7.csv")

clusters = np.unique(np.array(list(psrcs_found["cluster"])))
for i in clusters:
	c_snr = "/home/scratch/cromero/MUSTANG2/Reductions/"+direc_list+i+"/Jy_"+i+reduc_list+".fits"
	hdu_snr = fits.open(c_snr)[0]
	wcs = WCS(hdu_snr.header)
	found_cluster = psrcs_found.loc[psrcs_found["cluster"]==i]
	first_cluster = all_first.loc[all_first["cluster"]==i]
	sky_first = SkyCoord(ra=np.array(first_cluster["RA_2"])*u.degree,dec=np.array(first_cluster["DEC_2"])*u.degree,frame="icrs")
	x_first,y_first = wcs.world_to_pixel(sky_first)
	fig = plt.figure()
	img_snr = hdu_snr.data
	plt.subplot(projection=wcs)
	plt.imshow(img_snr,origin="lower",cmap="bwr")
	plt.title(str(i))	
	#plt.scatter(x_first,y_first,color="none",marker="s",edgecolor="blue")
	plt.scatter(np.array(found_cluster["x"]),np.array(found_cluster["y"]),color="none",edgecolor="green")	
	plt.savefig("psrc_img/"+str(i)+"_psrc_first.png")
	plt.close(fig)
