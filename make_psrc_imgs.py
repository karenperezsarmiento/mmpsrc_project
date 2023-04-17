import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u

reduc_list = ["_2asp_pca3_qm2_fitel_0f11-to-25f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1","_2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1"]
direc_list = ["ACT_Sources_2023_0f11-to-25f5_PCA3/","ACT_Sources_2023_0f09-to-35f5_PCA0/"]
code = ["19","20"]
theta_1 = 3.0
theta_2 = 5.0
all_vlass = pd.read_csv("all_vlass_sources_in_cluster_maps_rad_6arcmin.csv")
all_first = pd.read_csv("all_first_sources_in_cluster_maps_rad_6arcmin.csv")

for k in range(len(reduc_list)):
	psrc_found = pd.read_csv("psrc_lists/uncleaned_psrcs_sigma_"+reduc_list[k]+"_4.0sigma.csv")
	psrc_found = psrc_found.loc[(psrc_found["theta_1"]==theta_1)&(psrc_found["theta_2"]==theta_2)]
	paths_to_fits = np.array(pd.read_csv("reductions_lists/"+reduc_list[k]+"_files.txt"))
	matched = pd.read_csv("psrc_lists/matched_uncleaned_psrcs_sigma_"+reduc_list[k]+"_4.0sigma.csv")
	matched = matched.loc[(matched["theta_1"]==theta_1)&(matched["theta_2"]==theta_2)]
	for i in paths_to_fits:
		cluster = i[0][2:20]
		psrcs_cluster = psrc_found.loc[psrc_found["cluster"]==cluster]
		matched_cluster = matched.loc[matched["cluster"]==cluster]
		vlass_cluster = all_vlass.loc[all_vlass["cluster"]==cluster]
		first_cluster = all_first.loc[all_first["cluster"]==cluster]
		c_snr = "/home/scratch/cromero/MUSTANG2/Reductions/"+direc_list[k]+i[0][2:]
		print(c_snr)
		#c_map = re.sub("snr","map",c_snr)
		hdu_snr = fits.open(c_snr)[0]
		wcs = WCS(hdu_snr.header)
		sky_vlass = SkyCoord(ra=np.array(vlass_cluster["RA_2"])*u.degree,dec=np.array(vlass_cluster["DEC_2"])*u.degree, frame="icrs")
		sky_first = SkyCoord(ra=np.array(first_cluster["RA_2"])*u.degree,dec=np.array(first_cluster["DEC_2"])*u.degree, frame="icrs")
		x_vlass,y_vlass = wcs.world_to_pixel(sky_vlass)
		x_first,y_first = wcs.world_to_pixel(sky_first)
		fig = plt.figure()
		img_snr = hdu_snr.data
		plt.subplot(projection=wcs)
		plt.imshow(img_snr,origin="lower")
		plt.title(cluster+" red "+code[k])
		plt.scatter(np.array(psrcs_cluster["x"]),np.array(psrcs_cluster["y"]),color="none",edgecolor="black")
		plt.scatter(np.array(matched_cluster["x"]),np.array(matched_cluster["y"]),color="none",edgecolor="red")
		plt.scatter(x_vlass,y_vlass,color="none",marker="s",edgecolor="blue")
		plt.scatter(x_first,y_first,color="none",marker="s",edgecolor="blue")
		plt.savefig("psrc_img/"+cluster+"reduc_"+code[k]+".png")
		plt.close(fig)

