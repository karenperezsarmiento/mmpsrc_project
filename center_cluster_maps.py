import numpy as np
import pandas as pd
from astropy.io import fits
from astropy import wcs

reduc = "_2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1"
reduction = reduc+"_files.txt"
reduction_list = "/users/ksarmien/mmpsrc_project/reductions_lists/"+reduction

cluster_list = np.array(pd.read_csv(reduction_list,header=None))
cluster_name = []
cluster_ra = []
cluster_dec = []

for i in cluster_list:
	cluster = i[0][2:20]
	c_snr = "/home/scratch/cromero/MUSTANG2/Reductions/ACT_Sources_2023_0f09-to-35f5_PCA0/"+cluster+i[0][20:]
	hdu_snr = fits.open(c_snr)[0]
	world = wcs.WCS(hdu_snr.header)	
	crpixx =hdu_snr.header["CRPIX1"] #central PIX IN X
	crpixy =hdu_snr.header["CRPIX2"] #central PIX IN Y
	crval_ra = hdu_snr.header["CRVAL1"]
	crval_dec = hdu_snr.header["CRVAL2"]
	cluster_name = np.append(cluster_name,cluster)
	cluster_ra = np.append(cluster_ra,crval_ra)
	cluster_dec = np.append(cluster_dec,crval_dec)
	
central_arr = np.column_stack((cluster_name,cluster_ra,cluster_dec))
df_central = pd.DataFrame(central_arr, columns = ["cluster","ra","dec"])
df_central.to_csv("/users/ksarmien/mmpsrc_project/catalogs/cluster_map_centers.csv")
