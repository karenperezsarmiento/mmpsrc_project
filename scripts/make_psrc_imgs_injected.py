import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import re

reduc_list = "injected.txt"
direc_list = "/home/scratch/sdicker/AGBT21B_298/"
code = "20"
theta_1 = 2.0
theta_2 = 6.0
fake = pd.read_csv("../psrc_lists/fake_pnt_amplitudes.csv")
injected_found = pd.read_csv("../psrc_lists/injected_psrcs_2_6_2aspcmsubqm2_fitel_0f05-to-49f5Hz_qc_0p6rr_M_PdoCals_dt10.csv")

ff = np.unique(np.array(pd.read_csv("../reductions_lists/injected.txt",header=None)))
for i in ff:
	c_snr = direc_list+i[2:]
	hdu_snr = fits.open(c_snr)[0]
	wcs = WCS(hdu_snr.header)
	c = i[5:23]
	project = re.search("AGBT\d+B_\d+_\d+",i).group(0)
	scanno = int(re.search("_s\d+_",i).group(0)[2:-1])
	session = i[89:-15]
	print(c)
	print(session)
	fake_cluster = fake.loc[(fake["Cluster"]==c)&(fake["project"]==project)&(fake["scanno"]==scanno)]
	print(project)
	print(scanno)
	print(len(fake_cluster))
	found_cluster = injected_found.loc[(injected_found["cluster"]==c)&(injected_found["project"]==project)&(injected_found["scanno"]==scanno)]
	print(len(found_cluster))
	fig =  plt.figure()
	img_snr = hdu_snr.data
	plt.subplot(projection=wcs)
	plt.imshow(img_snr,origin="lower",cmap="bwr")
	plt.title(str(c)+" project:"+project+" scan #:"+str(scanno))	
	#plt.scatter(x_first,y_first,color="none",marker="s",edgecolor="blue")
	plt.scatter(np.array(fake_cluster["X_loc_pixels"]),np.array(fake_cluster["y_loc_pixels"]),color="none",edgecolor="green")	
	plt.scatter(np.array(found_cluster["x"]),np.array(found_cluster["y"]),color="none",edgecolor="black")    
	plt.savefig("../psrc_img_injected/"+str(c)+"_"+str(session)+"_injected_found.png")
	plt.close(fig)
