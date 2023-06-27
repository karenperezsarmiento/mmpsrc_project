import numpy as np
import pandas as pd
from astropy.io import fits
import re
from astropy.convolution import Gaussian2DKernel, convolve
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
import scipy.optimize as opt

#reduc_list = ["_2aspcmsubqm2_fitel_0f11-to-25f1Hz_qc_1p2rr_M_PdoCals_dt20_snr_iter1_files.txt","_2asp_pca3_qm2_fitel_0f11-to-25f5Hz_qc_1p2rr_M_PdoCals_dt20_snr_iter1_files.txt"]
red_codes = pd.read_csv("reductions_code.csv")
reduc_list = list(red_codes["reduction"])
codes = np.array(list(red_codes["code"]))
dir_map = "/home/scratch/cromero/MUSTANG2/Reductions/ACT_Sources_2023_0f09-to-35f5_PCA0"


#all_clusters = np.array(pd.read_csv("gbt298_obs.csv")["Source"])
#red1_list = np.array(pd.read_csv("/users/ksarmien/Documents/pnt_src_project/reductions_lists/_2asp_pca3_qm2_fitel_0f11-to-25f5Hz_qc_1p2rr_M_PdoCals_dt20_snr_iter1_files.txt",header=None))
#red2_list = np.array(pd.read_csv("/users/ksarmien/Documents/pnt_src_project/reductions_lists/_2aspcmsubqm2_fitel_0f09-to-41f1Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1_files.txt",header=None))

cluster_name = []
cluster_rms_map = []
cluster_rms_snr = []
cluster_avg_snr = []
cluster_ntail = []
reduction_type = []
gaussian_fit_snr = np.repeat(0,3)

def gaussian(xdata,amp,mean,std):
	return amp*np.exp(-1*(xdata-mean)**2/(2*std**2))

def minimize(pars,xdata,hist_data):
	res = hist_data - gaussian(xdata, *pars)
	return np.sqrt(np.mean(res**2)) 


for k in [20]:
	cluster_list = np.array(pd.read_csv("/users/ksarmien/mmpsrc_project/reductions_lists/"+reduc_list[k],header=None))
	current_code = codes[k]
	reduc = re.sub("_files.txt","",reduc_list[k])
	for i in cluster_list:
		c = i[0]
		substitute1 = dir_map
		substitute2 = "./ACT-CLJ\d+.\d+(\+|-)\d+/Jy_"
		substitute3 = reduc + ".fits"
		substitute4 = "/ACT_Sources_2023_0f09-to-35f5_PCA0"
		key=re.sub(substitute1,"",str(c))
		key=re.sub(substitute2,"",key)
		key=re.sub(substitute3,"",key)
		cluster =re.sub(substitute4,"",key)
		#print(cluster)
		#if current_code<18:
		#	c_snr = "/home/scratch/cromero/MUSTANG2/Reductions/"+cluster+i[0][20:]
		#	c_map = re.sub('_snr_iter1.fits','_map_iter1.fits',c_snr)
		#else:
		#	c_snr = "/home/scratch/sdicker/pnt_source_reductions/"+cluster+i[0][20:]
		#	c_map = re.sub('_snr','_map',c_snr)
		c_snr = dir_map+c[1:]
		c_map = re.sub('_snr_iter1.fits','_map_iter1.fits',c_snr)
		hdu_snr = fits.open(c_snr)[0]
		hdu_map = fits.open(c_map)[0]
		hdu_weight = fits.open(c_map)[1]
		deg_per_pixel = np.abs(hdu_map.header["CD1_1"])
		arcmin_per_pixel = 60* deg_per_pixel	
		arcsec_per_pixel = 3600* deg_per_pixel	
		pix_per_degree = 1/deg_per_pixel
		pix_per_arcmin = 1/arcmin_per_pixel
		pix_per_arcsec = 1/arcsec_per_pixel
		m2_beam_arcsec = 9
		m2_beam_pix = (m2_beam_arcsec/(2*np.sqrt(2*np.log(2))))*pix_per_arcsec
		#kernel= Gaussian2DKernel(x_stddev = m2_beam_pix)
		img_map = hdu_map.data
		img_snr = hdu_snr.data
		img_weight = hdu_weight.data
		histo = np.histogram(img_snr[img_weight>0.0].ravel(),bins=np.linspace(-7.0,7.0,100))[0]
		x_histo = np.linspace(-7.0,7.0,100)[:-1] + np.diff(np.linspace(-7.0,7.0,100))/2
		p = opt.minimize(minimize,[np.max(histo),0,1],args=(x_histo,histo)).x
		gaussian_fit_snr = np.vstack((gaussian_fit_snr,p))
		#res = convolve(img_map,kernel)
		res = gaussian_filter(img_map,m2_beam_pix,order=0)
		x = img_map.shape[1]
		y = img_map.shape[0]
		x_arr = np.arange(np.floor(-1*x/2),np.floor(x/2))
		y_arr = np.arange(np.floor(-1*y/2),np.floor(y/2))
		x_mesh,y_mesh = np.meshgrid(x_arr,y_arr)
		r = np.sqrt(x_mesh**2+y_mesh**2)
		rad = 2*pix_per_arcmin #2 arcmin radius in pixels
		mask = r<rad
		mask = mask.astype(int)
		rms_central_map = np.sqrt(np.sum((mask*res)**2)/np.sum(mask))
		cluster_rms_map = np.append(cluster_rms_map,rms_central_map)
		avg_snr = np.sum(mask*img_snr)/np.sum(mask)
		rms_central_snr = np.sqrt(np.sum((mask*img_snr)**2)/np.sum(mask))
		cluster_rms_snr = np.append(cluster_rms_snr,rms_central_snr)
		cluster_avg_snr = np.append(cluster_avg_snr,avg_snr)
		cluster_name = np.append(cluster_name,cluster)
		snr_flatten = img_snr.flatten()
		bins = np.arange(-7.0,7.0,0.25)
		hist = np.histogram(snr_flatten,bins=bins)
		negative_snr_filter = (hist[1]<=-5.0)
		negative_snr_filter = negative_snr_filter[:-1]
		negative_tail_count = np.sum(hist[0][negative_snr_filter])
		negative_tail_percent = 100* negative_tail_count/np.sum(hist[0])
		cluster_ntail = np.append(cluster_ntail,negative_tail_percent)
		reduction_type = np.append(reduction_type,codes[k])


gaussian_fit_snr = gaussian_fit_snr[1:]
red = np.column_stack((cluster_name,cluster_rms_map,cluster_rms_snr,cluster_avg_snr,cluster_ntail,gaussian_fit_snr,reduction_type))
df_red = pd.DataFrame(red, columns = ['Source','rms_map','rms_snr','avg_snr','ntail','amp_fit_snr','mean_fit_snr','std_fit_snr','red_type'])

df_red.to_csv("/users/ksarmien/mmpsrc_project/map_quality_tables/data_quality_code_20.csv")

