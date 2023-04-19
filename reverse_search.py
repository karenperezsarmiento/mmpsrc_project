import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.io import fits
from astropy import wcs
import re
import scipy.optimize as opt

dir_map = "/home/scratch/cromero/MUSTANG2/Reductions/ACT_Sources_2023_0f09-to-35f5_PCA0"
reduc = "_2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1"
reduc_map = re.sub("snr","map",reduc)
reduction = reduc+"_files.txt"
reduction_list = "/users/ksarmien/mmpsrc_project/reductions_lists/"+reduction

all_vlass = pd.read_csv("/users/ksarmien/mmpsrc_project/all_vlass_sources_in_cluster_maps_rad_6arcmin.csv")
all_first = pd.read_csv("/users/ksarmien/mmpsrc_project/all_first_sources_in_cluster_maps_rad_6arcmin.csv")

def twoD_Gaussian(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    (x,y) = xdata_tuple
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
    return g.ravel()

def minimize(pars,xdata_tuple,map_ravel):
    res = map_ravel - twoD_Gaussian(xdata_tuple, *pars)
    return np.sqrt(np.mean(res**2))

clusters_unique = np.array(pd.read_csv(reduction_list,header=None))
print(len(clusters_unique))
print(clusters_unique.shape)
ra_arr = []
dec_arr =[]
pixx_arr = []
pixy_arr = []
p_map_arr = np.zeros(7)
p_snr_arr = np.zeros(7)
cluster_arr =[]
vlass_name_arr = []
vlass_id_arr = []

for k in clusters_unique:
	cluster = k[0][2:20]
	filtered = all_vlass.loc[all_vlass["cluster"]==cluster]
	print(cluster)
	c_snr = dir_map+"/"+cluster +"/Jy_"+cluster+reduc+".fits"
	c_map = dir_map+"/"+cluster +"/Jy_"+cluster+reduc_map+".fits"
	hdu_snr = fits.open(c_snr)[0]
	img_snr = hdu_snr.data
	hdu_map = fits.open(c_map)[0]
	img_map = hdu_map.data
	world = wcs.WCS(hdu_map.header)
	xlen,ylen=img_map.shape
	x = np.linspace(0,xlen-1,xlen)
	y = np.linspace(0,ylen-1,ylen)
	x,y=np.meshgrid(x,y)
	img_map_ravel = img_map.ravel()
	img_snr_ravel = img_snr.ravel()
	for i in range(len(filtered)):
		cluster = filtered.iloc[i]["cluster"]
		cluster_arr = np.append(cluster_arr,cluster)
		vlass_name_arr = np.append(vlass_name_arr,filtered.iloc[i]["Component_name"])
		vlass_id_arr = np.append(vlass_id_arr,filtered.iloc[i]["Component_id"])
		ra = filtered.iloc[i]["RA_2"]
		dec = filtered.iloc[i]["DEC_2"]
		ra_arr = np.append(ra_arr,ra)
		dec_arr = np.append(dec_arr,dec)
		pixx,pixy = world.wcs_world2pix(ra,dec,1)
		pixx_arr = np.append(pixx_arr,pixx)
		pixy_arr = np.append(pixy_arr,pixy)
		if (pixx<xlen)&(pixy<ylen):
			p = opt.minimize(minimize,[img_map[int(pixy),int(pixx)],pixx,pixy,3,3,0,0],args=((x,y),img_map_ravel)).x
			p_map_arr = np.vstack((p_map_arr,p))
			p_snr = opt.minimize(minimize,[img_snr[int(pixy),int(pixx)],pixx,pixy,3,3,0,0],args=((x,y),img_snr_ravel)).x
			p_snr_arr = np.vstack((p_snr_arr,p_snr))
		else:
			p_map_arr = np.vstack((p_map_arr,np.zeros(7)))
			p_snr_arr = np.vstack((p_snr_arr,np.zeros(7)))
		

p_map_arr = p_map_arr[1:]
p_snr_arr = p_snr_arr[1:]
#p_map_arr = p_map_arr.T
#p_snr_arr = p_snr_arr.T

#print(cluster_arr.shape)
#print(pixx_arr.shape)
#print(p_map_arr.shape)
#print(p_snr_arr.shape)

df = np.vstack((cluster_arr,vlass_name_arr,vlass_id_arr,ra_arr,dec_arr,pixx_arr,pixy_arr,p_map_arr.T,p_snr_arr.T))

#print(df.shape)

df = pd.DataFrame(df.T, columns=['cluster','Component_name','Component_id','ra_vlass','dec_vlass','pix_x','pix_y','amp_map','x_map','y_map','sigma_x_map','sigma_y_map','theta_map','offset_map','amp_snr','x_snr','y_snr','sigma_x_snr','sigma_y_snr','theta_snr','offset_snr'])

df = df.astype(dtype={'cluster':str,'Component_name':str,'Component_id':float,'ra_vlass':float,'dec_vlass':float,'pix_x':float,'pix_y':float,'amp_map':float,'x_map':float,'y_map':float,'sigma_x_map':float,'sigma_y_map':float,'theta_map':float,'offset_map':float,'amp_snr':float,'x_snr':float,'y_snr':float,'sigma_x_snr':float,'sigma_y_snr':float,'theta_snr':float,'offset_snr':float})

df_joined = df.merge(all_vlass,how="left",left_on=["cluster","Component_name","ra_vlass","dec_vlass"],right_on=["cluster","Component_name","RA_2","DEC_2"])

df_joined.to_csv("reverse_search"+reduc+".csv")

