from skimage import data, feature
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import matplotlib.pyplot as plt
import numpy as np
from astropy.convolution import convolve
from astropy.coordinates import SkyCoord
import astropy.coordinates as coord
from astropy import wcs
import astropy.units as u
from scipy import optimize
import pandas as pd
import time
import os
import re
from astropy.coordinates import Angle
from astroquery.simbad import Simbad
from astroquery.nrao import Nrao as Nrao
from astroquery.ned import Ned
import scipy.optimize as opt
from astropy.convolution import Gaussian2DKernel
from skimage import filters
from scipy.ndimage import gaussian_filter
from scipy import stats
import psrc_lib as pslib

#RUN THIS TO FIND ALL THE CLUSTER MAPS WITH THESE DATA REDUCTION PARAMS
#find . -type f -name *_2asp_pca3_qm2_fitel_0f11-to-25f5Hz_qc_1p2rr_M_PdoCals_dt20_snr_iter1.fits > /home/users/ksarmien/Documents/clusters_substructure/out.txt

dir_map = "/home/scratch/sdicker/AGBT21B_298"
reduc = "_2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1"
reduction = "injected.txt"
reduction_list = "/users/ksarmien/mmpsrc_project/reductions_lists/"+reduction
nsigma = 5.

def point_srcs(clustername,theta1,theta2,nsigma):    
    snr_original,snr_scaled,factor=load_and_scale(clustername)
    signal_map,world,crval_ra,crval_dec = load_map(clustername)
    hits_map = load_hits_map(clustername)
    noise_map = load_noise(clustername)
    central_coord = SkyCoord(ra = crval_ra*u.degree,dec = crval_dec*u.degree,frame="icrs")
    th = np.std(snr_original)*nsigma #five sigma? 
    substitute1 = dir_map
    substitute2 = "/Jy_"
    substitute3 = reduc + ".fits"
    substitute4 = "/ACT_Sources_2023_0f09-to-35f5_PCA0"
    key=re.sub(substitute1,"",clustername)
    key=re.sub(substitute2,"",key)
    key=re.sub(substitute3,"",key)
    cluster =re.sub(substitute4,"",key)
    running = True
    snr_copy = snr_original.copy()
    blob_list_tot = np.empty([1,5])
    param_list_tot = np.empty([1,10])
    while running:
        print("iter")
        blob_list,running = psrc_finder(snr_copy,theta1,theta2,th)
        if running:
            blob_list_tot = np.vstack((blob_list_tot,blob_list))
            param_list,snr_copy = fitting_blobs(signal_map,snr_copy,blob_list)
            param_list_tot = np.vstack((param_list_tot,param_list))
    blob_list_tot = blob_list_tot[1:]
    param_list_tot = param_list_tot[1:]
    if len(blob_list_tot) > 0:
        cluster_list = list(np.repeat(cluster,len(blob_list_tot)))
        coords = coords_blobs(blob_list_tot,world,central_coord)
        ps_tot = np.empty([1,4])
        for src in blob_list_tot:
            ps_val = snr_scaled[int(src[1]),int(src[0])]
            ps_mask = bool(ps_val == 0.0)
            ps_noise = np.sum(noise_map[int(src[1])-2:int(src[1])+2,int(src[0])-2:int(src[0])+2])/16
            ps_hits = np.sum(hits_map[int(src[1])-2:int(src[1])+2,int(src[0])-2:int(src[0])+2])/16
            ps = np.array([ps_val,ps_mask,ps_noise,ps_hits])
            ps_tot = np.vstack((ps_tot,ps))
        ps_tot = ps_tot[1:]
        result = np.column_stack((np.array(cluster_list),blob_list_tot,coords,param_list_tot,ps_tot))
    else:
        result = None
    return result

all_snr_files = []
with open(reduction_list) as f:
    for line in f:
        l = dir_map+line[1:]
        l = l[:-1]
        all_snr_files.append(l)
theta1 = [2]
theta2 = [6]
psrc_list = np.empty([1,23])

t0=time.time()
for i in all_snr_files:
    print(i)
    res= point_srcs(i,theta1,theta2,nsigma)
    try:
        length = len(res)>0
        psrc_list = np.vstack((psrc_list,res))
    except TypeError:
        pass
t1=time.time()
t = t1-t0
print(t)

df_psrcs = pd.DataFrame(psrc_list[1:],columns = ['cluster', 'x', 'y','sigma_dog','theta_1','theta_2', 'ra_deg', 'dec_deg', 'dist_center_radians','amp_fit', 'x_center_fit', 'y_center_fit', 'sigma','int_flux_Jy','int_flux_err_Jy','amp_snr','x_snr','y_snr','sigma_snr','snr','masked','noise_ps','hits_ps'])
df_psrcs = df_psrcs.astype(dtype={'cluster':str,'x':float,'y':float,'sigma_dog':float,'theta_1':float,'theta_2':float,'ra_deg':float,'dec_deg':float,'dist_center_radians':float,'amp_fit':float,'x_center_fit':float,'y_center_fit':float,'sigma':float,'int_flux_Jy':float,'int_flux_err_Jy':float,'amp_snr':float,'x_snr':float,'y_snr':float,'sigma_snr':float,'snr':float,'masked':float,'noise_ps':float,'hits_ps':float})

filename_1 = "/users/ksarmien/mmpsrc_project/psrc_lists/injected_2_6_psrcs_"+reduc+"_"+str(nsigma)+"sigma.csv"
df_psrcs.to_csv(filename_1,index=False)

