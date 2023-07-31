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
import argparse as ap

#RUN THIS TO FIND ALL THE CLUSTER MAPS WITH THESE DATA REDUCTION PARAMS
#find . -type f -name *_2asp_pca3_qm2_fitel_0f11-to-25f5Hz_qc_1p2rr_M_PdoCals_dt20_snr_iter1.fits > /home/users/ksarmien/Documents/clusters_substructure/out.txt

parser = ap.ArgumentParser(description="Point source finder for M2 SZ maps")
parser.add_argument("-d","--dir_map",type=str,help="Directory of M2 maps")
parser.add_argument("-r","--reduc",type=str,help="Imaging params of M2 maps")
parser.add_argument("-re","--reduction",type=str,help="List of files with the maps")
parser.add_argument("-ns","--nsigma",type=float,default=4.0,help="Sigma threshold for point source finder")
parser.add_argument("-t1","--theta1",nargs="+",type=int,help="List of min Gaussian kernels")
parser.add_argument("-t2","--theta2",nargs="+",type=int,help="List of max Gaussian kernels")
parser.add_argument("-o","--outfile",type=str,help="Name of outfile")
args = parser.parse_args()
dir_map = args.dir_map
reduc = args.reduc
reduction = args.reduction
reduction_list = "/users/ksarmien/mmpsrc_project/reductions_lists/"+reduction
nsigma = args.nsigma
theta1 = args.theta1
theta2 = args.theta2
outfile = args.outfile

#dir_map = "/home/scratch/cromero/MUSTANG2/Reductions/ACT_Sources_2023_0f09-to-35f5_PCA0"
#reduc = "_2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1"
#reduction = reduc+"_files.txt"
#reduction_list = "/users/ksarmien/mmpsrc_project/reductions_lists/"+reduction
#nsigma = 4.

#red_codes = pd.read_csv("/users/ksarmien/mmpsrc_project/map_quality_tables/reductions_code.csv")
code = 20.0

def load_and_scale(cluster_file):
    mapD=cluster_file
    hdu_map = fits.open(mapD)[0]
    img_map = hdu_map.data
    hdu_weight_map = fits.open(mapD)[1]
    temp=np.copy(img_map)
    weight_map = hdu_weight_map.data
    pos_weight_map = weight_map[weight_map >0.]
    weight_median = np.median(pos_weight_map)
    temp[weight_map< 0.2*weight_median]= 0.0
    sd = np.std(temp)
    scale_factor = 1./sd
    scaled_snr = scale_factor * temp
    return img_map,scaled_snr,scale_factor

def load_map(clustername):
    filename = re.sub('snr','map',clustername)
    mapD=filename
    hdu_map = fits.open(mapD)[0]
    img_map = hdu_map.data
    world = wcs.WCS(hdu_map.header)
    crval_ra = hdu_map.header["CRVAL1"]
    crval_dec = hdu_map.header["CRVAL2"]
    cd_pix = np.abs(hdu_map.header["CD1_1"])
    return img_map,world,crval_ra,crval_dec,cd_pix


def load_noise(clustername):
    mapD = clustername
    hdu = fits.open(mapD)[1]
    img = hdu.data
    return img

def load_hits_map(clustername):
    filename = re.sub('snr','map',clustername)
    hdu = fits.open(filename)[1]
    img = hdu.data
    return img

def load_weight_mask(cluster_file):
    mapD = cluster_file
    hdu_weight_map = fits.open(mapD)[1]
    weight_map = hdu_weight_map.data
    pos_weight_map = weight_map[weight_map >0.]
    weight_median = np.median(pos_weight_map)
    weight_map[weight_map > 0.2*weight_median] = 0.0
    return weight_map

def load_snr(clustername):
    mapD = clustername
    hdu_map = fits.open(mapD)[0]
    img_map = hdu_map.data
    world = wcs.WCS(hdu_map.header)
    return img_map,world


maps_w_stripes=["ACT-CLJ0248.0-0331","ACT-CLJ2048.2+0238","ACT-CLJ0248.2-0216","ACT-CLJ0248.3+0122","ACT-CLJ0248.3-0337","ACT-CLJ0248.7-0019","ACT-CLJ0248.9-0328","ACT-CLJ0248.7-0019","ACT-CLJ0248.9-0328","ACT-CLJ0250.1+0008","ACT-CLJ0254.4-0538","ACT-CLJ0257.0-0300","ACT-CLJ0259.4-0524","ACT-CLJ0301.5-0248","ACT-CLJ1258.4+0043","ACT-CLJ1315.0+0121"]

def rms(dog_map,snr_map,signal_map):
    rms_dog = np.std(dog_map.flatten())
    rms_snr = np.std(snr_map.flatten())
    rms_signal = np.std(signal_map.flatten())
    return rms_dog, rms_snr,rms_signal

def twoD_Gaussian_elliptical(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    (x,y) = xdata_tuple
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
    return g.ravel()

def twoD_Gaussian(xdata_tuple, amplitude, xo, yo, sigma):
    (x,y) = xdata_tuple
    xo = float(xo)
    yo = float(yo)
    g = amplitude*np.exp( - ((x-xo)**2 + (y-yo)**2)/(2*sigma**2))
    return g.ravel()

def minimize(pars,xdata_tuple,map_ravel):
    res = map_ravel - twoD_Gaussian(xdata_tuple, *pars)
    return np.sqrt(np.mean(res**2))

def get_depth_map(clusterfile):
    clusterfile = re.sub("snr","map",clusterfile) #make sure I'm loading map and NOT snr
    mapD = clusterfile #load signal map
    hdu_map = fits.open(mapD)[0]
    img_map = hdu_map.data
    hdu_weight_map = fits.open(mapD)[1] #load weight map
    weight_map = hdu_weight_map.data
    pixel_per_degree = hdu_map.header["CD1_1"] #number of pixels there are in onre degree
    crpixx =hdu_map.header["CRPIX1"] #central PIX IN X
    crpixy =hdu_map.header["CRPIX2"] #central PIX IN Y
    central_weight = weight_map[int(crpixx),int(crpixy)] #get the weight at the center
    weight_map_scaled = central_weight / weight_map  #### or is it weight_map/central_weight??????
    pixel_per_degree = np.abs(pixel_per_degree)
    pixel_per_arcsec = pixel_per_degree * 3600. #fnumber of pixels per arcsec
    pixels_beam = pixel_per_arcsec * 9. # #number of pixels in 9 arcsec (fwhm for Mustang)
    sigma_pixels_beam = pixels_beam/(2*np.sqrt(2*np.log(2))) #convert fwhm to number of pixel
    smoothed = gaussian_filter(img_map,sigma_pixels_beam,order=0) #convolve with a gaussian kernel
    #code to mask outer regions of map (leave only central region)
    y,x = img_map.shape
    x = np.linspace(-1*(int(x/2)),int(x/2)-1,x)
    y = np.linspace(-1*(int(y/2)),int(y/2)-1,y)
    x,y = np.meshgrid(x,y)
    dist = np.sqrt(x**2 + y**2)
    central = smoothed[dist<sigma_pixels_beam]
    rms_central = np.std(central)
    depth_map = rms_central * np.sqrt(weight_map_scaled)
    return rms_central,depth_map


def psrc_finder(data,theta1,theta2,th):
    blob_list = np.empty([1,5])
    for i in range(len(theta1)):
        for j in range(len(theta2)):
            if theta2[j]>theta1[i]:
                blobs=feature.blob_dog(data,theta1[i],theta2[j],threshold=th)
                x_coord = np.array(blobs[:,0],int)
                y_coord = np.array(blobs[:,1],int)
                sigma_dog = np.array(blobs[:,2],float)
                theta1_arr = np.repeat(theta1[i],len(x_coord))
                theta2_arr = np.repeat(theta2[j],len(y_coord))
                try:
                    b = np.array(list(zip(y_coord,x_coord,sigma_dog,theta1_arr,theta2_arr)))
                    blob_list = np.vstack((blob_list,b))
                except ValueError:
                    pass
    blob_list = blob_list[1:]
    if len(blob_list)>0:
        running = True
        return blob_list,running
    else:
        running=False
        return None,running


def coords_blobs(blob_list,world,central_coord):
    coords_ra_dec=world.wcs_pix2world(blob_list[:,:2],1)
    coords_ra_dec = np.array(coords_ra_dec)
    sky_coord_blob = SkyCoord(ra=coords_ra_dec[:,0]*u.degree,dec=coords_ra_dec[:,1]*u.degree,frame="icrs")
    sep = np.array(central_coord.separation(sky_coord_blob).radian)
    return np.column_stack((coords_ra_dec,sep))

def fitting_blobs(signal,snr,blob_list,cd_pix):
    area_pix = (np.abs(cd_pix)*3600)**2
    area_beam = 140.0
    xlen,ylen=signal.shape
    x = np.linspace(0,xlen-1,xlen)
    y = np.linspace(0,ylen-1,ylen)
    x,y=np.meshgrid(x,y)
    param_list = np.empty([1,11])
    signal_ravel = signal.ravel()
    signal_copy = signal_ravel.copy()
    snr_ravel = snr.ravel()
    for src in blob_list:
        fit = opt.minimize(minimize,[signal[int(src[1]),int(src[0])],src[0],src[1],3],args=((x,y),signal_copy))
        fit_snr = opt.minimize(minimize,[snr[int(src[1]),int(src[0])],src[0],src[1],3],args=((x,y),snr_ravel))
        p = fit.x
        p_snr = fit_snr.x
        g_r_snr = twoD_Gaussian((x,y),p_snr[0],p_snr[1],p_snr[2],p_snr[3])
        g_r = twoD_Gaussian((x,y),p[0],p[1],p[2],p[3])
        int_flux = np.sum(g_r)*area_pix/area_beam
        flux_err = np.sqrt(np.diag(fit.hess_inv)[0]**2)
        int_snr = np.sum(g_r_snr)*area_pix/area_beam
        g_r_snr = g_r_snr.reshape(xlen,ylen)
        snr = snr - g_r_snr
        p = np.append(p,int_flux)
        p = np.append(p,flux_err)
        p = np.append(p,p_snr)
        p = np.append(p,int_snr)
        param_list = np.vstack((param_list,p))
    param_list = param_list[1:]
    return param_list,snr

def point_srcs(clustername,theta1,theta2,nsigma):
    snr_original,snr_scaled,factor=load_and_scale(clustername)
    signal_map,world,crval_ra,crval_dec,cd_pix = load_map(clustername)
    hits_map = load_hits_map(clustername)
    noise_map = load_noise(clustername)
    central_coord = SkyCoord(ra = crval_ra*u.degree,dec = crval_dec*u.degree,frame="icrs")
    th = np.std(snr_original)*nsigma #five sigma? 
    substitutes = [dir_map,"/ACT-CLJ\d+.\d+(\+|-)\d+/Jy_",".fits","/ACT_Sources_2023_0f09-to-35f5_PCA0","/Jy_","/ACT_Sources_2023_0f09-to-35f5_PCA0","Jy_","_2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1","_2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20AGBT21B","_2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20AGBT22B","_snr_iter1"]
    cluster = clustername
    for s in substitutes:
        cluster = re.sub(s,"",cluster)
    print(cluster)
    running = True
    snr_copy = snr_original.copy()
    blob_list_tot = np.empty([1,5])
    param_list_tot = np.empty([1,11])
    iters = 0
    while running and iters<1:
        print("iter")
        blob_list,running = psrc_finder(snr_copy,theta1,theta2,th)
        iters += 1
        if running:
            blob_list_tot = np.vstack((blob_list_tot,blob_list))
            param_list,snr_copy = fitting_blobs(signal_map,snr_copy,blob_list,cd_pix)
            param_list_tot = np.vstack((param_list_tot,param_list))
    blob_list_tot = blob_list_tot[1:]
    param_list_tot = param_list_tot[1:]
    if len(blob_list_tot) > 0:
        cluster_list = list(np.repeat(cluster,len(blob_list_tot)))
        coords = coords_blobs(blob_list_tot,world,central_coord)
        ps_tot = np.empty([1,4])
        for src in blob_list_tot:
            ps_val = snr_original[int(src[1]),int(src[0])]
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
psrc_list = np.empty([1,24])

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

df_psrcs = pd.DataFrame(psrc_list[1:],columns = ['cluster', 'x', 'y','sigma_dog','theta_1','theta_2', 'ra_deg', 'dec_deg', 'dist_center_radians','amp_fit', 'x_center_fit', 'y_center_fit', 'sigma','int_flux_Jy','int_flux_err_Jy','amp_snr','x_snr','y_snr','sigma_snr','int_snr','snr','masked','noise_ps','hits_ps'])
df_psrcs = df_psrcs.astype(dtype={'cluster':str,'x':float,'y':float,'sigma_dog':float,'theta_1':float,'theta_2':float,'ra_deg':float,'dec_deg':float,'dist_center_radians':float,'amp_fit':float,'x_center_fit':float,'y_center_fit':float,'sigma':float,'int_flux_Jy':float,'int_flux_err_Jy':float,'amp_snr':float,'x_snr':float,'y_snr':float,'sigma_snr':float,'int_snr':float,'snr':float,'masked':float,'noise_ps':float,'hits_ps':float})

df_quality = pd.read_csv("/users/ksarmien/mmpsrc_project/map_quality_tables/data_quality_code_20.csv")
df_quality = df_quality.loc[df_quality["red_type"]==code]
df_psrcs = pd.merge(df_psrcs,df_quality,how="left",left_on="cluster",right_on="Source")

df_psrcs.to_csv(outfile,index=False)
