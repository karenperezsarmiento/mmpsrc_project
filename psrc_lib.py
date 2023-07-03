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

#RUN THIS TO FIND ALL THE CLUSTER MAPS WITH THESE DATA REDUCTION PARAMS
#find . -type f -name *_2asp_pca3_qm2_fitel_0f11-to-25f5Hz_qc_1p2rr_M_PdoCals_dt20_snr_iter1.fits > /home/users/ksarmien/Documents/clusters_substructure/out.txt

dir_map = "/home/scratch/cromero/MUSTANG2/Reductions/ACT_Sources_2023_0f09-to-35f5_PCA0"
reduc = "_2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1"
reduction = reduc+"_files.txt"
reduction_list = "/users/ksarmien/mmpsrc_project/reductions_lists/"+reduction
nsigma = 4.

red_codes = pd.read_csv("/users/ksarmien/mmpsrc_project/map_quality_tables/reductions_code.csv")
code = np.array(pd.to_numeric(red_codes["code"][red_codes["reduction"] == reduction]))[0]

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
    return img_map,world,crval_ra,crval_dec


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

def point_srcs(clustername,theta1,theta2,nsigma):    
    #blob_dict={}
    blob_list = np.empty([1,5])
    snr_original,snr_scaled,factor=load_and_scale(clustername)
    signal_map,world,crval_ra,crval_dec = load_map(clustername)
    hits_map = load_hits_map(clustername)
    noise_map = load_noise(clustername)
    central_coord = SkyCoord(ra = crval_ra*u.degree,dec = crval_dec*u.degree,frame="icrs")
    #rms_central,depth_map = get_depth_map(clustername)
    #nt,psd,wf = check_bad_weather(snr_original)
    th = np.std(snr_original)*nsigma #five sigma? 
    for i in range(len(theta1)): 
        for j in range(len(theta2)):
            if theta2[j]>theta1[i]:
                blobs=feature.blob_dog(snr_original,theta1[i],theta2[j],threshold=th)
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
    substitute1 = dir_map
    substitute2 = "/ACT-CLJ\d+.\d+(\+|-)\d+/Jy_"
    substitute3 = reduc + ".fits"
    substitute4 = "/ACT_Sources_2023_0f09-to-35f5_PCA0"
    key=re.sub(substitute1,"",clustername)
    key=re.sub(substitute2,"",key)
    key=re.sub(substitute3,"",key)
    cluster =re.sub(substitute4,"",key)
    #param = 11
    blob_list = blob_list[1:]
    xlen,ylen=signal_map.shape
    x = np.linspace(0,xlen-1,xlen)
    y = np.linspace(0,ylen-1,ylen)
    x,y=np.meshgrid(x,y)
    signal_map_ravel = signal_map.ravel()
    snr_ravel = snr_original.ravel()
    if len(blob_list) > 0:
        cluster_list = list(np.repeat(cluster,len(blob_list)))
        #param_map_list = list(np.repeat(param,len(blob_list)))
        coords_ra_dec=world.wcs_pix2world(blob_list[:,:2],1)
        coords_ra_dec = np.array(coords_ra_dec)
        sky_coord_blob = SkyCoord(ra=coords_ra_dec[:,0]*u.degree,dec=coords_ra_dec[:,1]*u.degree,frame="icrs")
        sep = np.array(central_coord.separation(sky_coord_blob).radian) 
        param_list = np.empty([1,14])
        for src in blob_list:
            ps_val = snr_scaled[int(src[1]),int(src[0])]
            ps_mask = bool(ps_val == 0.0)
            ps_noise = np.sum(noise_map[int(src[1])-2:int(src[1])+2,int(src[0])-2:int(src[0])+2])/16
            ps_hits = np.sum(hits_map[int(src[1])-2:int(src[1])+2,int(src[0])-2:int(src[0])+2])/16
            #initial_guess = (11,src[0],src[1],3,3,0,0)
            #ps_depth = depth_map[int(src[1]),int(src[0])]
            #ps_nt = nt
            #ps_psd = psd
            #ps_wf = wf
            try:
                #p = opt.minimize(minimize,[signal_map[int(src[1]),int(src[0])],src[0],src[1],3,3,0,0],args=((x,y),signal_map_ravel)).x
                fit = opt.minimize(minimize,[signal_map[int(src[1]),int(src[0])],src[0],src[1],3],args=((x,y),signal_map_ravel))
                p = fit.x
                flux_err = np.sqrt(np.sum(np.diag(fit.hess_inv)))
                #g_r = twoD_Gaussian_elliptical((x,y),p[0],p[1],p[2],p[3],p[4],p[5],p[6])
                g_r = twoD_Gaussian((x,y),p[0],p[1],p[2],p[3])
                int_flux = np.sum(g_r)
                #p_snr = opt.minimize(minimize,[snr_original[int(src[1]),int(src[0])],src[0],src[1],3,3,0,0],args=((x,y),snr_ravel)).x
                p_snr = opt.minimize(minimize,[snr_original[int(src[1]),int(src[0])],src[0],src[1],3],args=((x,y),snr_ravel)).x
            except RuntimeError:
                p = np.repeat(0,4)
                int_flux = 0
                flux_err = 0
                p_snr = np.repeat(0,4)
            #p = gaussianFit(signal_map,src[1],src[0])
            #p = np.append(p,inj_bool)
            p = np.append(p, int_flux)
            p = np.append(p,flux_err)
            p = np.append(p,p_snr)
            p = np.append(p,ps_val)
            p = np.append(p,ps_mask)
            p = np.append(p,ps_noise)
            p = np.append(p,ps_hits)
            #p = np.append(p,rms_central)
            #p = np.append(p,ps_depth)
            #p = np.append(p,ps_nt)
            #p =  np.append(p,ps_psd)
            #p = np.append(p,ps_wf)
            param_list = np.vstack((param_list,p))

            
        param_list = param_list[1:]
        #print("param list "+str(param_list.shape))
        #print("cluster list "+str(np.array(cluster_list).shape)) 
        #print("blob list "+str(blob_list.shape))
        #print("coords list"+str(coords_ra_dec.shape))
	#ident = np.array(list(zip(cluster_list,param_map_list)))
        result = np.column_stack((np.array(cluster_list),blob_list,coords_ra_dec,sep,param_list))
        #print(param_map_list)
    else:
        result = None
    #print(time.time())
    return result

all_snr_files = []
with open(reduction_list) as f:
    for line in f:
        l = dir_map+line[1:]
        l = l[:-1]
        all_snr_files.append(l)
theta1 = [2,3,4,5,6,7]
theta2 = [2,3,4,5,6,7]
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

df_quality = pd.read_csv("/users/ksarmien/mmpsrc_project/map_quality_tables/data_quality_code_20.csv")
df_quality = df_quality.loc[df_quality["red_type"]==code]
df_psrcs = pd.merge(df_psrcs,df_quality,how="left",left_on="cluster",right_on="Source")

filename_1 = "/users/ksarmien/mmpsrc_project/psrc_lists/uncleaned_psrcs_sigma_"+reduc+"_"+str(nsigma)+"sigma.csv"
df_psrcs.to_csv(filename_1,index=False)

