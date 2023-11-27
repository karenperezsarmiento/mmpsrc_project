import pandas as pd
import numpy as np
from astropy.io import fits
import astropy.coordinates as coord
from astropy import wcs
import scipy.optimize as opt
from astropy.coordinates import Angle
import re

def twoD_Gaussian(xdata_tuple, amplitude, xo, yo, sigma):
    (x,y) = xdata_tuple
    xo = float(xo)
    yo = float(yo)
    g = amplitude*np.exp( - ((x-xo)**2 + (y-yo)**2)/(2*sigma**2))
    return g.ravel()

def minimize(pars,xdata_tuple,map_ravel):
    res = map_ravel - twoD_Gaussian(xdata_tuple, *pars)
    return np.sqrt(np.mean(res**2))

direc = "/home/scratch/sdicker/AGBT21B_298/real_maps/"
reduc = "_2aspcmsubqm2_fitel_0f05-to-49f5Hz_qc_0p6rr_M_PdoCals_dt10_snr.fits"
df = pd.read_csv("../psrc_lists/unfiltered_psrcs_2aspcmsubqm2_fitel_0f05-to-49f5Hz_qc_0p6rr_M_PdoCals_dt10_snr_5.0sigma.csv")

theta_1 = 2.0
theta_2 = 6.0
df = df.loc[(df["theta_1"]==theta_1)&(df["theta_2"]==theta_2)]
area_beam = 140.0

res = np.empty((1,10))
for i in range(len(df)):
	fits_snr = direc+"Jy_"+np.array(df["cluster"])[i]+reduc
	fits_map = re.sub("snr","map",fits_snr)
	hdu_map = fits.open(fits_map)[0]
	hdu_snr = fits.open(fits_snr)[0]
	img_map = hdu_map.data
	img_snr = hdu_snr.data
	world = wcs.WCS(hdu_map.header)
	crval_ra = hdu_map.header["CRVAL1"]
	crval_dec = hdu_map.header["CRVAL2"]
	cd_pix = np.abs(hdu_map.header["CD1_1"])
	area_pix = (np.abs(cd_pix)*3600)**2
	s_0 = np.array(df["x"])[i]
	s_1 = np.array(df["y"])[i]
	signal_ravel = img_map.ravel()
	snr_ravel= img_snr.ravel()
	xlen,ylen=img_map.shape
	x = np.linspace(0,xlen-1,xlen)
	y = np.linspace(0,ylen-1,ylen)
	x,y=np.meshgrid(x,y)
	fit = opt.minimize(minimize,[img_map[int(s_1),int(s_0)],s_0,s_1,3],args=((x,y),signal_ravel))
	p = fit.x
	g_r = twoD_Gaussian((x,y),p[0],p[1],p[2],p[3])
	int_flux = np.sum(g_r)*area_pix/area_beam
	fit_snr = opt.minimize(minimize,[img_snr[int(s_1),int(s_0)],s_0,s_1,3],args=((x,y),snr_ravel))
	p_snr = fit_snr.x
	g_r_snr = twoD_Gaussian((x,y),p_snr[0],p_snr[1],p_snr[2],p_snr[3])
	int_snr = np.sum(g_r_snr)*area_pix/area_beam
	res_i = np.hstack((p,int_flux,p_snr,int_snr))
	res = np.vstack((res,res_i))

res = res[1:]

df[["amp_fit_new","x_center_fit_new","y_center_fit_new","sigma_new","int_flux_Jy_new","amp_snr_new","x_snr_new","y_snr_new","sigma_snr_new","int_snr_new"]] = res

df.to_csv("../psrc_lists/new_fitting_no_offset.csv")


	

