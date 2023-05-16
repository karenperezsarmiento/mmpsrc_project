import numpy as np
import scipy
import pandas as pd

def twoD_Gaussian(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    (x,y) = xdata_tuple
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
    return g.ravel()

df_psrcs = pd.read_csv("psrc_lists/uncleaned_psrcs_sigma__2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1_4.0sigma.csv")

xlen = 324
ylen = 324
x = np.linspace(0,xlen-1,xlen)
y = np.linspace(0,ylen-1,ylen)
x,y=np.meshgrid(x,y)

flux_1 = []
flux_2 = []
cluster_arr = []
ra = []
dec = []

for i in range(len(df_psrcs)):
	p = np.array(list(df_psrcs.iloc[i][['amp_fit', 'x_center_fit', 'y_center_fit', 'sigma_x_fit','sigma_y_fit','theta','offset']]))
	g_r_1 = np.sum(twoD_Gaussian((x,y),p[0],p[1],p[2],p[3],p[4],p[5],p[6]))
	g_r_2 = np.sum(twoD_Gaussian((x,y),p[0],p[1],p[2],p[3],p[4],p[5],0))
	flux_1 = np.append(flux_1,g_r_1)
	flux_2 = np.append(flux_2,g_r_2)
	cluster_arr = np.append(cluster_arr,df_psrcs.iloc[i]["cluster"])
	ra = np.append(ra,df_psrcs.iloc[i]["ra_deg"])
	dec = np.append(dec,df_psrcs.iloc[i]["dec_deg"])

df_flux = pd.DataFrame(np.array(list(zip(cluster_arr,ra,dec,flux_1,flux_2))),columns=["cluster","ra_deg","dec_deg","flux_offset","flux_no_offset"])

df_flux.to_csv("fluxes_text.csv")
