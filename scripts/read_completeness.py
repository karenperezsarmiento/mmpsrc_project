from scipy.io import readsav
import numpy as np
from scipy import interpolate

wkdir='/home/scratch/sdicker/AGBT21B_298/'
inFile='All_results_sorted.dat'

def read_comp_bias():
	complete=readsav(wkdir+inFile,python_dict=True)#read in the data
	#fields amps / sim_rings, the amplitudes and simulation ring boundries
	#found, bias - completeness and bias as a function of flux (average of all radius)
	#found_rad, bias_rad - completeness and bias as function of flux and radius bin
	#names - string array of names of radius bins
	radius=np.zeros(len(complete['sim_rings'])-1)
	for nn in range(len(complete['sim_rings'])-1):
	    radius[nn]=(complete['sim_rings'][nn]+complete['sim_rings'][nn+1])/2.
	bias=interpolate.RegularGridInterpolator((complete['amps'],radius),complete['bias_rad'],bounds_error=False,fill_value=None)#divide fluxes by this to get true debosted value
	completeness=interpolate.RegularGridInterpolator((complete['amps'],radius),complete['found_rad'],bounds_error=False,fill_value=None)#(flux,rad)
	#print(np.min(complete["amps"]))
	#print(np.max(complete["amps"]))
	#print(np.min(radius))
	#print(np.max(radius))
	return bias,completeness
