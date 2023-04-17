import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from psrc_lib import twoD_Gaussian

dir_map = "/home/scratch/cromero/MUSTANG2/Reductions/ACT_Sources_2023_0f09-to-35f5_PCA0"
reduc = "_2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1"
reduction = reduc+"_files.txt"
reduction_list = "/users/ksarmien/mmpsrc_project/reductions_lists/"+reduction
nsigma = 4.0

filename_1 = "/users/ksarmien/mmpsrc_project/psrc_lists/uncleaned_psrcs_sigma_"+reduc+"_"+str(nsigma)+"sigma.csv"

psrc_df = pd.read_csv(filename_1)
clusterlist = pd.read_csv(reduction_list)

for i in clusterlist:
	l = 
