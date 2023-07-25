#!/bin/sh
source activate clusters_proj
python psrc_lib.py -d /home/scratch/cromero/MUSTANG2/Reductions/ACT_Sources_2023_0f09-to-35f5_PCA0 -r _2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1 -re _2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1_files.txt -ns 4.0 -t1 2 3 4 5 6 7 -t2 2 3 4 5 6 7 -o /users/ksarmien/mmpsrc_project/psrc_lists/uncleaned_psrcs_sigma__2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1_4.0sigma.csv 
