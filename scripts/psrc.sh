#!/bin/sh
source activate clusters_proj
python psrc_lib.py -d /home/scratch/sdicker/AGBT21B_298/real_maps -r _2aspcmsubqm2_fitel_0f05-to-49f5Hz_qc_0p6rr_M_PdoCals_dt10_snr -re _2aspcmsubqm2_fitel_0f05-to-49f5Hz_qc_0p6rr_M_PdoCals_dt10_snr_files.txt -ns 5.0 -t1 2 3 4 5 6 7 -t2 2 3 4 5 6 7 -o /users/ksarmien/mmpsrc_project/psrc_lists/unfiltered_psrcs_2aspcmsubqm2_fitel_0f05-to-49f5Hz_qc_0p6rr_M_PdoCals_dt10_snr_5.0sigma.csv 
