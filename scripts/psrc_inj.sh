#!/bin/sh
source activate clusters_proj
python psrc_lib.py -d /home/scratch/sdicker/AGBT21B_298/ -r _2aspcmsubqm2_fitel_0f05-to-49f5Hz_qc_0p6rr_M_PdoCals_dt10AGBT21B_298 -re injected.txt -ns 4.0 -t1 2 -t2 6 -o /users/ksarmien/mmpsrc_project/psrc_lists/injected_psrcs_2_6_2aspcmsubqm2_fitel_0f05-to-49f5Hz_qc_0p6rr_M_PdoCals_dt10.csv 
