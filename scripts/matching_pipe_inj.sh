#!/bin/sh
topcat -stilts tmatch2 matcher=sky params=18 join=all1 find=all out=/users/ksarmien/mmpsrc_project/psrc_lists/injected_found.csv in1=/users/ksarmien/mmpsrc_project/psrc_lists/fake_pnt_amplitudes.csv in2=/users/ksarmien/mmpsrc_project/psrc_lists/injected_psrcs_2_6_2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1_4.0sigma.csv values1='RA DEC' values2='ra_deg dec_deg'
