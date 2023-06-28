import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

m_vlass = pd.read_csv("psrc_lists/matched_9arcsec_1and2_best1_vlass.csv").add_suffix("_vlass")
m_first = pd.read_csv("psrc_lists/matched_9arcsec_1and2_best1_first.csv").add_suffix("_first")
m_gleam = pd.read_csv("psrc_lists/matched_100arcsec_1and2_best1_gleam.csv").add_suffix("_gleam")
m_tgss = pd.read_csv("psrc_lists/matched_25arcsec_1and2_best1_tgss.csv").add_suffix("_tgss")
#m_lotss = pd.read_csv("psrc_lists/matched_9arcsec_1and2_best1_lotss.csv").add_suffix("_lotss")
m_2mass = pd.read_csv("psrc_lists/matched_9arcsec_1and2_best1_2mass.csv").add_suffix("_2mass")
m_wise = pd.read_csv("psrc_lists/matched_9arcsec_1and2_best1_wise.csv").add_suffix("_wise")

blind = pd.read_csv("psrc_lists/uncleaned_psrcs_sigma__2aspcmsubqm2_fitel_0f09-to-35f5Hz_qc_0p6rr_M_PdoCals_dt20_snr_iter1_4.0sigma.csv")
all_matched = blind
tabs_list = [m_vlass,m_first,m_gleam,m_tgss,m_2mass,m_wise]
suffixes_right=["_vlass","_first","_gleam","_tgss","_2mass","_wise"]
cols_to_match = ["cluster","ra_deg","dec_deg","theta_1","theta_2"]

for i in range(len(tabs_list)):
    right_on_cols = [s + suffixes_right[i] for s in cols_to_match]
    all_matched = all_matched.merge(tabs_list[i],how="outer",left_on=cols_to_match,right_on=right_on_cols)

f_0_w1 = 309.540 #Jy
f_0_w2 = 171.787 #Jy
f_0_w3 = 31.674#Jy
f_0_w4 = 8.363 #Jy
all_matched["w1_flux_wise"] = f_0_w1 * 10**(-1*all_matched["w1mpro_wise"]/2.5) * 1000 #mJy
all_matched["w2_flux_wise"] = f_0_w2 * 10**(-1*all_matched["w2mpro_wise"]/2.5) * 1000 #mJy
all_matched["w3_flux_wise"] = f_0_w3 * 10**(-1*all_matched["w3mpro_wise"]/2.5) * 1000 #mJy
all_matched["w4_flux_wise"] = f_0_w4 * 10**(-1*all_matched["w4mpro_wise"]/2.5) * 1000 #mJy

f_0_2mass_j = 1590 #Jy
f_0_2mass_h = 1020 #Jy
f_0_2mass_k = 667 # Jy
all_matched["j_flux_2mass"] = f_0_2mass_j * 10**(-1*all_matched["j_m_2mass"]/2.5) * 1000 #mJy
all_matched["h_flux_2mass"] = f_0_2mass_h * 10**(-1*all_matched["h_m_2mass"]/2.5) * 1000 #mJy
all_matched["k_flux_2mass"] = f_0_2mass_k * 10**(-1*all_matched["k_m_2mass"]/2.5) * 1000 #mJy

all_matched.to_csv("psrc_lists/all_matched.csv")

