import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

theta_1 = 2.0
theta_2 = 6.0
inner_rad = (3./60)*np.pi/180
max_rad = (4./60)*np.pi/180
snr_lower = 5.0
snr_higher = 7.0
noise_cut = 0.00065

full = pd.read_csv("../psrc_lists/unfiltered_psrcs_2aspcmsubqm2_fitel_0f05-to-49f5Hz_qc_0p6rr_M_PdoCals_dt10_snr_5.0sigma.csv")
matched = pd.read_csv("../psrc_lists/all_matched.csv",low_memory=False)

def cut_cat(df,theta_1,theta_2,inner_rad,snr_lower,snr_higher,max_rad):
    df = df.loc[(df["theta_1"]==theta_1)&(df["theta_2"]==theta_2)]
    inner = df.loc[(df["dist_center_radians"]<inner_rad)&(df["peak_snr"]>=snr_lower)]
    outer = df.loc[(df["dist_center_radians"]>=inner_rad)&(df["dist_center_radians"]<max_rad)&(df["peak_snr"]>=snr_higher)]
    final = pd.concat([inner,outer])
    return final


final = cut_cat(full,theta_1,theta_2,inner_rad,snr_lower,snr_higher,max_rad)
final_matched = cut_cat(matched,theta_1,theta_2,inner_rad,snr_lower,snr_higher,max_rad)

final.to_csv("../psrc_lists/final_cat_snr_5_3arcmin_snr_7.csv")
final_matched.to_csv("../psrc_lists/final_matched_snr_5_3arcmin_snr_7.csv")

final_paper = final[["cluster","ra_deg","dec_deg","dist_center_radians","amp_fit","int_flux_Jy","int_flux_err_Jy","peak_snr","peak_flux_mJy_adj"]]
final_matched_paper = final_matched[["cluster","ra_deg","dec_deg","dist_center_radians","amp_fit","int_flux_Jy","int_flux_err_Jy","peak_snr","peak_flux_mJy_adj","dist_center_radians","Component_name_vlass","Component_id_vlass","RA_vlass","DEC_vlass","Total_flux_vlass","RA_first","DEC_first","FINT_first","Name_gleam","RAJ2000_gleam","DEJ2000_gleam","int_flux_wide_gleam","Source_name_tgss","RA_tgss","DEC_tgss","Total_flux_tgss","ra_wise","dec_wise","w1_flux_wise","w2_flux_wise","w3_flux_wise","w4_flux_wise","ra_2mass","dec_2mass","j_flux_2mass","h_flux_2mass","k_flux_2mass"]]

final_paper.to_csv("../psrc_lists/PAPER_final_cat_snr_5_3arcmin_snr_7.csv")

final_matched_paper.to_csv("../psrc_lists/PAPER_final_matched_snr_5_3arcmin_snr_7.csv")
