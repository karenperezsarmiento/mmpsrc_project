import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.optimize as opt

all_matched = pd.read_csv("psrc_lists/all_matched.csv")
f_all_matched = all_matched.loc[(all_matched["theta_1"]==2.0)&(all_matched["theta_2"]==6.0)]

freqs = [3e9,1.4e9,76e6,84e6,92e6,99e6,107e6,115e6,122e6,130e6,143e6,151e6,158e6,166e6,174e6,181e6,189e6,197e6,204e6,212e6,220e6,150e6,144e6,2.4275e14,1.8038e14,1.3886e14,8.901e13,6.492e13,2.4813e13,1.3508e13]
markers=["+","s","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","x","v","^","*","*","*","8","8","8","8"]
colors = ["orange" for i in range(len(markers))]
freqs_fit = np.log10(np.array([1.4e9,3e9,90e9],dtype=float))
sed_arr= np.empty([1,8])

def linear_func(x,m,b):
    return m*x+b

def minimize(pars,log10_freqs,log10_fluxes):
    res = log10_fluxes - linear_func(log10_freqs,*pars)
    return np.sqrt(np.mean(res**2))

for i in range(len(f_all_matched)):
    m2_flux = f_all_matched.iloc[i]["int_flux_Jy"]*1000
    cat_fluxes = np.array(list(f_all_matched.iloc[i][["Total_flux_vlass","FINT_first","int_flux_076_gleam","int_flux_084_gleam","int_flux_092_gleam","int_flux_099_gleam","int_flux_107_gleam","int_flux_115_gleam","int_flux_122_gleam","int_flux_130_gleam","int_flux_143_gleam","int_flux_151_gleam","int_flux_158_gleam","int_flux_166_gleam","int_flux_174_gleam","int_flux_181_gleam","int_flux_189_gleam","int_flux_197_gleam","int_flux_204_gleam","int_flux_212_gleam","int_flux_220_gleam","Total_flux_tgss","Total_flux_lotss","j_m_2mass","h_m_2mass","k_m_2mass","w1_flux_wise","w2_flux_wise","w3_flux_wise","w4_flux_wise"]]))
    cat_fluxes_fit = np.log10(np.append(np.array(list(f_all_matched.iloc[i][["FINT_first","Total_flux_vlass"]]),dtype=float),m2_flux))
    if np.sum(np.isnan(cat_fluxes_fit))<2:
        p = opt.minimize(minimize,[-0.7,0.0],args=(freqs_fit[~np.isnan(cat_fluxes_fit)],cat_fluxes_fit[~np.isnan(cat_fluxes_fit)])).x
    else:
        p = np.repeat(999,2)
    row=np.append(np.array(list(f_all_matched.iloc[i][["cluster","ra_deg","dec_deg","FINT_first","Total_flux_vlass","int_flux_Jy"]])),p)
    sed_arr = np.vstack((sed_arr,row))
    print(i)
    if len(cat_fluxes) == np.sum(np.isnan(cat_fluxes)):
        pass
    else:
        fig = plt.figure()
        plt.scatter(90e9,m2_flux)
        for k in range(len(freqs)):
            plt.scatter(freqs[k],cat_fluxes[k],marker=markers[k],c=colors[k])
        if p[0]!=999:
            plt.plot(10**freqs_fit[~np.isnan(cat_fluxes_fit)],10**(p[0]*freqs_fit[~np.isnan(cat_fluxes_fit)]+p[1]),linestyle='dashed',c="r")
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel("Frequency (Hz)")
        plt.ylabel("Flux density (mJy)")
        title = str(f_all_matched.iloc[i]["cluster"])+"\n RA: "+str(f_all_matched.iloc[i]["ra_deg"])+" DEC: "+str(f_all_matched.iloc[i]["dec_deg"])
        filename=str(f_all_matched.iloc[i]["cluster"])+"_ra_"+str(f_all_matched.iloc[i]["ra_deg"])+"_dec_"+str(f_all_matched.iloc[i]["dec_deg"])
        plt.title(title)
        plt.savefig("sed_psrcs/"+filename+".png")
        plt.close(fig)

sed_arr = sed_arr[1:]
df_sed = pd.DataFrame(sed_arr,columns=["cluster","ra_deg","dec_deg","FINT_first","Total_flux_vlass","int_flux_Jy","alpha","sed_offset"])

df_sed.to_csv("spectral_indices.csv")
