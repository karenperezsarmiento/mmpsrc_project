#!/bin/sh
source activate clusters_proj
python center_cluster_maps.py
python data_quality.py
python make_maps.py
python noise_all.py
python snr_all.py
python psrc_lib.py
python reverse_search.py
/bin/bash /users/ksarmien/mmpsrc_project/matching_pipe.sh
python match_all_cats.py
python final_selection_psrcs.py
python make_psrc_imgs.py
python sed_plots.py
python make_injected_maps.py
python psrc_injected_lib.py
python make_psrc_imgs_injected.py 
