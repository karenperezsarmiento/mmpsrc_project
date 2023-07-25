#!/bin/sh
source activate clusters_proj
python center_cluster_maps.py
python data_quality.py
python make_maps.py
python noise_all.py
python snr_all.py
/bin/bash /users/ksarmien/mmpsrc_project/scripts/psrc.sh
python reverse_search.py
/bin/bash /users/ksarmien/mmpsrc_project/scripts/matching_pipe.sh
python match_all_cats.py
python final_selection_psrcs.py
python make_psrc_imgs.py
python sed_plots.py
