#!/bin/sh
source activate clusters_proj
python make_injected_maps.py
/bin/bash /users/ksarmien/mmpsrc_project/scripts/psrc_inj.sh
/bin/bash /users/ksarmien/mmpsrc_project/script/matching_pipe_inj.sh
python make_psrc_imgs_injected.py 
