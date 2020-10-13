#!/bin/bash

rclone --transfers=2 sync \
raul_dropbox:01_norskov/00_git_repos/PROJ_IrOx_Active_Learning_OER \
$norskov_research/00_git_repos/PROJ_IrOx_Active_Learning_OER \
--exclude ".git/**" \
--exclude "*voro_temp*" \
--exclude "00_abx_al_runs/out_data/**" \
--exclude "01_abx_al_runs_new/out_data/**" \
--exclude "__old__/**" \
--verbose \

# rclone --transfers=2 sync \
# raul_dropbox:01_norskov/00_git_repos/PROJ_IrOx_Active_Learning_OER \
# $norskov_research/00_git_repos/PROJ_IrOx_Active_Learning_OER \
# --exclude ".git/**" \
# --exclude "*voro_temp*" \
# --exclude "*01_abx_al_runs_new/out_data*" \
# --verbose \
  

# --exclude "190321_new_job_df/**" \
# --exclude "bulk_systems/**" \
