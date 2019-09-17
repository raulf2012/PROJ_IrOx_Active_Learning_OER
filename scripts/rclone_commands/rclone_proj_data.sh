#!/bin/bash

rclone --transfers=2 sync \
raul_dropbox:01_norskov/PROJECT_DATA/04_IrOx_surfaces_OER \
$norskov_research/PROJECT_DATA/04_IrOx_surfaces_OER \
--exclude "190321_new_job_df/**" \
--exclude "181226_new_job_df/**" \
--exclude "190103_new_job_df/**" \
--exclude "190115_new_job_df/**" \
--exclude "190313_new_job_df/**" \
--exclude "190315_new_job_df/**" \
--exclude "190318_new_job_df/**" \
--exclude "bulk_systems/**" \
--verbose \
