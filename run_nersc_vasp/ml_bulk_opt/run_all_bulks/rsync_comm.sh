rclone copy \
raul_dropbox:01_norskov/00_git_repos/PROJ_IrOx_Active_Learning_OER/run_nersc_vasp/ml_bulk_opt/run_all_bulks \
. \
--exclude "__temp__/**" \
--exclude "__old__/**" \
--exclude "out_data/**"
