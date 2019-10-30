import os
import sys

proj_data_save_dir = os.path.join(
    "01_norskov/PROJECT_DATA",
    "04_IrOx_surfaces_OER/ml_bulk_irox_dft/iro2")

ignore_ids = [
    "088",
    "156",
    "164",
    "174",
    "242",
    "255",
    "370",
    ]


input_dict = {
    "proj_data_save_dir": proj_data_save_dir,
    "root_dir": "/global/cscratch1/sd/flores12/IrOx_Project_temp_190510/ml_bulk_calculations/iro2_calcs",
    "ignore_ids": ignore_ids,
    }
