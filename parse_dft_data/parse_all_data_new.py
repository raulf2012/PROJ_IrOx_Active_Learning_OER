# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_IrOx_Active_Learning_OER]
#     language: python
#     name: conda-env-PROJ_IrOx_Active_Learning_OER-py
# ---

# # Parse extra data
# ---
#
# I'm coming back to this after a whle
#
# There was only really 1 extra row that I needed to get into the dataframe for processing
#
# I'll do that here and then read it in here:
#
# PROJ_IrOx_Active_Learning_OER/workflow

read_from_PROJ_DATA = False
read_from_PROJ_DATA = True

# +
# | - Import Modules
import os
print(os.getcwd())
import sys

sys.path.append(".")

import numpy as np
import pandas as pd

from datetime import datetime
startTime = datetime.now()

from methods import job_maint
from dft_job_automat.job_analysis import DFT_Jobs_Analysis
from dft_job_automat.job_types_classes.dft_methods import DFT_Methods
#__|
# -

# #############################################################################
# Project Data
# from proj_data_irox import data_dir


# + active=""
#
#
#

# +
# # #########################################################
# import pickle; import os
# path_i = os.path.join(
#     data_dir,
#     # "df_master.pickle",
#     "job_dataframe.pickle"
#     )

# with open(path_i, "rb") as fle:
#     df_prev = pickle.load(fle, encoding="latin1")
# # #########################################################

# df_prev.columns.tolist()
# -

['Job',
 'adsorbate',
 'bulk_system',
 'coverage',
 'coverage_type',
 'facet',
 'job_type',
 'layers',
 'max_revision',
 'ooh_direction',
 'path',
 'priority',
 'revision_number',
 'success',
 'surface_type',
 'elec_energy',
 'init_atoms',
 'atoms_object',
 'incar']

# +
# | - Script Parameters
parse_data = True
parse_all_rev = False

maint_data = True
cross_check_jobs = True

parallel_exec = False
#__|

from job_dirs import dir_list

# +
dir_list_new = [
    "norskov_research_storage/nersc/IrOx_Project_temp_190510/03_OER_Calc/IrO2/100/01_O_covered/02_ooh/02_face_down_2"
    ]

dir_list = []
for dir_i in dir_list_new:
    dir_new = os.path.join(
        os.environ["gdrive"],
        dir_i)

    dir_list.append(dir_new)
# -

cwd_orig = os.getcwd()

if read_from_PROJ_DATA:
    import pickle; import os
    path_i = os.path.join(
        os.environ["PROJ_DATA"],
        "04_IrOx_surfaces_OER/PROJECT_COMPUTED_OUT_DATA/PROJ_IrOx_Active_Learning_OER",
        "parse_dft_data",
        "out_data/df_data_new.pickle")
    with open(path_i, "rb") as fle:
        df_m = pickle.load(fle)

# +
# assert False
# -

if parse_data and not read_from_PROJ_DATA:
    # | - Instantiate Classes
    dft_inst = DFT_Methods(
        methods_to_run=[
            "elec_energy",
            "init_atoms",
            "atoms_object",
            # "incar",
            # "outcar"
            ],
        DFT_code="VASP",
        )

    Jobs = DFT_Jobs_Analysis(
        indiv_dir_lst=dir_list,
        working_dir=".",
        folders_exist=True,
        load_dataframe=False,
        job_type_class=dft_inst,
        parse_all_revisions=parse_all_rev,
        parallel_exec=parallel_exec,
        )

    df_all = Jobs.data_frame
    df_m = Jobs.filter_early_revisions(Jobs.data_frame)
    #__|
os.chdir(cwd_orig)

# Pickling data ###########################################
import os; import pickle
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)
with open(os.path.join(directory, "df_data_new.pickle"), "wb") as fle:
    pickle.dump(df_m, fle)
# #########################################################

# +
# # #########################################################
# import pickle; import os
# path_i = os.path.join(
#     os.environ["PROJ_irox"],
#     "parse_dft_data/out_data",
#     "df_data_new.pickle")
# with open(path_i, "rb") as fle:
#     df_m = pickle.load(fle)
# # #########################################################

# pd.concat([
# df_m, df_m   
    
# ], axis=0)
# -

print(20 * "# # ")
print("All done!")
assert False

# # Read Previous DataFrame and Combine Data

# + active=""
#
#
#
#
#
#
# -

if maint_data:
    # | - Job Maintance
    print(25 * "*")

    tally = {"successes": 0, "failures": 0, "running": 0, "pending": 0}

    for Job_i in Jobs.Job_list:
        path_i = Job_i.full_path
        job_i_params = Job_i.job_params

        print(path_i)

        tally = job_maint(
            0,
            path_i,
            job_i_params,
            {"jobs_man_list": [Jobs]},
            tally,
            file_ops=False,
            )
    #__|

# +
if cross_check_jobs:
    # | - NEW | Parse for Job Folders w/o dir_list
    from dft_job_automat.job_analysis import (
        parse_job_dirs,
        compare_parsed_and_user_job_dirs,
        )

    rt_1 = os.path.join(os.environ["wd"], "IrOx_Project")
    dirs_to_parse = [
        os.path.join(rt_1, "01_surface_calcs"),
        os.path.join(rt_1, "02_surface_coverage"),
        os.path.join(rt_1, "03_OER_Calc"),
        os.path.join(rt_1, "07_diff_coverages_term"),
        ]

    parsed_dir_list = parse_job_dirs(dirs_to_parse)
    compare_parsed_and_user_job_dirs(parsed_dir_list, dir_list)

    for path_i in parsed_dir_list:
        files_i = os.listdir(path_i)
    #__|

print(datetime.now() - startTime)
