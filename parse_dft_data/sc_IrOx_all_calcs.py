#!/usr/bin/env python

"""IrOx Project.

Author: Raul A. Flores
"""

# | - Import Modules
from datetime import datetime
startTime = datetime.now()

import os
# from os.path import join as join

import sys
sys.path.append(".")

from methods import job_maint
from dft_job_automat.job_analysis import DFT_Jobs_Analysis
from dft_job_automat.job_types_classes.dft_methods import DFT_Methods
#__|

# | - Script Parameters
parse_data = True
parse_all_rev = False

maint_data = True
cross_check_jobs = True

parallel_exec = False
#__|

from job_dirs import dir_list

if parse_data:
    # | - Instantiate Classes
    dft_inst = DFT_Methods(
        methods_to_run=[
            # "elec_energy",
            # "init_atoms",
            # "atoms_object",
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
