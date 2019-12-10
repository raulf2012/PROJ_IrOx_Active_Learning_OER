# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.0
#   kernelspec:
#     display_name: Python [conda env:PROJ_IrOx_Active_Learning_OER]
#     language: python
#     name: conda-env-PROJ_IrOx_Active_Learning_OER-py
# ---

# + {"Collapsed": "false"}
import os

import time

import numpy as np

import itertools

from methods import run_al_i

from multiprocessing import Pool
from functools import partial

# + {"Collapsed": "false"}
t0 = time.time()

# + {"Collapsed": "false"}
stoich_i = "AB2"
verbose = False

gp_settings = {
    "noise": 0.02542,
    "sigma_l": 1.0049,
    "sigma_f": 5.19,
    "alpha": 0.018,
    }

# duplicate_analysis = False
# acquisition_method = "gp_ucb"
# acquisition_method = "random"

# #############################################################################
runs_list = list(range(5))
acquisition_methods = ["gp_ucb", "random"]
duplicate_analysis = [True, False]

# TEST SETTINGS # #############################################################
runs_list = list(range(1))
acquisition_methods = ["gp_ucb"]
duplicate_analysis = [True]

# + {"Collapsed": "false"}
from misc_modules.misc_methods import GetFriendlyID

# GetFriendlyID()

# + {"Collapsed": "false"}
variables_dict = dict(
    stoich_i=stoich_i,
    # acquisition_method=acquisition_method,
    # duplicate_analysis=duplicate_analysis,
    verbose=verbose,
    gp_settings=gp_settings,
    )

def run_al_meth(
    input_dict,
    stoich_i=None,
    # acquisition_method=None,
    verbose=None,
    gp_settings=None,
    name_i=None,
    save_dir=None,
    # duplicate_analysis=None,
    ):
    """
    """
    i = input_dict['i']
    acquisition_method = input_dict['acquisition_method']
    duplicate_analysis = input_dict['duplicate_analysis' ]
    seed = input_dict['seed']

    # #########################################################################
    i_str = str(i).zfill(2)

    
    print(80 * "#")
    print(i_str, 77 * "#")

    name_i = "TEST_AL_" + GetFriendlyID()
    save_dir = stoich_i + "/" + acquisition_method + "_" + str(duplicate_analysis)
    run_al_i(
        stoich_i=stoich_i,
        verbose=verbose,
        gp_settings=gp_settings,
        name_i=name_i,
        save_dir_extra=save_dir,
        acquisition_method=acquisition_method,
        duplicate_analysis=duplicate_analysis,
        seed=seed,
        )

# + {"Collapsed": "false"}
# for i in range(num_runs):
#     run_al_meth(i, **variables_dict)

# + {"Collapsed": "false"}
input_list = []
for i in itertools.product(acquisition_methods, duplicate_analysis, runs_list):
    data_dict_i = dict(
        acquisition_method=i[0],
        duplicate_analysis=i[1],
        i=i[2],
        seed=np.random.randint(0, 1000)
        )
    input_list.append(data_dict_i)
    
# input_list

# + {"Collapsed": "false"}
traces_all = Pool().map(
    partial(
        run_al_meth,  # METHOD
        **variables_dict,  # KWARGS
        ),
    input_list,
    )

# + {"Collapsed": "false"}
print("Notebook runtime (s):", time.time() - t0)
print("Notebook runtime (s):", (time.time() - t0) / 60)

# + {"Collapsed": "false", "active": ""}
#
#
#

# + {"jupyter": {"source_hidden": true}, "Collapsed": "false"}
# input_list = []
# for i in range(num_runs):
#     data_dict_i = dict(
#         i=i,
#         seed=np.random.randint(0, 1000)
#         )
#     print(data_dict_i)
#     input_list.append(data_dict_i)
# # input_list

# + {"jupyter": {"source_hidden": true}, "Collapsed": "false"}
# for i in range(1):
#     i_str = str(i).zfill(2)
#     print(i_str)

#     print(80 * "#")
#     print(i_str, 77 * "#")

#     name_i = "AL_" + i_str
#     save_dir = stoich_i + "/" + acquisition_method + "_" + str(duplicate_analysis)
#     # print(save_dir)
#     run_al_i(
#         stoich_i=stoich_i,
#         verbose=verbose,
#         gp_settings=gp_settings,
#         name_i=name_i,
#         save_dir_extra=save_dir,
#         acquisition_method=acquisition_method,
#         duplicate_analysis=duplicate_analysis,
#         )
#     print("")

# + {"jupyter": {"source_hidden": true}, "Collapsed": "false"}
# AL_i_list = [i for i in ALBulkOpt.al_gen_dict.values()]

# traces_all = Pool().map(partial(
#     get_trace_j,  # METHOD

#     # KWARGS
#     **get_trace_j_kwargs,
#     ), AL_i_list)
