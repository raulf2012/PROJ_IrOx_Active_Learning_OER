# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.7
#   kernelspec:
#     display_name: Python [conda env:research-new]
#     language: python
#     name: conda-env-research-new-py
# ---

# # Import Modules

# + {"jupyter": {"source_hidden": true}}
# %%capture

import os
import sys

import copy

import time
import pickle

import numpy as np
import pandas as pd

import gpflow

from sklearn.decomposition import PCA

# #############################################################################
from catlearn.regression.gaussian_process import GaussianProcess
from catlearn.preprocess.clean_data import (
    clean_infinite,
    clean_variance,
    clean_skewness)
from catlearn.preprocess.scaling import standardize

# #############################################################################
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "python_classes/active_learning"))
from active_learning import (
    ALBulkOpt,
    ALGeneration,
    RegressionModel,
    FingerPrints,
    CandidateSpace,
    )

from al_analysis import ALAnalysis, ALAnimation

# #############################################################################
from IPython.display import display
# -

# # Script Inputs

# +
stoich_i = "AB2"
verbose = False
num_gen_stop = 5
# num_gen_stop = 10

gp_settings = {
    "noise": 0.02542,
    "sigma_l": 1.0049,
    "sigma_f": 5.19,
    "alpha": 0.018,
    }
# -

# # Read Data

# +
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling"))
from ml_methods import get_data_for_al

out_dict = get_data_for_al(
    stoich="AB2", verbose=False,
    drop_too_many_atoms=True,
#     drop_too_many_atoms=False,
    )

df_bulk_dft = out_dict["df_bulk_dft"]
df_bulk_dft = df_bulk_dft[df_bulk_dft["source"] == "raul"]

# +
# df_bulk_dft = df_bulk_dft[["atoms", "energy_pa"]]
df_bulk_dft = df_bulk_dft[["atoms", "dH"]]
df_bulk_dft.columns.values[1] = "y_real"

df_features_pre = out_dict["df_features_pre"]
df_features_post = out_dict["df_features_post"]

df_ids = out_dict["df_ids"]


df_static_irox = out_dict["df_static_irox"]
df_dij = out_dict["df_dij"]

# +
ids_w_dft = df_bulk_dft.index

# TEMP | Reduce size of candidate space
np.random.seed(8)
# ids_w_dft = np.sort(np.random.choice(np.sort(ids_w_dft), size=200))
ids_w_dft = list(set(ids_w_dft))

df_bulk_dft = df_bulk_dft.loc[ids_w_dft]

df_features_pre = df_features_pre.loc[ids_w_dft]
df_features_post = df_features_post.loc[ids_w_dft]

# + {"jupyter": {"source_hidden": true}}
color_list = [
    "rgb(202,88,66)",
    "rgb(71,189,198)",
    "rgb(210,70,147)",
    "rgb(120,181,66)",
    "rgb(157,99,201)",
    "rgb(81,163,108)",
    "rgb(189,104,138)",
    "rgb(131,128,57)",
    "rgb(101,130,203)",
    "rgb(209,154,68)",
    ]

ids_top_ten = [
    '64cg6j9any',
    'n36axdbw65',
    'clc2b1mavs',
    'ck638t75z3',
    'mkbj6e6e9p',
    'b49kx4c19q',
    '85z4msnl6o',
    'bpc2nk6qz1',
    '926dnunrxf',
    'mwmg9p7s6o',

    # "6r716sxr9t",
    # "n36axdbw65",
    # "clc2b1mavs",
    # "ck638t75z3",
    # "mkbj6e6e9p",
    # "vp7fvs6q81",
    # "85z4msnl6o",
    # "bpc2nk6qz1",
    # "926dnunrxf",
    # "mwmg9p7s6o",
    ]

id_color_dict = dict(zip(
    ids_top_ten,
    # df_bulk_dft.sort_values("y_real").iloc[0:10].index,
    color_list,
    ))


ids_top_system_duplicates = [
"6r716sxr9t",
# "6avov5cy64",
# "cfcivdxrc2",
# "m29j648g6i",
# "vunhmsbrml",
"9yz2mt8hbh",
# "nazu9q9l9h",
# "64cg6j9any",
# "b46enqnq8e",
]

# "6r716sxr9t",
# "9yz2mt8hbh",
# "64cg6j9any"

color_list = [
"green",
"blue",
# "rgb(202,88,66)",
# "rgb(71,189,198)",
# "rgb(210,70,147)",
# "rgb(120,181,66)",
# "rgb(157,99,201)",
# "rgb(81,163,108)",
# "rgb(189,104,138)",
# "rgb(131,128,57)",
# "rgb(101,130,203)",
]

id_color_dict = dict(zip(
    ids_top_system_duplicates,
    color_list,
    ))

print(id_color_dict)

# + {"active": ""}
#
#
#
#
#
#
#
# -

# # CCF Class

# + {"jupyter": {"source_hidden": true}}
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "python_classes"))
from ccf_similarity.ccf import CCF

d_thresh = 0.02
CCF = CCF(
    df_dij=df_dij,
    d_thresh=d_thresh)

# +
RM = RegressionModel(
    opt_hyperparameters=True,
    gp_settings_dict=gp_settings,
    verbose=verbose,
    )

FP = FingerPrints(
    df_features_pre,
    df_features_post=df_features_post,
    pca_mode="num_comp",  # 'num_comp' or 'perc'
    pca_comp=11,
    pca_perc=None,
    verbose=verbose,
    )

CS = CandidateSpace(
    Y_data=df_bulk_dft,
    Y_key="y_real",
    FingerPrints=FP,
    )

# +
name_i = "AL_" + stoich_i + "_" + str(num_gen_stop).zfill(2)
print("name:", name_i, "\n")
AL = ALBulkOpt(
    CandidateSpace=CS,
    RegressionModel=RM,
    DuplicateFinder=CCF,  # Optional
    num_seed_calcs=11,
    acquisition_bin=50,
    stop_mode="num_generations",
#     stop_mode=None,
    stop_num_generations=num_gen_stop,
    name="TEST__acq_10",
    verbose=verbose,
    )

run_al = True
if run_al:
    AL.run_AL()

    AL.__save_state__()

# + {"active": ""}
#
#
#

# +
# #############################################################################
import pickle; import os
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow",
    "dev_new_al_class/out_data",
    "TEST__acq_10.pickle")
    # "TEST.pickle")
    # "TEST_small.pickle")

with open(path_i, "rb") as fle:
    AL = pickle.load(fle)

ALAnim = ALAnimation(ALBulkOpt=AL, verbose=True)

if True:
    ALAnim.create_animation(
        # duration_long=1000 * 0.5,
        # duration_short=800 * 0.5,
        duration_long=1000 * 3,
        duration_short=800 * 3,
        serial_parallel="parallel",  # 'serial' or 'parallel'
        marker_color_dict=id_color_dict,
        )

# +
# assert False

# +
# #############################################################################
self = AL

duplicate_ids = AL.duplicate_ids

duplicate_swap_dict = AL.duplicate_swap_dict
# #############################################################################

# +
"64cg6j9any" in list(set(AL.duplicate_ids))

# AL.duplicate_ids

# +
# #############################################################################
gen_i = 1

AL_i = AL.al_gen_dict[gen_i]
self = AL_i

model = self.model
# #############################################################################

# AL_i.duplicate_swap_dict

model.loc[
    ["64cg6j9any", "9yz2mt8hbh", "6r716sxr9t"]]

# model

"6r716sxr9t" in AL.duplicate_ids
# -



# +
# indices_that_are_duplicates
# indices_that_are_duplicates
# .extend(

# indices_that_are_duplicates_i

# + {"active": ""}
#
#
#
#

# + {"active": ""}
#               y_real         y       err  acquired  gen_acquired
# id_unique                                                       
# 6avov5cy64 -7.045677 -7.045334  0.000591      True           1.0
# cfcivdxrc2 -7.045449 -7.046304  0.000632      True           1.0
# m29j648g6i -7.045449 -7.046259  0.000622      True           1.0
# vunhmsbrml -7.045582 -7.045695  0.000579      True           1.0
# nazu9q9l9h -7.045528 -7.045367  0.000581      True           1.0
# b46enqnq8e -7.047508 -7.046924  0.000705      True           4.0
#
# 64cg6j9any -7.047516 -7.046930  0.000701      True           4.0
# 6r716sxr9t -7.040906 -7.041141  0.001500      True           0.0
# 9yz2mt8hbh -7.047426 -7.046936  0.000704      True           1.0
#
# y_real             -7.04752
# y                  -7.04693
# err             0.000701333
# acquired               True
# gen_acquired              4
# Name: 64cg6j9any, dtype: object
#
# 64cg6j9any
# b46enqnq8e

# + {"active": ""}
#
#
#
# -

[i for i in df_bulk_dft.index.tolist() if "64cg" in i]
