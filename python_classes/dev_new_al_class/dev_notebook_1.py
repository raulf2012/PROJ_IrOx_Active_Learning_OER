# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python [conda env:research-new]
#     language: python
#     name: conda-env-research-new-py
# ---

# # Import Modules

# + jupyter={}
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
    # ALGeneration,
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
# stoich_i = "AB2"
# verbose = False
# num_gen_stop = 5
# # num_gen_stop = 10

# gp_settings = {
#     "noise": 0.02542,
#     "sigma_l": 1.0049,
#     "sigma_f": 5.19,
#     "alpha": 0.018,
#     }


# # #############################################################################
# color_list = [
#     "rgb(202,88,66)",
#     "rgb(71,189,198)",
#     "rgb(210,70,147)",
#     "rgb(120,181,66)",
#     "rgb(157,99,201)",
#     "rgb(81,163,108)",
#     "rgb(189,104,138)",
#     "rgb(131,128,57)",
#     "rgb(101,130,203)",
#     "rgb(209,154,68)",
#     ]

# ids_top_ten = [
#     '64cg6j9any',
#     'n36axdbw65',
#     'clc2b1mavs',
#     'ck638t75z3',
#     'mkbj6e6e9p',
#     'b49kx4c19q',
#     '85z4msnl6o',
#     'bpc2nk6qz1',
#     '926dnunrxf',
#     'mwmg9p7s6o',
#     ]

# id_color_dict = dict(zip(
#     ids_top_ten,
#     color_list,
#     ))
# -

# # Read Data

# +
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling"))
from ml_methods import get_data_for_al

out_dict = get_data_for_al(
    stoich=stoich_i, verbose=False,
    drop_too_many_atoms=True,
#     drop_too_many_atoms=False,
    )

df_bulk_dft = out_dict["df_bulk_dft"]
df_bulk_dft = df_bulk_dft[df_bulk_dft["source"] == "raul"]

# df_bulk_dft = df_bulk_dft[["atoms", "energy_pa"]]
df_bulk_dft = df_bulk_dft[["atoms", "dH"]]
df_bulk_dft.columns.values[1] = "y_real"

df_features_pre = out_dict["df_features_pre"]
df_features_post = out_dict["df_features_post"]

df_ids = out_dict["df_ids"]


df_static_irox = out_dict["df_static_irox"]
df_dij = out_dict["df_dij"]
# -

# # Filter to candidates w/ DFT energy

# +
ids_w_dft = df_bulk_dft.index

# TEMP | Reduce size of candidate space
# np.random.seed(8)
# ids_w_dft = np.sort(np.random.choice(np.sort(ids_w_dft), size=200))
ids_w_dft = list(set(ids_w_dft))

df_bulk_dft = df_bulk_dft.loc[ids_w_dft]

df_features_pre = df_features_pre.loc[ids_w_dft]
df_features_post = df_features_post.loc[ids_w_dft]

# + active=""
#
#
#
#
# -

# # CCF Class

# +
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
    pca_comp=10,
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
    # num_seed_calcs=11,
    num_seed_calcs=5,
    acquisition_bin=10,
#     stop_mode="num_generations",
    stop_mode=None,
    stop_num_generations=num_gen_stop,
    name="TEST__0",
    verbose=verbose,
    acquisition_method="gp_ucb",
    )

run_al = True
if run_al:
    AL.run_AL()
    AL.duplicate_system_history_analysis()
    AL.__save_state__()

# +
# df = AL.CandidateSpace.FingerPrints.df_train
# num_data_points = df.shape[0]
# if num_data_points < pca_comp:
#     pca_comp = num_data_points + 1
# -

assert False

# + active=""
#
#
#

# +
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow",
    "dev_new_al_class/out_data",

    "TEST__acq_5_all.pickle"
    # "TEST__acq_10.pickle"
    )

with open(path_i, "rb") as fle:
    AL = pickle.load(fle)

ALAnim = ALAnimation(
    ALBulkOpt=AL,
    marker_color_dict=id_color_dict,
    verbose=True)

if False:
    ALAnim.create_animation(
        # duration_long=1000 * 0.5,
        # duration_short=800 * 0.5,
        duration_long=1000 * 4,
        duration_short=800 * 4,
        serial_parallel="parallel",  # 'serial' or 'parallel'
        # marker_color_dict=id_color_dict,
        )

# +
from al_analysis import ALPerformance

ALPerf = ALPerformance(
    ALBulkOpt=AL,
    verbose=False,
    )
ALPerf.num_sys_discovered()

df = ALPerf.num_sys_discovered_df

# df

# +
import chart_studio.plotly as py
import plotly.graph_objs as go
import os


trace = go.Scatter(
    x=df["num_dft"],
    y=df["num_ids_discovered"],
    mode="markers",
    marker=dict(
        symbol="circle",
        color="grey",
        size=14,
        line=dict(
            color='black',
            width=2
            )
        ),
    )

data = [trace]

fig = go.Figure(data=data)
fig.show()

# + active=""
#
#
#
# -

# # SANDBOX

# +
# #############################################################################
self = AL

al_gen_dict = AL.al_gen_dict
duplicate_ids = AL.duplicate_ids

duplicate_swap_dict = AL.duplicate_swap_dict
# #############################################################################
# -

CandidateSpace = AL.CandidateSpace
FingerPrints = CandidateSpace.FingerPrints
df_pre = FingerPrints.df_pre
all_indices = df_pre.index.tolist()

# +
# #############################################################################
gen_i = 0

AL_i = al_gen_dict[gen_i]
self = AL_i

model = self.model
prev_duplicate_ids = self.prev_duplicate_ids
indices_that_are_duplicates = self.indices_that_are_duplicates
# #############################################################################
# -

assert False

# # Performance # of structures vs DFT calcs

# + active=""
#
#
#

# +
# # Percent of total number of structures to track
# perc_of_structs = 10

# top_idslast_gen = list(al_gen_dict.keys())[-1]


# AL_last = al_gen_dict[last_gen]
# model = AL_last.model


# # #############################################################################
# num_candidates_init = model.shape[0]
# num_track_structs = round(num_candidates_init * (perc_of_structs * 0.01))



# model_tmp = model[model["duplicate"] == False]
# model_tmp = model_tmp.sort_values("y_real")
# model_tmp = model_tmp.iloc[0:num_track_structs]

# top_ids = model_tmp.index.tolist()

# top_ids_static = copy.deepcopy(top_ids)

# swap_histories = ALAnim.swap_histories

# duplicates_of_top_ids = []
# for id_i in top_ids:
#     if id_i in swap_histories.keys():
#         swap_lists = swap_histories.get(id_i, "TEMP")

#         swap_ids_i = []
#         for gen_j, swap_list_j in swap_lists.items():
#             swap_ids_i.extend(swap_list_j)

#         duplicates_of_top_ids.extend(swap_ids_i)

# duplicates_of_top_ids = list(set(duplicates_of_top_ids))

# print(len(duplicates_of_top_ids + top_ids))

# print(len(set(duplicates_of_top_ids + top_ids)))

# top_ids_w_dupl = list(set(duplicates_of_top_ids + top_ids))

# new_swap_dict = dict()
# for id_i, swap_history_i in swap_histories.items():
#     for gen_j, swap_list_j in swap_history_i.items():
#         for swap_id in swap_list_j:
#             # #################################################################
#             if swap_id in new_swap_dict.keys():
#                 if new_swap_dict[swap_id] != id_i:
#                     print("This id corresponds to more than 1 final id")

#             new_swap_dict[swap_id] = id_i
#             # #################################################################

# data_list_master = []
# for gen_i, AL_i in al_gen_dict.items():
#     data_dict_i = dict()

#     model_i = AL_i.model
#     model_tmp = model_i[
#         (model_i["acquired"] == True) & \
#         (model_i["duplicate"] == False)
#         ]

#     # Number of DFT experiments
#     num_dft_calcs = model_i[model_i["acquired"] == True].shape[0]
#     data_dict_i["num_dft"] = num_dft_calcs

#     for id_i in model_tmp.index:
#         if id_i in top_ids:
#             # print("id_i directly in top ids")
#             top_ids.remove(id_i)

#         final_swap_id = new_swap_dict.get(id_i, None)
#         if final_swap_id is not None and final_swap_id in top_ids:
#             # print("DIFJIDS")

#             top_ids.remove(final_swap_id)

#     num_ids_disc = len(top_ids_static) - len(top_ids)
#     data_dict_i["num_ids_discovered"] = num_ids_disc

#     # print(num_ids_disc)

#     data_list_master.append(data_dict_i)

# df = pd.DataFrame(data_list_master)

# df

# +
# ALAnim.swap_histories


# ALAnim.__get_color_dict__(gen_i=0)

# al_gen_dict

# ALAnim.swap_histories
