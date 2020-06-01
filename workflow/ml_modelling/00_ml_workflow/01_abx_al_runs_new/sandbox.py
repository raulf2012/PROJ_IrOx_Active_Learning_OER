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

# + jupyter={} endofcell="--"
# | - Import Modules
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
from active_learning.al_bulkopt import ALBulkOpt
from active_learning.active_learning import (
    RegressionModel,
    FingerPrints,
    CandidateSpace,
    )
from active_learning.al_analysis import ALAnalysis, ALAnimation

# #############################################################################
from IPython.display import display
# -

# from inputs import (
#     stoich_i,
#     verbose,
#     gp_settings,
#     name_i,
#     )
#__|
# --

stoich_i = "AB3"
verbose = True

# +
# def run_al_i(
#     stoich_i=None,
#     verbose=None,
#     gp_settings=None,
#     name_i=None,
#     save_dir_extra=None,
#     acquisition_method=None,
#     duplicate_analysis=None,
#     seed=None,
#     ):
#     """
#     """

# | - run_al_i
# # Read Data

# + jupyter={}
# # # + {"jupyter": {"source_hidden": true}}
# sys.path.insert(0, os.path.join(
#     os.environ["PROJ_irox"],
#     "workflow/ml_modelling"))
# from ml_methods import get_data_for_al

# out_dict = get_data_for_al(
#     stoich=stoich_i, verbose=False,
#     drop_too_many_atoms=True,
# #     drop_too_many_atoms=False,
#     )

# df_bulk_dft = out_dict["df_bulk_dft"]
# df_bulk_dft = df_bulk_dft[df_bulk_dft["source"] == "raul"]

# # df_bulk_dft = df_bulk_dft[["atoms", "energy_pa"]]
# df_bulk_dft = df_bulk_dft[["atoms", "dH"]]
# df_bulk_dft.columns.values[1] = "y_real"

# df_features_pre = out_dict["df_features_pre"]
# df_features_post = out_dict["df_features_post"]

# df_ids = out_dict["df_ids"]


# df_static_irox = out_dict["df_static_irox"]
# df_dij = out_dict["df_dij"]
# # -

# # # Filter to candidates w/ DFT energy

# # # + {"jupyter": {"source_hidden": true}}
# ids_w_dft = df_bulk_dft.index

# # TEMP | Reduce size of candidate space
# # np.random.seed(8)
# # ids_w_dft = np.sort(np.random.choice(np.sort(ids_w_dft), size=200))
# ids_w_dft = list(set(ids_w_dft))
# # print("ids_w_dft:", ids_w_dft)

# df_bulk_dft = df_bulk_dft.loc[ids_w_dft]

# df_features_pre = df_features_pre.loc[ids_w_dft]
# df_features_post = df_features_post.loc[ids_w_dft]

# # # + {"active": ""}
# #
# #
# #
# #
# # -

# # # CCF Class

# # # # + {"jupyter": {"source_hidden": true}}
# # sys.path.insert(0, os.path.join(
# #     os.environ["PROJ_irox"],
# #     "python_classes"))
# # from ccf_similarity.ccf import CCF

# # d_thresh = 0.02
# # CCF = CCF(
# #     df_dij=df_dij,
# #     d_thresh=d_thresh)

# # # # + {"jupyter": {"source_hidden": true}}
# # RM = RegressionModel(
# #     opt_hyperparameters=True,
# #     gp_settings_dict=gp_settings,
# #     verbose=verbose,
# #     )

# # FP = FingerPrints(
# #     df_features_pre,
# #     df_features_post=df_features_post,
# #     pca_mode="num_comp",  # 'num_comp' or 'perc'
# #     pca_comp=10,
# #     pca_perc=None,
# #     verbose=verbose,
# #     )

# # CS = CandidateSpace(
# #     Y_data=df_bulk_dft,
# #     Y_key="y_real",
# #     FingerPrints=FP,
# #     )

# # # # +
# # AL = ALBulkOpt(
# #     CandidateSpace=CS,
# #     RegressionModel=RM,
# #     DuplicateFinder=CCF,  # Optional
# #     duplicate_analysis=duplicate_analysis,
# #     # num_seed_calcs=11,
# #     num_seed_calcs=5,
# #     acquisition_bin=5,
# #     # stop_mode="num_generations",
# #     stop_mode=None,
# #     stop_num_generations=3,
# #     name=name_i,
# #     save_dir_extra=save_dir_extra,
# #     verbose=verbose,
# #     # acquisition_method="gp_ucb",
# #     acquisition_method=acquisition_method,
# #     seed=seed,
# #     )

# # run_al = True
# # if run_al:
# #     AL.run_AL()
# #     AL.duplicate_system_history_analysis()
# #     AL.__save_state__()

# #__|
# -

sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling"))
from ml_methods import get_data_for_al
from ml_methods import get_ml_dataframes

stoich = stoich_i
verbose = verbose
drop_too_many_atoms = True

# +
# def get_data_for_al(
# stoich="AB2",
# verbose=True,
# drop_too_many_atoms=True,
# ):
# """
# """

# | - get_data_for_al
sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import (ids_to_discard__too_many_atoms_path)

# | - Get all necessary dfs
df_dict = get_ml_dataframes(
    names=[
        "bulk_dft_data_path",
        "unique_ids_path",
        # "prototypes_data_path",
        "static_irox_structures_path",
        # "static_irox_structures_kirsten_path",
        # "oqmd_irox_data_path",
        "df_features_pre_opt_path",
        "df_features_pre_opt_kirsten_path",
        "df_features_post_opt_path",
        # "oer_bulk_structures_path",
        # "df_ccf_path",
        "df_dij_path",
        # "ids_to_discard__too_many_atoms_path",
        ],
    )

df_ids = df_dict.get("unique_ids", None)
df_bulk_dft = df_dict.get("bulk_dft_data", None)
df_features_pre = df_dict.get("df_features_pre_opt", None)
# df_features_pre = df_dict.get("df_features_pre_opt_kirsten", None)
df_features_post = df_dict.get("df_features_post_opt", None)

df_dij = df_dict.get("df_dij", None)

print("ISDFIODISFIDS*F*SDF*SDYUGFSODIUFG")
print("6fcdbh9fz2 in df_bulk_dft", "6fcdbh9fz2" in df_bulk_dft.index)


df_static_irox = df_dict.get("static_irox_structures", None)
#__|

# | - Filter ids to user specifications
df_ids = df_ids[
    (df_ids["stoich"] == stoich) & \
    (df_ids["source"] != "oqmd") & \
    (df_ids["source"] != "raul_oer") & \
    [True for i in range(len(df_ids))]]
ids = df_ids["unique_ids"]
#__|

# #####################################################

# | - DFT dataframe
df_i = df_bulk_dft

# print("isidfjisdjifjsidjf8yu2894h90832uy4908tyu98023wht0982quj098gtfujw3e")
# print(df_i.index.shape)
# print(df_i.index.unique().shape)

# Common ids between user ids and df
common_ids = list(set(df_i.index) & set(ids))

ids_not_in__df_i = [i for i in ids if i not in common_ids]

df_i = df_i.loc[common_ids]

if verbose:
    print("len(ids):", len(ids))
    print("len(common_ids)", len(common_ids))
    print("len(ids_not_in__bulk_dft_data):", len(ids_not_in__df_i))
    print("\n", "df_i.shape: ", df_i.shape, sep="")

df_i = df_i[df_i.source == "raul"]

df_bulk_dft = df_i

print("ISDFIODISFIDS*F*SDF*SDYUGFSODIUFG")
print("6fcdbh9fz2 in df_bulk_dft", "6fcdbh9fz2" in df_bulk_dft.index)

# print("TEMP TEMP TEMP 89ihsjdgf", "6dzhcimdxs" in df_bulk_dft.index)
#__|

# | - Featurs pre-DFT
df_i = df_features_pre

# Common ids between user ids and df
common_ids = list(set(df_i.index) & set(ids))

ids_not_in__df_i = [i for i in ids if i not in common_ids]

df_i = df_i.loc[common_ids]

if verbose:
    print("len(ids):", len(ids))
    print("len(common_ids)", len(common_ids))
    print("len(ids_not_in__bulk_dft_data):", len(ids_not_in__df_i))
    print("\n", "df_i.shape: ", df_i.shape, sep="")

df_features_pre = df_i

df_features_pre = df_features_pre["voronoi"]
#__|

# | - Features post-DFT
df_i = df_features_post

# Common ids between user ids and df
common_ids = list(set(df_i.index) & set(ids))

ids_not_in__df_i = [i for i in ids if i not in common_ids]

df_i = df_i.loc[common_ids]

if verbose:
    print("len(ids):", len(ids))
    print("len(common_ids)", len(common_ids))
    print("len(ids_not_in__bulk_dft_data):", len(ids_not_in__df_i))
    print("\n", "df_i.shape: ", df_i.shape, sep="")

df_features_post = df_i

# Only use post-DFT features from my data set
df_features_post = \
    df_features_post[df_features_post["data"]["source"] == "raul"]

df_features_post = df_features_post["voronoi"]
#__|


# | - Dropping certain rows
all_ids = list(set(
    df_bulk_dft.index.tolist() + \
    df_features_pre.index.tolist() + \
    df_features_post.index.tolist() ))


ids_to_drop = []

# #########################################################################
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/processing_bulk_dft/static_prototypes_structures/out_data",
    "ids_to_discard__proto_dupl.pickle")
with open(path_i, "rb") as fle:
    ids_to_discard__proto_dupl = pickle.load(fle)
    # ids_to_drop.extend(ids_to_discard__proto_dupl)
# #########################################################################

if drop_too_many_atoms:
    # #####################################################################
    with open(ids_to_discard__too_many_atoms_path, "rb") as fle:
        ids_to_drop__too_many_atoms = pickle.load(fle)
        ids_to_drop.extend(ids_to_drop__too_many_atoms)

        # ids_to_drop = ids_to_drop__too_many_atoms
        # ids_to_drop = [i for i in ids_to_drop if i in all_ids]
    # #####################################################################


ids_to_drop = [i for i in ids_to_drop if i in all_ids]

# print("in ids to drop", "6fcdbh9fz2" in ids_to_drop)

df_features_pre = df_features_pre.drop(
    labels=ids_to_drop, axis=0)

df_i = df_features_post
df_features_post = df_i.loc[
    df_i.index.intersection(
        df_features_pre.index.unique()
        ).unique()
    ]


# print("TEMP TEMP TEMP 89ihsjdgf", "6dzhcimdxs" in df_bulk_dft.index)

df_bulk_dft = df_bulk_dft.loc[
    df_bulk_dft.index.intersection(
        df_features_pre.index
        ).unique()
    ]

print("ISDFIODISFIDS*F*SDF*SDYUGFSODIUFG")
print("6fcdbh9fz2 in df_bulk_dft", "6fcdbh9fz2" in df_bulk_dft.index)

# print("TEMP TEMP TEMP 89ihsjdgf", "6dzhcimdxs" in df_bulk_dft.index)

df_static_irox = df_static_irox.loc[
    df_static_irox.index.intersection(
        df_features_pre.index
        ).unique()
    ]

ids_static = df_dij.index.intersection(df_static_irox["static_id"])
ids_completed_post_dft = \
    df_dij.index.intersection(df_features_pre.index)


# print("TEMP TEMP TEMP", "6fcdbh9fz2" in df_dij.index)

ids_dij = ids_static.tolist() + ids_completed_post_dft.tolist()
df_dij = df_dij.loc[ids_dij, ids_dij]

# print("TEMP TEMP TEMP", "6fcdbh9fz2" in df_dij.index)

#__|

out_dict = dict()

out_dict["df_features_post"] = df_features_post
out_dict["df_features_pre"] = df_features_pre
out_dict["df_bulk_dft"] = df_bulk_dft

# TEMP
out_dict["df_ids"] = df_ids
out_dict["df_dij"] = df_dij
out_dict["df_static_irox"] = df_static_irox


# return(out_dict)
#__|

# +
# print("ISDFIODISFIDS*F*SDF*SDYUGFSODIUFG")
# print("6fcdbh9fz2 in df_bulk_dft", "6fcdbh9fz2" in df_bulk_dft.index)

df_bulk_dft.loc["6fcdbh9fz2"]

# + active=""
#
#
#
#
#
#
#
# -

stoich = "AB3"
verbose = True
drop_too_many_atoms = True

# +
# # + {"jupyter": {"source_hidden": true}}
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

# +
"6fcdbh9fz2" in df_bulk_dft.index

# df_bulk_dft

# + jupyter={}
# df_bulk_dft = df_bulk_dft[df_bulk_dft["source"] == "raul"]

# # df_bulk_dft = df_bulk_dft[["atoms", "energy_pa"]]
# df_bulk_dft = df_bulk_dft[["atoms", "dH"]]
# df_bulk_dft.columns.values[1] = "y_real"

# df_features_pre = out_dict["df_features_pre"]
# df_features_post = out_dict["df_features_post"]

# df_ids = out_dict["df_ids"]


# df_static_irox = out_dict["df_static_irox"]
# df_dij = out_dict["df_dij"]
# # -

# # # Filter to candidates w/ DFT energy

# # # + {"jupyter": {"source_hidden": true}}
# ids_w_dft = df_bulk_dft.index

# # TEMP | Reduce size of candidate space
# # np.random.seed(8)
# # ids_w_dft = np.sort(np.random.choice(np.sort(ids_w_dft), size=200))
# ids_w_dft = list(set(ids_w_dft))
# # print("ids_w_dft:", ids_w_dft)

# df_bulk_dft = df_bulk_dft.loc[ids_w_dft]

# df_features_pre = df_features_pre.loc[ids_w_dft]
# df_features_post = df_features_post.loc[ids_w_dft]
# -

assert False

df_bulk_dft.head()

df_features_pre

df_features_post
