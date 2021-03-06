# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_irox]
#     language: python
#     name: conda-env-PROJ_irox-py
# ---

# # New ML Active Learning Workflow
# ---
#
# A model that predicts the mean (~ -6.05 eV/atom) has a MAE of ~0.3 eV/atom)

# # Import Modules

# +
import os
print(os.getcwd())
import sys

import pickle

import numpy as np
import pandas as pd

# #########################################################
# Python Utils
import itertools
import time

# #########################################################
# Project Imports
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow/190611_new_workflow/02_gaus_proc"))
# -

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import read_from_PROJ_DATA

# +
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "python_classes"))

from active_learning.al_bulkopt import ALBulkOpt
from active_learning.active_learning import (
    RegressionModel,
    FingerPrints,
    CandidateSpace,
    )
# -

# # Script Inputs

# + jupyter={}
# stoich_i = "AB2"
stoich_i = "AB3"

# gp_model = gp_model_gpflow
# gp_model = gp_model_catlearn

aqs_bin_size = 5

# output_key = "form_e_chris"
output_key = "energy_pa"

verbosity_level = 6  # 1-10 scale
verbose = True

params_dict = {
    "noise": [0.0001],
    "sigma_l": [10.],
    "sigma_f": [5],
    "alpha": [0.1],
    }

c = list(itertools.product(*params_dict.values()))
df_gp_params = pd.DataFrame(c, columns=params_dict.keys())

gp_settings = df_gp_params.iloc[0].to_dict()

gp_settings = {
    "noise": 0.02542,
    "sigma_l": 1.0049,
    "sigma_f": 5.19,
    "alpha": 0.018,
    }
# -

# n_fold_cv = 15
n_fold_cv = 20
pca_comp = 11

# # Reading Data

# +
sys.path.insert(0,
    os.path.join(os.environ["PROJ_irox"], "workflow/ml_modelling"))
from ml_methods import get_ml_dataframes

DF_dict = get_ml_dataframes()

df_dft = DF_dict["df_dft_final_final"]
df_feat_pre = DF_dict["df_features_pre_opt"]
df_feat_post = DF_dict["df_features_post_opt"]

df_ids = DF_dict['unique_ids']
# -

# # Preprocessing Dataframes

# +
df_dft = df_dft[df_dft.stoich == stoich_i]

df_feat_post = df_feat_post[df_feat_post.data.source == "raul"]
df_feat_post = df_feat_post.drop(columns=["data"])

# #########################################################
df_feat_post = df_feat_post.loc[df_dft.index]
df_feat_pre = df_feat_pre.loc[df_dft.index]

# #########################################################
df_feat_post = df_feat_post["voronoi"]
df_feat_pre = df_feat_pre["voronoi"]

# + active=""
#
#
#
#
#
#
#
#
#
#
#
#
#
# -

if read_from_PROJ_DATA:
    path_i = os.path.join(
        os.environ["PROJ_DATA"], "04_IrOx_surfaces_OER",
        "PROJECT_COMPUTED_OUT_DATA/PROJ_IrOx_Active_Learning_OER",
        "workflow/ml_modelling/opt_mae_err_gp_model",
        "out_data/" + stoich_i + "_data.pickle")
    with open(path_i, "rb") as fle:
        df_m = pickle.load(fle)

if not read_from_PROJ_DATA:
    data_dict_list = []
    for pca_comp_i in range(1, 40, 1):
    # for pca_comp_i in range(1, 4, 1):
        print("pca_comp_i:", pca_comp_i)
        data_dict_j = dict()

        data_dict_j["pca_comp"] = pca_comp_i

        # #####################################################
        fold_size = int(df_dft.shape[0] / n_fold_cv)
        # Shuffling training data
        df_dft = df_dft.sample(
            n=None,
            frac=1.,
            replace=False,
            axis=None)
        # print("n_fold_cv * fold_size:", n_fold_cv * fold_size)
        ids_0 = df_dft.index[:n_fold_cv * fold_size]
        folds = np.split(ids_0, n_fold_cv)
        ids_leftover = df_dft.index[n_fold_cv * fold_size:]
        if ids_leftover.shape[0] > 0:
            folds.append(ids_leftover)
        folds = np.array(folds)



        # #####################################################
        data_dict_list_j = []
        for i_cnt, fold_i in enumerate(folds):
            data_dict_i = dict()

            row_i = df_gp_params.iloc[0]

            df_train_dft = df_dft.drop(
                labels=fold_i,
                axis=0)

            df_train_feat = df_feat_post.loc[df_train_dft.index]
            df_test_feat = df_feat_post.loc[fold_i]






            # #################################################
            FP = FingerPrints(
                # df_feat_pre,
                df_feat_post,
                df_features_post=df_feat_post,
                pca_mode="num_comp",  # 'num_comp' or 'perc'
                # pca_comp=10,
                pca_comp=pca_comp_i,
                pca_perc=None,
                # verbose=verbose,
                verbose=False,
                )

            CS = CandidateSpace(
                Y_data=df_dft,
                Y_key="dH",
                FingerPrints=FP,
                )

            # FP.clean_data(df_feat_post, df_feat_post)
            FP.clean_data(df_train_feat, df_test_feat)
            FP.pca_analysis()

            df_train = FP.df_train
            df_test = FP.df_test

            RM = RegressionModel(
                df_train=df_train,
                train_targets=df_train_dft.dH,
                # train_targets=CS.Y_data_series,
                # train_targets=CS.Y_data_series.drop(labels=leave_out_ids),
                df_test=df_test,
                opt_hyperparameters=True,
                gp_settings_dict=gp_settings,
                uncertainty_type='regular',
                verbose=verbose,
                )

            RM.run_regression()

            model = pd.concat([
                CS.Y_data_series,
                RM.model,
                ], axis=1, sort=False)

            model_i = model[~model["y"].isna()]

            mae = np.abs(model_i.dH - model_i.y).mean()
            print(mae)

            data_dict_i["mae"] = mae

            data_dict_list_j.append(data_dict_i)

        # #####################################################
        df_i = pd.DataFrame(data_dict_list_j)
        mae_ave = df_i.mae.mean()

        data_dict_j["mae_ave"] = mae_ave

        data_dict_list.append(data_dict_j)

    df_m = pd.DataFrame(data_dict_list)

# Pickling data ###########################################
import os; import pickle
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)
with open(os.path.join(directory, stoich_i + "_data.pickle"), "wb") as fle:
    pickle.dump(df_m, fle)
# #########################################################

print(20 * "# # ")
print("All done!")
assert False

# + active=""
#
#
#
#
#
#
# -

# # Plotting (Run this cell to plot)

# +
import plotly.graph_objs as go

trace = go.Scatter(
    x=df_m.pca_comp,
    y=df_m.mae_ave,
    mode="markers",
    )
data = [trace]

fig = go.Figure(data=data)
fig.show()

# + active=""
#
#
#
#
#
#
