# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_IrOx_Active_Learning_OER]
#     language: python
#     name: conda-env-PROJ_IrOx_Active_Learning_OER-py
# ---

# +
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

stoich_i = "AB3"
verbose = True
name_i = "TEMP"
save_dir_extra=None
acquisition_method=None
duplicate_analysis=None
seed=None

# + endofcell="--"
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

# # + {"jupyter": {"source_hidden": true}}
ids_w_dft = df_bulk_dft.index

# TEMP | Reduce size of candidate space
# np.random.seed(8)
# ids_w_dft = np.sort(np.random.choice(np.sort(ids_w_dft), size=200))
ids_w_dft = list(set(ids_w_dft))
# print("ids_w_dft:", ids_w_dft)

df_bulk_dft = df_bulk_dft.loc[ids_w_dft]

df_features_pre = df_features_pre.loc[ids_w_dft]
df_features_post = df_features_post.loc[ids_w_dft]
# --

# +
# def chunks(l, n):
#     """Yield successive n-sized chunks from l."""
#     for i in range(0, len(l), n):
#         yield(l[i:i + n])

# models_list = []
# for i_cnt, i in enumerate(chunks(ids_w_dft, 100)):
    
# leave_out_ids = i


# df_post_i = df_features_post.drop(labels=leave_out_ids)
# df_features_post
# leave_out_ids

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

df_train = CS.FingerPrints.df_post
# .drop(labels=leave_out_ids)
#     df_test = CS.FingerPrints.df_post.loc[leave_out_ids]
df_test = CS.FingerPrints.df_pre

FP.clean_data(df_train, df_test)
FP.pca_analysis()

df_train = FP.df_train
df_test = FP.df_test



gp_settings = {
    "noise": 0.02542,
    "sigma_l": 1.0049,
    "sigma_f": 5.19,
    "alpha": 0.018,
    }

RM = RegressionModel(
    df_train=df_train,
    train_targets=CS.Y_data_series,
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


# models_list.append(model_i)
# -

# # Plotting

# +
# models_list[0]

# model_master = pd.concat(models_list, axis=0, sort=False)

# +
import chart_studio.plotly as py
import plotly.graph_objs as go

import os

x_array = model.y
y_array = model.y_real

trace = go.Scatter(
    x=x_array,
    y=y_array,
    mode="markers",

    marker=dict(
        symbol="circle",
        color='blue',

        # color=z,
        # colorscale='Viridis',
        # colorbar=dict(thickness=20),

        size=4,
        line=dict(
            color='black',
            width=1,
            )
        ),

    # line=dict(
    #     color="firebrick",
    #     width=2,
    #     dash="dot",
    #     ),

    # error_y={
    #     "type": 'data',
    #     "array": [0.4, 0.9, 0.3, 1.1],
    #     "visible": True,
    #     },

    )
trace_xy = go.Scatter(x=[-3, 5], y=[-3, 5], mode="lines")
data = [trace, trace_xy]

layout = go.Layout(
    xaxis=dict(range=[-3.5, 6]),
    yaxis=dict(range=[-3.5, 6])
    )

fig = go.Figure(data=data, layout=layout)
fig.show()

# +
model_master = model
model_master["err_pred_real"] = np.abs(model_master["y_real"] - model_master["y"])
model_master["err_pred_real"].mean()

# model_master
# -1.658857 - -1.784430	
# 0.12557299999999993

# + active=""
# 0.20012585885222797
# -

# Pickling data ######################################################
import os; import pickle
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)
with open(os.path.join(directory, stoich_i + "_" + "pre_dft_cv_data.pickle"), "wb") as fle:
    pickle.dump(model, fle)
# #####################################################################
