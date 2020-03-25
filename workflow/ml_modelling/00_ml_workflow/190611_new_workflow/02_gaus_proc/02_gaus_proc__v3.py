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

# # New ML Active Learning Workflow
# ---
#
# 32 (Too big, not computed),
#
# 226 (Not computed),

# # Import Modules

# +
# %%capture

# | - OUT_OF_SIGHT
import os
import sys

import json
import pickle

import time
t0_init = time.time()

import itertools

import pandas as pd
import numpy as np

from ase.visualize import view

import chart_studio.plotly as py
import plotly.graph_objs as go

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import (
    bulk_dft_data_path, unique_ids_path,
    df_features_pre_opt_path,
    df_features_pre_opt_kirsten_path,
    df_features_post_opt_path,
    ids_to_discard__too_many_atoms_path,
    )

from gp_methods import gp_model_gpflow, gp_model_catlearn

from methods import get_trace_j
from plotting.my_plotly import my_plotly_plot

import pprint
pp = pprint.PrettyPrinter()

from gp_methods import gp_workflow, job_aquisition, test_al_conv

sys.path.insert(0,
    os.path.join(os.environ["PROJ_irox"], "workflow/ml_modelling"))

from ml_methods import create_mixed_df

from ase_modules.ase_methods import view_in_vesta

sys.path.insert(0, "../04_final_ml_plots")
from layout import get_layout

layout = get_layout(model=None)

from methods import plot_model
# -

from gp_methods import random_job_aquisition

# # Script Inputs

# +
stoich_i = "AB3"
# stoich_i = "AB2"

custom_name = "regular"
# custom_name = "random"

# gp_model = gp_model_gpflow
gp_model = gp_model_catlearn

# aqs_bin_size = 5
aqs_bin_size = 10

# output_key = "form_e_chris"
output_key = "energy_pa"

verbosity_level = 6  # 1-10 scale

# +
params_dict = {

    # "noise": [0.02542],
    # "sigma_l": [0.0049],
    # "sigma_f": [5.19],
    # "alpha": [0.018],

    # "noise": [0.000001],
    # "sigma_l": [0.001],
    # "sigma_f": [0.5],
    # "alpha": [0.5],

    # Good for AB3
    "noise": [0.02542],
    "sigma_l": [1.0049],
    "sigma_f": [5.19],
    "alpha": [0.018],

    # "noise": [0.02542],
    # "sigma_l": [10.0049],
    # "sigma_f": [5.19],
    # "alpha": [0.018],

    }

# noise = 0.0042  # Regularisation parameter.
# sigma_l = 6.3917  # Length scale parameter.
# sigma_f = 0.5120  # Scaling parameter.
# alpha = 0.3907  # Alpha parameter.


c = list(itertools.product(*params_dict.values()))
df_gp_params = pd.DataFrame(c, columns=params_dict.keys())
# -

dir_i = "out_data";
if not os.path.exists(dir_i):
    os.makedirs(dir_i)

# # Read Data

# +
with open(bulk_dft_data_path, "rb") as fle:
    df_bulk_dft = pickle.load(fle)

with open(df_features_pre_opt_path, "rb") as fle:
    df_features_pre = pickle.load(fle)

# with open(df_features_pre_opt_kirsten_path, "rb") as fle:
#     df_features_pre = pickle.load(fle)

with open(df_features_post_opt_path, "rb") as fle:
    df_features_post = pickle.load(fle)

df_ids = pd.read_csv(unique_ids_path)

df_ids = df_ids[
    (df_ids["stoich"] == stoich_i) & \
    (df_ids["source"] != "oqmd") & \
    (df_ids["source"] != "raul_oer") & \
    [True for i in range(len(df_ids))]]
# -
# # Filtering dataframes to the correct stoicheometry

# # TEMP DROP DUPLICATE and OUTLIER SYSTEMS

# ## Collect ids to ignore list

# +
df_ids[df_ids["unique_ids"] == 'zwmivrzazu']

# 'zwmivrzazu' in df_ids

# df_ids["unique_ids"].tolist()

# +
# # df_ids.index.unique().tolist()

# print("len(ids_to_drop__duplicates):", len(ids_to_drop__duplicates))
# ids_to_drop__duplicates = [i for i in ids_to_drop__duplicates if i in df_ids["unique_ids"].unique().tolist()]
# print("len(ids_to_drop__duplicates):", len(ids_to_drop__duplicates))

# +
# #############################################################################
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/ccf_similarity_analysis/out_data",
    "all_ids_to_elim_1.pickle")
with open(path_i, "rb") as fle:
    ids_to_drop__duplicates = pickle.load(fle)
    ids_to_drop__duplicates = ids_to_drop__duplicates[stoich_i]

    print("len(ids_to_drop__duplicates):", len(ids_to_drop__duplicates))
    ids_to_drop__duplicates = [i for i in ids_to_drop__duplicates if i in df_ids["unique_ids"].unique().tolist()]
    print("len(ids_to_drop__duplicates):", len(ids_to_drop__duplicates))

# #############################################################################
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/visualizing_data/out_data",
    "outlier_features.json")
with open(path_i, 'r') as f:
    ids_to_drop__outliers = json.load(f)

    print("len(ids_to_drop__outliers):", len(ids_to_drop__outliers))
    ids_to_drop__outliers = [i for i in ids_to_drop__outliers if i in df_ids["unique_ids"].unique().tolist()]
    print("len(ids_to_drop__outliers):", len(ids_to_drop__outliers))


# #############################################################################
with open(ids_to_discard__too_many_atoms_path, "rb") as fle:
    ids_to_drop__too_many_atoms = pickle.load(fle)

    print("len(ids_to_drop__too_many_atoms):", len(ids_to_drop__too_many_atoms))
    ids_to_drop__too_many_atoms = [i for i in ids_to_drop__too_many_atoms if i in df_ids["unique_ids"].unique().tolist()]
    print("len(ids_to_drop__too_many_atoms):", len(ids_to_drop__too_many_atoms))

# #############################################################################
ids_to_drop = [] + \
    ids_to_drop__too_many_atoms
    # ids_to_drop__duplicates + \
    # ids_to_drop__outliers + \

print("len(ids_to_drop):", len(ids_to_drop))
ids_to_drop = list(set(ids_to_drop))
print("len(ids_to_drop):", len(ids_to_drop))

abx_ids = df_ids[df_ids["stoich"] == stoich_i]["unique_ids"].tolist()
ids_to_drop = [i for i in ids_to_drop if i in abx_ids]
print("len(ids_to_drop):", len(ids_to_drop))

# +
# #############################################################################
# Filter ids ##################################################################
# df_ids = df_ids[
#     (df_ids["stoich"] == stoich_i) & \
#     (df_ids["source"] != "oqmd") & \
#     (df_ids["source"] != "raul_oer") & \
#     [True for i in range(len(df_ids))]
#     ]

print("df_ids.shape:", df_ids.shape)
# IDS TO DROP
df_ids = df_ids[~df_ids["unique_ids"].isin(ids_to_drop)]
print("df_ids.shape:", df_ids.shape)
unique_ids = df_ids["unique_ids"].tolist()

# #############################################################################
# Training Features ###########################################################
index_filter = np.intersect1d(df_features_post.index, unique_ids)
df_features_post = df_features_post.loc[index_filter]

# #############################################################################
# Training Features ###########################################################
index_filter = np.intersect1d(df_bulk_dft.index, unique_ids)
print("df_bulk_dft.shape:", df_bulk_dft.shape)
df_bulk_dft = df_bulk_dft.loc[index_filter]
print("df_bulk_dft.shape:", df_bulk_dft.shape)

# #############################################################################
# Test Features ###############################################################
index_filter = np.intersect1d(df_features_pre.index, unique_ids)
df_features_pre = df_features_pre.loc[index_filter]

# #############################################################################
# Filter training data ########################################################
df_features_post = \
    df_features_post[df_features_post["data"]["source"] != "chris"]
df_bulk_dft = df_bulk_dft[df_bulk_dft["source"] != "chris"]

# #############################################################################
df_post = df_features_post["voronoi"]
df_pre = df_features_pre["voronoi"]
# -

# # Get initial computed ids randomly

# +
all_ids = df_features_pre.index.unique()
print("len(all_ids):", len(all_ids))

computed_ids = df_bulk_dft.index.unique()
computed_ids = np.random.choice(computed_ids, size=20)
computed_ids = list(set(computed_ids))[0:11]
print("len(computed_ids):", len(computed_ids))

# TEMP | Use all training data initially **************************************
computed_ids = df_bulk_dft.index.tolist()
#__|

# +
# assert False
# -

# # Running GP models

# +
# layout["xaxis"]["range"] = None
# layout["yaxis"]["range"] = None
# layout["showlegend"] = True

# +
# df_features_post=df_train, df_test=df_test_tmp,
# df_bulk_dft=df_bulk_dft_i, df_bulk_dft_all=df_bulk_dft,
# df_ids=df_ids, gp_model=gp_model_catlearn,
# opt_hyperparameters=True, gp_params=gp_params_i,
# y_train_key="energy_pa", verbose=False, pca_comp=11,
# pca_mode="num_comp")

# df_train
# df_test_tmp
# df_bulk_dft_i
# df_bulk_dft
# df_ids
# -

df_bulk_dft.shape

# +
data_dict = dict()
for al_iter_i in range(80):
# for al_iter_i in range(5):
    # | - GP AL Iteration ******************************************************
    # *************************************************************************
    t0 = time.time(); num_training = str(len(computed_ids)).zfill(3)
    al_iter_i_str = str(al_iter_i).zfill(3)
    print(al_iter_i_str, " | ", num_training + " " + 68 * "#"); print(80 * "#")
    row_i = df_gp_params.iloc[0]


    # #########################################################################
    # #########################################################################
    df_test_tmp = create_mixed_df(
        all_ids, computed_ids,
        df_post, df_pre, verbose=False)


    # #########################################################################
    # Filter training data by 'computed' ids ##################################
    computed_ids = list(set(computed_ids))
    df_bulk_dft_i = df_bulk_dft.loc[computed_ids]
    df_train = df_post.loc[computed_ids]


    # #########################################################################
    # Running GP Model ########################################################
    gp_params_i = row_i.to_dict()
    out = gp_workflow(
        df_features_post=df_train, df_test=df_test_tmp,
        df_bulk_dft=df_bulk_dft_i, df_bulk_dft_all=df_bulk_dft,
        df_ids=df_ids, gp_model=gp_model_catlearn,
        opt_hyperparameters=True, gp_params=gp_params_i,
        y_train_key="energy_pa", verbose=True, pca_comp=11,
        pca_mode="num_comp")

    model_i = out["model"]; model_inst = out["model_inst"]


    # #########################################################################
    # Job Aquisition ##########################################################    
    if custom_name == "random":
        aquisition_out_dict = random_job_aquisition(
            model_i, aqs_bin_size=aqs_bin_size,
            df_bulk_dft_all=df_bulk_dft,
            y_train_key="energy_pa")
    else:
        aquisition_out_dict = job_aquisition(
            model_i, aqs_bin_size=aqs_bin_size,
            df_bulk_dft_all=df_bulk_dft,
            y_train_key="energy_pa")

    new_ids_to_compute = aquisition_out_dict["new_ids_to_compute"]
    computed_ids += new_ids_to_compute


    # #########################################################################
    # Test for AL Convergence #################################################
    al_converged = test_al_conv(model_i)
    if al_converged:
        print("AL CONVERGED!!!!!!!!!!!!!!!!!!!!!!!")


    # #########################################################################
    # SAVE DATA ###############################################################
    out_dict = {**out,
        "gp_instance": model_inst, "computed_ids": computed_ids,
        "aquisition_data": aquisition_out_dict, "al_converged": al_converged}
    data_dict[al_iter_i] = out_dict


    file_name_i = "data_dict_" + stoich_i + "_" + custom_name + ".pickle"
    with open(os.path.join(dir_i, file_name_i), "wb") as fle:
        pickle.dump(data_dict, fle)


    # #########################################################################
    # BREAKING LOOP WHEN DFT DATA RUNS OUT ####################################
    if len(new_ids_to_compute) == 0:
        print("NO MORE DFT DATA AVAILABLE")
        break


    # #########################################################################
    # Printing loop time info #################################################
    if verbosity_level > 5:
        print("Loop time (sec):", (time.time() - t0), "\n");

    #__| **********************************************************************

# +
# #########################################################################
# CREATE FIGURE ###########################################################
# COMBAK | Don't need to create the figure in this step | TODO
# fig_i = plot_model(
#     model_i, layout=get_layout(model_i), model_i=model_i,
#     df_bulk_dft=df_bulk_dft, name=al_iter_i_str, custom_text=num_training)
# -

assert False

# +
# df_train.shape

df_bulk_dft_i.shape
df_train.shape

# computed_ids

# + {"active": ""}
#
#
#
#
# -

# # ID's that were needed but are not computed

# +
ids_needed_all = []
for key, val in data_dict.items():
    ids_needed = val["aquisition_data"]["ids_needed_but_not_avail"]
#     ids_needed = val["ids_needed_but_not_avail"]
    if not val["al_converged"]:
        ids_needed_all += ids_needed

print("IDs that were needed during AL but were not computed",
    "\n", list(set(ids_needed_all)))
# -

# # Plotting last AL generation

# +
# gen_i = data_dict[
#     list(data_dict.keys())[-1]]

# gen_i["fig"]
# -

# # Plotting the Log Marginal Likelihood for each model

# +
lml_list = []
for key, val in data_dict.items():
    model = val["model"]
    gp_instance = val["gp_instance"]

    # #########################################################################
    gp_instance.kernel_list

    # if type(gp_instance.log_marginal_likelihood) == list:
    try:
        lml_i = gp_instance.log_marginal_likelihood[0]
    except:
    # else:
        lml_i = gp_instance.log_marginal_likelihood
    lml_list.append(lml_i)


trace = go.Scatter(
    y=lml_list,
    mode="markers")
data = [trace]

fig = go.Figure(data=data)
# fig.show()
# -
print("Notebook executed in ", time.time() - t0_init, "(s)")

# + {"active": ""}
#
#
#
#

# +
# #########################################################################
# #########################################################################
computed_ids = []
df_test_tmp = create_mixed_df(
    all_ids, computed_ids,
    df_post, df_pre, verbose=False)

# #########################################################################
# Filter training data by 'computed' ids ##################################
# computed_ids = list(set(computed_ids))

# df_bulk_dft_i = df_bulk_dft.loc[computed_ids]
df_bulk_dft_i = df_bulk_dft

# df_train = df_post.loc[computed_ids]
df_train = df_post


intersected_ids = list(set(df_bulk_dft_i.index) & set(df_train.index))
df_train = df_train.loc[intersected_ids]
df_bulk_dft_i = df_bulk_dft_i.loc[intersected_ids]

# #########################################################################
# Running GP Model ########################################################
gp_params_i = row_i.to_dict()
out = gp_workflow(
    df_features_post=df_train, df_test=df_test_tmp,
    df_bulk_dft=df_bulk_dft_i, df_bulk_dft_all=df_bulk_dft,
    df_ids=df_ids, gp_model=gp_model_catlearn,
    opt_hyperparameters=True, gp_params=gp_params_i,
    y_train_key="energy_pa", verbose=False, pca_comp=11,
    pca_mode="num_comp")

model_i = out["model"]; model_inst = out["model_inst"]

model_i = model_i[~model_i["energy_pa"].isnull()]

# +
# (abs(model_i["prediction_unstandardized"] - model_i["energy_pa"]) ** 2)
y_pred = model_i["prediction_unstandardized"]
y_test = model_i["energy_pa"]

from sklearn.metrics import mean_squared_error
rms = np.sqrt(mean_squared_error(y_test, y_pred))
print(rms)

# 0.5370373324014677

# +
import chart_studio.plotly as py
import plotly.graph_objs as go

import os

# x_array = [0, 1, 2, 3]
# y_array = [0, 1, 2, 3]

x_array = model_i["prediction_unstandardized"]
y_array = model_i["energy_pa"]

trace = go.Scatter(
    x=x_array,
    y=y_array,
    mode="markers",

    marker=dict(
        symbol="circle",
        color='LightSkyBlue',
        size=14,
        line=dict(
            color='MediumPurple',
            width=2
            )
        ),

    line=dict(
        color="firebrick",
        width=2,
        dash="dot",
        ),

    )


trace_0 = go.Scatter(
    x=[-7, -4],
    y=[-7, -4],
    mode="lines",

#     marker=dict(
#         symbol="circle",
#         color='LightSkyBlue',
#         size=14,
#         line=dict(
#             color='MediumPurple',
#             width=2
#             )
#         ),

    line=dict(
        color="firebrick",
        width=2,
        dash="dot",
        ),

    )
data = [trace, trace_0]

fig = go.Figure(data=data)
fig.show()
