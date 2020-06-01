# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
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

# # Import Modules

# + jupyter={}
import os
import sys

sys.path.insert(0,
    os.path.join(
        os.environ["PROJ_irox"],
        "workflow/ml_modelling/00_ml_workflow",
        "190611_new_workflow/02_gaus_proc"))

from gp_methods import gp_workflow

from IPython.display import display

# +
# %%capture
# | - OUT_OF_SIGHT
import os
import sys

import pickle
# import time

# import itertools

import pandas as pd
import numpy as np

# import chart_studio.plotly as py
# import plotly.graph_objs as go

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import (
    bulk_dft_data_path, unique_ids_path,
    df_features_pre_opt_path,
    df_features_post_opt_path)

from gp_methods import gp_model_gpflow, gp_model_catlearn

# from methods import get_trace_j
# from plotting.my_plotly import my_plotly_plot

import pprint
pp = pprint.PrettyPrinter()

from gp_methods import gp_workflow, job_aquisition, test_al_conv

sys.path.insert(0,
    os.path.join(os.environ["PROJ_irox"], "workflow/ml_modelling"))

from ml_methods import create_mixed_df
# -

# # Script Inputs

# +
# stoich_i = "AB2"
stoich_i = "AB3"

# gp_model = gp_model_gpflow
gp_model = gp_model_catlearn

aqs_bin_size = 10

# output_key = "form_e_chris"
output_key = "energy_pa"

verbosity_level = 6  # 1-10 scale

# +
with open(bulk_dft_data_path, "rb") as fle:
    df_bulk_dft = pickle.load(fle)

with open(df_features_pre_opt_path, "rb") as fle:
    df_features_pre = pickle.load(fle)

with open(df_features_post_opt_path, "rb") as fle:
    df_features_post = pickle.load(fle)

df_ids = pd.read_csv(unique_ids_path)
# -

df_bulk_dft.loc[
    ['mony9ibt9r', '8avd8d7rbe', 'zwnung6s71', '7f7svsnpvg', 'x48d9oc591']
    ]

# +
ids_to_drop = ['mony9ibt9r', '8avd8d7rbe', 'zwnung6s71', '7f7svsnpvg', 'x48d9oc591']

import json
directory = "out_data"
if not os.path.exists(directory):
    os.makedirs(directory)
with open(os.path.join(directory, 'outlier_features.json'), 'w') as outfile: 
    json.dump(ids_to_drop, outfile, indent=2)

# +
# #############################################################################
# Filter ids ##################################################################
df_ids = df_ids[
    (df_ids["stoich"] == stoich_i) & \
    (df_ids["source"] != "oqmd") & \
    # (df_ids["source"] != "raul")
    [True for i in range(len(df_ids))]
    ]

# IDS TO DROP
df_ids = df_ids[~df_ids["unique_ids"].isin(ids_to_drop)]
unique_ids = df_ids["unique_ids"].tolist()

# TEMP | Not needed anymore, taken care of on line 9
# unique_ids = [x for x in unique_ids if x not in ids_to_drop]

# #############################################################################
# Training Features ###########################################################
index_filter = np.intersect1d(df_features_post.index, unique_ids)
df_features_post = df_features_post.loc[index_filter]

# #############################################################################
# Training Features ###########################################################
index_filter = np.intersect1d(df_bulk_dft.index, unique_ids)
df_bulk_dft = df_bulk_dft.loc[index_filter]

# #############################################################################
# Test Features ###############################################################
index_filter = np.intersect1d(df_features_pre.index, unique_ids)
df_features_pre = df_features_pre.loc[index_filter]

# #############################################################################
# Filter training data ########################################################
df_features_post = \
    df_features_post[df_features_post["data"]["source"] != "chris"]
df_bulk_dft = df_bulk_dft[df_bulk_dft["source"] != "chris"]

# +
# %%capture
sys.path.insert(0,
    os.path.join(os.environ["PROJ_irox"], "workflow/ml_modelling"))

all_ids = df_features_pre.index.unique()

computed_ids = df_bulk_dft.index.unique()
computed_ids = np.random.choice(computed_ids, size=10)
computed_ids = list(computed_ids)

# TEMP | Use all training data initially
computed_ids = df_bulk_dft.index.tolist()

df_post = df_features_post["voronoi"]
df_pre = df_features_pre["voronoi"]

# +
# t0 = time.time()
# num_training = str(len(computed_ids)).zfill(3)
# step_num = str(i_cnt).zfill(3); i_cnt_str = str(i_cnt).zfill(3)
# print(step_num, " | ", num_training + " " + 68 * "#"); print(80 * "#")
# row_i = df_gp_params.iloc[0]


# #########################################################################
df_test_tmp = create_mixed_df(
    all_ids, computed_ids,
    df_post, df_pre, verbose=False)
# df_test_tmp = df_pre


# #########################################################################
# #########################################################################
# computed_ids = [i for i in computed_ids if i in df_bulk_dft.index]

# computed_ids = df_bulk_dft.index.tolist()
computed_ids = list(set(computed_ids))


# TEMP | Needed because of nan in df_train
# computed_ids = [i for i in computed_ids if i in df_post.index]


df_bulk_dft_i = df_bulk_dft.loc[computed_ids]
# df_train = df_post.loc[computed_ids]

df_train = df_post.loc[computed_ids]

# #########################################################################
# Running GP Model ########################################################
# gp_params_i = row_i.to_dict()

out = gp_workflow(
    df_features_post=df_train, df_test=df_test_tmp,
    df_bulk_dft=df_bulk_dft_i, df_bulk_dft_all=df_bulk_dft,
    df_ids=df_ids, gp_model=gp_model_catlearn,
    opt_hyperparameters=True, gp_params=None,
    y_train_key="energy_pa", run_gp=False)
model_i = out["model"]; model_inst = out["model_inst"]

# + active=""
# train_data.shape: (155, 271)
# train_data.shape: (155, 254)

# +
out.keys()

train_x = out["train_x"]
train_y = out["train_y"]
train_y_standard = out["train_y_standard"]

test_x = out["test_x"]

# +
df_m = pd.concat(
    [
        train_x,
        train_y,
        # train_y_standard,
        ],
    axis=1
    )

test_x = test_x.drop(
    # labels=df_m["id_unique"].unique().tolist(),
    labels=df_m.index.unique().tolist(),
    )

df_m = pd.concat([
    df_m,
    test_x,
    ])

df_m["id"] = df_m.index
df_m = df_m.sort_values("energy_pa")
# -

df_m = df_m.fillna(
    value=-4.51082304,
    method=None,
    axis=None,
    inplace=False,
    limit=None,
    downcast=None,
    # **kwargs,
    )

df_m = df_m.reset_index()

display(df_m.describe())
print("")
display(df_m.head())

# +
import plotly.express as px
df = df_m
fig = px.scatter_3d(
    df,
    x="PCA0",
    y="PCA1",
    z="PCA2",
#     z="energy_pa",
    color='energy_pa',
#     text='id',

    range_x=None,
    range_y=None,
    range_z=None,
    )


fig.show()

# +
fig.data[0]

df_m.sort_values("PCA1", ascending=False)[0:5]["index"].tolist()
