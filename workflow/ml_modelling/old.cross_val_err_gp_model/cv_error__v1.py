# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.2.1
#   kernelspec:
#     display_name: Python [conda env:PROJ_IrOx_Active_Learning_OER]
#     language: python
#     name: conda-env-PROJ_IrOx_Active_Learning_OER-py
# ---

# # New ML Active Learning Workflow
# ---
#
# A model that predicts the mean (~ -6.05 eV/atom) has a MAE of ~0.3 eV/atom)

# # Import Modules

# + {"jupyter": {"source_hidden": true}}
# %%capture
#| - OUT_OF_SIGHT
import os
import sys

import pickle
import time

import itertools

import pandas as pd
import numpy as np

import chart_studio.plotly as py
import plotly.graph_objs as go

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import (
    bulk_dft_data_path, unique_ids_path,
    df_features_pre_opt_path,
    df_features_post_opt_path)


from plotting.my_plotly import my_plotly_plot

import pprint
pp = pprint.PrettyPrinter()


sys.path.insert(0,
    os.path.join(os.environ["PROJ_irox"], "workflow/ml_modelling"))

sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow/190611_new_workflow/02_gaus_proc"))

from gp_methods import gp_workflow, job_aquisition, test_al_conv
from methods import get_trace_j
from gp_methods import gp_model_gpflow, gp_model_catlearn
from ml_methods import create_mixed_df
# -

# # Script Inputs

# + {"jupyter": {"source_hidden": true}}
# stoich_i = "AB2"
stoich_i = "AB3"

# gp_model = gp_model_gpflow
gp_model = gp_model_catlearn

aqs_bin_size = 5

# output_key = "form_e_chris"
output_key = "energy_pa"

verbosity_level = 6  # 1-10 scale

# +
params_dict = {
#     "noise": [0.02542],
#     "sigma_l": [0.0049],
#     "sigma_f": [5.19],
#     "alpha": [0.018],

    "noise": [0.0001],
    "sigma_l": [10.],
    "sigma_f": [5],
    "alpha": [0.1],

    }


c = list(itertools.product(*params_dict.values()))
df_gp_params = pd.DataFrame(c, columns=params_dict.keys())
# -

# # Read Data

# + {"jupyter": {"source_hidden": true}}
with open(bulk_dft_data_path, "rb") as fle:
    df_bulk_dft = pickle.load(fle)

with open(df_features_pre_opt_path, "rb") as fle:
    df_features_pre = pickle.load(fle)

with open(df_features_post_opt_path, "rb") as fle:
    df_features_post = pickle.load(fle)

df_ids = pd.read_csv(unique_ids_path)

print("df_ids.shape:", df_ids.shape)
# + {"jupyter": {"source_hidden": true}}
# iro2_indices_already_computed = df_bulk_dft[df_bulk_dft["stoich"] == "AB2"].index


# df_ids_tmp = df_ids[df_ids["stoich"] == "AB2"].set_index("unique_ids")
# df_ids_tmp.drop(iro2_indices_already_computed)[0:50]["id"].tolist()
# -

# # Systems to Discard

# + {"jupyter": {"source_hidden": true}}
import pickle
import os

path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/ccf_similarity_analysis/out_data",
    "all_ids_to_elim.pickle")
with open(path_i, "rb") as fle:
    ids_to_drop = pickle.load(fle)

# + {"jupyter": {"source_hidden": true}}
import pickle
import os
import json


path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/ccf_similarity_analysis/out_data",
    "all_ids_to_elim.pickle")
with open(path_i, "rb") as fle:
    ids_to_drop__duplicates = pickle.load(fle)
    
# #############################################################################
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/visualizing_data/out_data",
    "outlier_features.json")
with open(path_i, 'r') as f:
    ids_to_drop__outliers = json.load(f)

# #############################################################################
ids_to_drop = ids_to_drop__outliers + ids_to_drop__duplicates
print("len(ids_to_drop):", len(ids_to_drop))
ids_to_drop = list(set(ids_to_drop))
print("len(ids_to_drop):", len(ids_to_drop))
# -

# # Filtering dataframes to the correct stoicheometry

# +
# #############################################################################
# Filter ids ##################################################################
print("df_ids.shape:", df_ids.shape)
df_ids = df_ids[
    (df_ids["stoich"] == stoich_i) & \
    (df_ids["source"] != "oqmd") & \
    # (df_ids["source"] != "raul")
    [True for i in range(len(df_ids))]
    ]

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
df_bulk_dft = df_bulk_dft.loc[index_filter]
print("df_bulk_dft.shape:", df_bulk_dft.shape)

# #############################################################################
# Test Features ###############################################################
index_filter = np.intersect1d(df_features_pre.index, unique_ids)
df_features_pre = df_features_pre.loc[index_filter]

# #############################################################################
# Filter training data ########################################################
# df_features_post = \
#     df_features_post[df_features_post["data"]["source"] != "chris"]
# df_bulk_dft = df_bulk_dft[df_bulk_dft["source"] != "chris"]

# +
# %%capture
sys.path.insert(0,
    os.path.join(os.environ["PROJ_irox"], "workflow/ml_modelling"))

all_ids = df_features_pre.index.unique()

computed_ids = df_bulk_dft.index.unique()
computed_ids = np.random.choice(computed_ids, size=10)
computed_ids = list(computed_ids)

df_post = df_features_post["voronoi"]
df_pre = df_features_pre["voronoi"]
# -

# # Preparing N-fold CV Folds

# +
print("df_bulk_dft:", df_bulk_dft.shape)

n_fold_cv = df_bulk_dft.shape[0]
n_fold_cv = 5


fold_size = int(df_bulk_dft.shape[0] / n_fold_cv)
print("fold_size:", fold_size)

# Shuffling training data
df_bulk_dft = df_bulk_dft.sample(
    n=None,
    frac=1.,
    replace=False,
    axis=None)

print("n_fold_cv * fold_size:", n_fold_cv * fold_size)

ids_0 = df_bulk_dft.index[:n_fold_cv * fold_size]
folds = np.split(ids_0, n_fold_cv)

ids_leftover = df_bulk_dft.index[n_fold_cv * fold_size:]

if ids_leftover.shape[0] > 0:
    folds.append(ids_leftover)

folds = np.array(folds)

# + {"active": ""}
#
#
#

# +
rows_list = []
models_inst_list = []
out_list = []
# for i_cnt, (name_i, row_i) in enumerate(df_bulk_dft.iterrows()):
for i_cnt, fold_i in enumerate(folds):
    #| - GP AL Iteration ******************************************************
    # *************************************************************************
    # *************************************************************************
    print("")
    t0 = time.time()
    num_training = str(len(fold_i)).zfill(3)
    step_num = str(i_cnt).zfill(3);
    print(step_num, " | ", num_training + " " + 68 * "#"); print(80 * "#")
    row_i = df_gp_params.iloc[0]


    df_bulk_dft_i = df_bulk_dft.drop(
        # labels=name_i,
        labels=fold_i,
        axis=0)

    df_train = df_post.loc[df_bulk_dft_i.index]
    # df_test_tmp = df_post.loc[[name_i]]
    df_test_tmp = df_post.loc[fold_i]


    # #########################################################################
    # Running GP Model ########################################################
    if True:
#     try:
        gp_params_i = row_i.to_dict()
        out = gp_workflow(
            df_features_post=df_train,
            df_test=df_test_tmp,
            df_bulk_dft=df_bulk_dft_i,
            df_bulk_dft_all=df_bulk_dft,

            df_ids=df_ids,
            gp_model=gp_model_catlearn,
            opt_hyperparameters=True,
            gp_params=gp_params_i,
            y_train_key="energy_pa",

            verbose=False,

            clean_variance_flag=True,
            clean_skewness_flag=True,
            clean_infinite_flag=True,
            standardize_data_flag=True,


            pca_comp=11,
            # pca_comp=11,
            pca_perc=0.99,
            pca_mode="num_comp",
            # pca_mode="perc",
            ); out_list.append(out)

        model_i = out["model"]; model_inst = out["model_inst"]


        models_inst_list.append(model_inst)
        test_row_i = model_i[model_i["prediction"].notnull()]
        rows_list.append(test_row_i)


        mae_i = abs(test_row_i["prediction_unstandardized"] - test_row_i["energy_pa"]).mean()
        print(
            "MAE_i: ",
            mae_i)

        print("")
        print("model_inst.regularization", model_inst.regularization)
        
        print("model_inst.kernel_list:", model_inst.kernel_list)

#         try:
#             print("kernel_list['slope']:", model_inst.kernel_list[0]["slope"])
#             print("kernel_list['scaling']", model_inst.kernel_list[0]["scaling"])
#             print("kernel_list['degree']", model_inst.kernel_list[0]["degree"][0])
#         except:
#             pass

#         try:
#             print("width:", model_inst.kernel_list[0]["width"])
#             print("scaling", model_inst.kernel_list[0]["scaling"])
#         except:
#             pass

        print(
            "model_inst.log_marginal_likelihood: ",
            model_inst.log_marginal_likelihood)

#     except:
#         print("Failed!!")
# -

# # Analyzing PCA

# +
pca = out["pca"]

columns_pre_pca = out["df_train_pre_pca"].columns

df_pca_comp = pd.DataFrame(
    pca.components_,
    columns=columns_pre_pca,
    # index = ['PC-1','PC-2']
    index=["PCA_" + str(i).zfill(2) for i in range(pca.n_components)],
    )

df_pca_comp = abs(df_pca_comp.T)

# df_pca_comp
# df_pca_comp.sort_values("PCA_00", ascending=False)
# -

# # Calc MAE

# +
df_all_pred = pd.concat(rows_list)

df_all_pred["error"] = abs(
    df_all_pred["prediction_unstandardized"] - df_all_pred["energy_pa"])

mae = df_all_pred["error"].mean()

print("MAE (eV): ", mae)

# + {"active": ""}
# MAE (eV):  0.14186230367245345
#
# MAE (eV):  0.1352335205391316
#
# MAE (eV):  0.12055492648019178
# -

df_all_pred.head()

# + {"active": ""}
#
#
#

# + {"active": ""}
# # Without any optimization
# MAE (eV):  0.08421070878633005
#
# # With only local opt
# MAE (eV):  0.08245272125264294
#
# # With Global opt (a lot failed)
# MAE (eV):  0.07163131081485057

# + {"active": ""}
# # Removing duplicates via CCF analysis (3 PCA components)
# MAE (eV):  0.1946592371070702
#
#
# # Removing duplicates via CCF analysis (99% variance captured w/ PCA)
# MAE (eV):  0.1341779368327114
# -

# # 5 PCA Components | Optmization

# +
# Optimizatino on
0.13447096075807743

# Global Optimization on
0.13244164782822093

# No optimization
0.13184483308309516

# +
mae_pca_dict = {
    30: 0.1831923511185789,
    20: 0.14531282803281226,
    15: 0.13394643100712159,

    13: 0.12004683884649776,
    12: 0.11351107652779499,
    11: 0.11091405471923632,

    # 10: 0.13401349790117334,    
    10: 0.1298666071766859,

    9: 0.13246529211025734,
    8: 0.13562396520300607,
    7: 0.13272655356364796,
    6: 0.1333824520990277,
    5: 0.13447096075807743,
    4: 0.1384055763748185,
    3: 0.15767416200944745,
    2: 0.22931482558271382,
    1: 0.21253287970745607,
    }

df = pd.DataFrame(mae_pca_dict, index=["mae"]).T

# +
import chart_studio.plotly as py
import plotly.graph_objs as go
import os

y_array = df["mae"]
x_array = df.index.tolist()

trace = go.Scatter(
    x=x_array,
    y=y_array,
    mode="markers",

    marker=dict(
        symbol="circle",
        color='LightSkyBlue',
#         colorscale='Viridis',
        colorbar=dict(thickness=20),
        size=20,
        line=dict(
            color='MediumPurple',
            width=2
            )
        ),
    )

data = [trace]

fig = go.Figure(data=data)
fig.show()
