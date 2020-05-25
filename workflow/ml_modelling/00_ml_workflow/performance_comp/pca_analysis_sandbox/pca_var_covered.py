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

# + [markdown] Collapsed="false"
# # Import Modules

# + Collapsed="false"
import os
print(os.getcwd())
import sys

import pickle

import numpy as np
import pandas as pd

import chart_studio.plotly as py
import plotly.graph_objs as go

from ccf_similarity.ccf import CCF

from active_learning.al_analysis import ALPerformance

from plotting.my_plotly import my_plotly_plot

# from layout import layout

# +
# from inputs import stoich_i

stoich_i = "AB3"

# +
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow/191102_new_workflow"))
from al_data import al_data_files_dict

files_list_gp_ucb = al_data_files_dict[stoich_i]["files_list_gp_ucb"]
files_list_random = al_data_files_dict[stoich_i]["files_list_random"]

# + [markdown] Collapsed="false"
# # Script Inputs

# + Collapsed="false" jupyter={"outputs_hidden": false}
perc_of_structs = 2.5

num_disc = 7

subdirs_list = ["gp_ucb", "random"]

shared_scatter_props = dict(
    mode="lines",

    )

data_path_root = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow",
    "191102_new_workflow/00_abx_al_runs/out_data",
    stoich_i,
    # "gp_ucb_False",
    )


# from inputs import top_ids_to_track_ab2, top_ids_to_track_ab3

# plot_guidlines = False

# if stoich_i == "AB2":
#     top_ids_to_track = top_ids_to_track_ab2
# elif stoich_i == "AB3":
#     top_ids_to_track = top_ids_to_track_ab3
# else:
#     print("ISDJIFSDJI")
# -

subdir=None
shared_scatter_props=None
unique_scatter_props=None
ALPerf_account_duplicates=True
top_ids_to_track=None
files_list=files_list_gp_ucb
color2=None

# +
# assert False

# + Collapsed="false"
# def process_data(
# subdir=None,
# shared_scatter_props=None,
# unique_scatter_props=None,
# ALPerf_account_duplicates=True,
# top_ids_to_track=None,
# files_list=None,
# color2=None,
# ):

# + Collapsed="false"
out_data_dict = dict()

# #############################################################################

data_dict = dict()
for file_i in files_list:
    # #########################################################################
    num = file_i.split("_")[-1].split(".")[0]

    # file_path_i =os.path.join(
    #     data_path_root, subdir, file_i)
    # COMBAK
    # with open(file_path_i, "rb") as fle:
    with open(file_i, "rb") as fle:
        AL_i = pickle.load(fle)

    data_dict[num] = AL_i
out_data_dict["AL_dict"] = data_dict

# +
pca = AL_i.CandidateSpace.FingerPrints.PCA

dir(pca)

# pca.explained_variance_ratio_
pca.explained_variance_ratio_.cumsum()

# + Collapsed="false"
# # #############################################################################
# df_list = []
# for num, AL in data_dict.items():
#     ALPerf = ALPerformance(
#         ALBulkOpt=AL,
#         verbose=False)
#     ALPerf.num_sys_discovered(
#         # perc_of_structs=perc_of_structs,
#         # account_duplicates=ALPerf_account_duplicates,

#         mode="user_specified",  # 'perc' or 'num'
#         # mode="perc",  # 'perc' or 'num'
#         perc_of_structs=perc_of_structs,
#         num_structs=None,
#         ids_to_track=top_ids_to_track,
#         account_duplicates=ALPerf_account_duplicates,

#         )

#     # #########################################################################
#     df = ALPerf.num_sys_discovered_df
#     df_list.append(df)


# df_m = pd.concat(
#     df_list,
#     axis=1,
#     keys=data_dict.keys(),
#     )

# # Checking that the x-axis series are all the same
# # Necessary if the different runs are to be averaged
# x_axis_series_list = []
# for i in data_dict.keys():
#     x_axis_series = df_m[i]["num_dft"].tolist()
#     x_axis_series_list.append(x_axis_series)
# all_x_axis_the_same = all(x_axis_series_list)
# assert all_x_axis_the_same is True, "ISFIDSIFJISDIfj"






# # #############################################################################
# # df_m.index = df_m[0, "num_dft"].tolist()
# df_m.index = x_axis_series


# # #############################################################################
# for i in data_dict.keys():
#     del df_m[i, "num_dft"]


# # TEMP
# # out_data_dict["df_perf"] = df_m
# # return(out_data_dict)


# columns_list = list(df_m.columns.levels[0])
# col = df_m.loc[: , columns_list[0]:columns_list[-1]]

# # col = df_m.loc[: , 0:list(df_m.columns.levels[0])[-1]]

# y_mean = col.mean(axis=1)
# y_std = col.std(axis=1)

# df_ave = pd.DataFrame()
# df_ave["y_mean"] = y_mean
# df_ave["y_std"] = y_std
# df_ave.index = df_m.index

# # Adding 0 to trace
# df_ave.loc[0] = [0, 0]
# df_ave = df_ave.sort_index()

# dx = df_ave.index.values[-1] - df_ave.index.values[-2]
# # last_data_point_ind = df_ave.index.values[-1] + dx
# last_data_point_ind = df_ave.index.values[-1] + 50
# df_ave.loc[last_data_point_ind] = [10, 0]

# # df_ave.loc[260] = [10, 0]
# df_ave = df_ave.sort_index()


# out_data_dict["df_perf"] = df_m

# traces = []
# # #############################################################################
# trace = go.Scatter(
#     x=df_ave.index.tolist(),
#     y=df_ave["y_mean"],
#     line=dict(
#         width=1.,
#         ),
#     )
# trace.update(**shared_scatter_props)
# trace.update(**unique_scatter_props)
# traces.append(trace)
# # #########################################################################
# trace = go.Scatter(
#     x=df_ave.index.tolist(),
#     y=df_ave["y_mean"] + df_ave["y_std"],
#     line=dict(
#         width=0.5,
#         ),
#     )
# trace.update(**shared_scatter_props)
# trace.update(**unique_scatter_props)
# traces.append(trace)

# trace = go.Scatter(
#     x=df_ave.index.tolist(),
#     y=df_ave["y_mean"] - df_ave["y_std"],
#     fill="tonexty",
#     line=dict(
#         width=0.5,
#         # color="red",
#         ),
#     )
# trace.update(**shared_scatter_props)
# trace.update(**unique_scatter_props)
# traces.append(trace)
# # #########################################################################

# trace = go.Scatter(
#     x=df_ave.index.tolist(),
#     y=df_ave["y_mean"] - df_ave["y_std"],
#     line=dict(
#         width=0.5,
#         color=color2,
#         ),
#     )
# trace.update(**shared_scatter_props)
# trace.update(**unique_scatter_props)
# traces.append(trace)

# trace = go.Scatter(
#     x=df_ave.index.tolist(),
#     y=df_ave["y_mean"] + df_ave["y_std"],
#     line=dict(
#         width=0.5,
#         color=color2,
#         ),
#     )
# trace.update(**shared_scatter_props)
# trace.update(**unique_scatter_props)
# traces.append(trace)

# # out_data_dict["trace"] = trace
# out_data_dict["trace"] = traces
# out_data_dict["df_ave"] = df_ave
