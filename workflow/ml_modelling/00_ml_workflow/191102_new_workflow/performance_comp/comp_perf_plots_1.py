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

# +
import os
import sys

import pickle

import pandas as pd

import chart_studio.plotly as py
import plotly.graph_objs as go

from ccf_similarity.ccf import CCF

from active_learning.al_analysis import ALPerformance

from plotting.my_plotly import my_plotly_plot

from layout import layout
# -

# # Script Inputs

# +
stoich_i = "AB3"

perc_of_structs = 2.5

subdirs_list = ["gp_ucb", "random"]

shared_scatter_props = dict(
    mode="lines",
    marker=dict(
        symbol="circle",
        size=6,
        opacity=0.5,
        line=dict(
            color='black',
            width=1,
            )
        ),

    line=dict(
        width=1.5,
        # dash="dash",
        ),

    error_y=dict(
        thickness=.5,
        width=1.0,
        ),
    )

data_path_root = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow",
    "191102_new_workflow/00_abx_al_runs/out_data",
    stoich_i,
    # "gp_ucb_False",
    )


from inputs import top_ids_to_track_ab2, top_ids_to_track_ab3
# -

if stoich_i == "AB2":
    top_ids_to_track = top_ids_to_track_ab2
elif stoich_i == "AB3":
    top_ids_to_track = top_ids_to_track_ab3
else:
    print("ISDJIFSDJI")


def process_data(
    subdir=None,
    shared_scatter_props=None,
    unique_scatter_props=None,
    ALPerf_account_duplicates=True,
    top_ids_to_track=None,
    ):
    """
    """
    #| - process_data
    out_data_dict = dict()

    # #############################################################################
    files_list = os.listdir(
        os.path.join(
            # dir_i,
            data_path_root,
            # "out_data",
            subdir))
    files_list = [i for i in files_list if "pickle" in i]
    files_list = [i for i in files_list if "AL_" in i]

    data_dict = dict()
    for file_i in files_list:
        # #########################################################################
        num = file_i.split("_")[-1].split(".")[0]

        file_path_i =os.path.join(
            data_path_root, subdir, file_i)
        with open(file_path_i, "rb") as fle:
            AL_i = pickle.load(fle)

        data_dict[num] = AL_i
    out_data_dict["AL_dict"] = data_dict

    # #############################################################################
    df_list = []
    for num, AL in data_dict.items():
        ALPerf = ALPerformance(
            ALBulkOpt=AL,
            verbose=False)
        ALPerf.num_sys_discovered(
            # perc_of_structs=perc_of_structs,
            # account_duplicates=ALPerf_account_duplicates,

            mode="user_specified",  # 'perc' or 'num'
            # mode="perc",  # 'perc' or 'num'
            perc_of_structs=perc_of_structs,
            num_structs=None,
            ids_to_track=top_ids_to_track,
            account_duplicates=ALPerf_account_duplicates,

            )

        # #########################################################################
        df = ALPerf.num_sys_discovered_df
        df_list.append(df)


    df_m = pd.concat(
        df_list,
        axis=1,
        keys=data_dict.keys(),
        )

    # Checking that the x-axis series are all the same
    # Necessary if the different runs are to be averaged
    x_axis_series_list = []
    for i in data_dict.keys():
        x_axis_series = df_m[i]["num_dft"].tolist()
        x_axis_series_list.append(x_axis_series)
    all_x_axis_the_same = all(x_axis_series_list)
    assert all_x_axis_the_same is True, "ISFIDSIFJISDIfj"






    # #############################################################################
    # df_m.index = df_m[0, "num_dft"].tolist()
    df_m.index = x_axis_series


    # #############################################################################
    for i in data_dict.keys():
        del df_m[i, "num_dft"]

        
    # TEMP
    # out_data_dict["df_perf"] = df_m
    # return(out_data_dict)
    
    
    columns_list = list(df_m.columns.levels[0])
    col = df_m.loc[: , columns_list[0]:columns_list[-1]]
 
    # col = df_m.loc[: , 0:list(df_m.columns.levels[0])[-1]]

    y_mean = col.mean(axis=1)
    y_std = col.std(axis=1)

    df_ave = pd.DataFrame()
    df_ave["y_mean"] = y_mean
    df_ave["y_std"] = y_std
    df_ave.index = df_m.index

    out_data_dict["df_perf"] = df_m

    traces = []
    # #############################################################################
    trace = go.Scatter(
        x=df_ave.index.tolist(),
        y=df_ave["y_mean"],
        # error_y={
        #     "array": df_ave["y_std"],
        #     },
        )
    trace.update(**shared_scatter_props)
    trace.update(**unique_scatter_props)
    traces.append(trace)
    # #########################################################################
    trace = go.Scatter(
        x=df_ave.index.tolist(),
        y=df_ave["y_mean"] + df_ave["y_std"],
        )
    trace.update(**shared_scatter_props)
    trace.update(**unique_scatter_props)
    traces.append(trace)

    trace = go.Scatter(
        x=df_ave.index.tolist(),
        y=df_ave["y_mean"] - df_ave["y_std"],
        fill="tonexty"
        )
    trace.update(**shared_scatter_props)
    trace.update(**unique_scatter_props)
    traces.append(trace)
    # #########################################################################
    # out_data_dict["trace"] = trace
    out_data_dict["trace"] = traces
    out_data_dict["df_ave"] = df_ave

    return(out_data_dict)
    #__|

data = []

# + {"active": ""}
#
#
#
#
# -

# # Random | w/ Duplicates

# +
#############################################################################
color_i = "rgb(100,100,100,0.5)"
out_data_dict_i = process_data(
    subdir="random_True",
    unique_scatter_props=dict(
        name="random w/ dupl",
        marker=dict(color=color_i),
        error_y=dict(
            color=color_i,
            ),
        ),
    shared_scatter_props=shared_scatter_props,
    ALPerf_account_duplicates=True,
    top_ids_to_track=top_ids_to_track,
    )
trace_i = out_data_dict_i["trace"]
# data.append(trace_i)
data.extend(trace_i)

df_perf = out_data_dict_i["df_perf"]
# df_perf
# -

# # Random | w/o Duplicates

# +
# # #############################################################################
# # color_i = "grey"
# color_i = "rgb(60,120,100,0.5)"
# out_data_dict_i = process_data(
#     subdir="random_False",
#     unique_scatter_props=dict(
#         name="random w/o dupl",
#         marker=dict(color=color_i),
#         error_y=dict(
#             color=color_i,
#             ),
#         ),
#     shared_scatter_props=shared_scatter_props,
#     ALPerf_account_duplicates=False,
#     top_ids_to_track=top_ids_to_track,
#     )
# trace_i = out_data_dict_i["trace"]
# # data.append(trace_i)
# data.extend(trace_i)

# df_perf = out_data_dict_i["df_perf"]
# # df_perf
# -

# # GP-UCB | w/ Duplicates

# +
# #############################################################################
# color_i = "red"
color_i = "rgba(255,0,0,0.5)"

out_data_dict_i = process_data(
    subdir="gp_ucb_True",
    unique_scatter_props=dict(
        name="gp_ucb w/ dupl",
        marker=dict(color=color_i),

        error_y=dict(
            color=color_i,
            ),

        ),
    shared_scatter_props=shared_scatter_props,
    ALPerf_account_duplicates=True,
    top_ids_to_track=top_ids_to_track,
    )
trace_i = out_data_dict_i["trace"]
# data.append(trace_i)
data.extend(trace_i)

df_perf = out_data_dict_i["df_perf"]
# df_perf

# +
# df_perf
# pifehohu
# geheneva
# nisoponi
# -

# # GP-UCB | w/o Duplicates

# +
# # #############################################################################
# color_i = "orange"
# out_data_dict_i = process_data(
#     subdir="gp_ucb_False",
#     unique_scatter_props=dict(
#         name="gp_ucb w/o dupl",
#         marker=dict(color=color_i),
#         error_y=dict(
#             color=color_i,
#             ),
#         ),
#     shared_scatter_props=shared_scatter_props,
#     ALPerf_account_duplicates=True,
#     top_ids_to_track=top_ids_to_track,
#     )
# trace_i = out_data_dict_i["trace"]
# # data.append(trace_i)
# data.extend(trace_i)

# df_perf = out_data_dict_i["df_perf"]
# # df_perf
# -

# # Plotting

# +
shared_shape_dict = dict(
    xref="x", yref="y",
    type="line",
    line=dict(
        color="black",
        width=2,
        dash="dot",
        ),
    )

shapes = tuple([

    go.layout.Shape(
        x0=50,  y0=-1,
        x1=50, y1=7,

        **shared_shape_dict),

    go.layout.Shape(
        x0=0,  y0=7,
        x1=50, y1=7,

        **shared_shape_dict),


    go.layout.Shape(
        x0=157,  y0=-1,
        x1=157, y1=7,

        **shared_shape_dict),

    go.layout.Shape(
        x0=0,  y0=7,
        x1=157, y1=7,

        **shared_shape_dict),


    ])

# +
layout["height"] = 37.795275591 * 7.12
layout["width"] = 37.795275591 * 6.3

layout["paper_bgcolor"] = "rgba(0,0,0,0)"
layout["plot_bgcolor"] = "rgba(0,0,0,0)"

fig = go.Figure(data=data, layout=layout)

fig.layout.update(
    shapes=shapes,
    xaxis=dict(range=[-0.8, 250]),
    yaxis=dict(range=[-0.3, 10.6]),
    )

fig = my_plotly_plot(
    figure=fig,
    plot_name=stoich_i + "_" + "al_performance",
    write_html=True,
    write_png=False,
    png_scale=6.0,
    write_pdf=True,
    write_svg=True,
    )


fig.layout.update(paper_bgcolor="white")
fig.show()
# -

# Pickling data ######################################################
import os; import pickle
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)
with open(os.path.join(directory, stoich_i + "_" + "fig_al_perf.pickle"), "wb") as fle:
    pickle.dump(fig, fle)
# #####################################################################

# +
# go.Layout?
# go.layout.Scene?
# go.layout.Shape?

# + {"active": ""}
#
#
#
#
#

# + {"jupyter": {"source_hidden": true}}
# shapes = tuple([

#     go.layout.Shape(
#         type="line",
#         xref="x",
#         yref="y",
#         x0=50,
#         y0=0,
#         x1=100,
#         y1=6,
#         line=dict(
#             color="DarkOrange",
#             width=3,
#             ),
#         ),

#     go.layout.Shape(
#         type="line",
#         xref="x",
#         yref="y",
#         x0=20,
#         y0=1,
#         x1=120,
#         y1=8,
#         line=dict(
#             color="DarkOrange",
#             width=3,
#             ),
#         ),

#     ])

# fig.layout.update(shapes=shapes)

# fig

# + {"jupyter": {"source_hidden": true}}
# # AL_i = 
# keys = list(out_data_dict_i["AL_dict"].keys())

# for key in keys:
#     # AL_i = out_data_dict_i["AL_dict"][key].al_gen_dict[49]
#     AL = out_data_dict_i["AL_dict"][key]
#     al_gen_dict = AL.al_gen_dict

#     print(len(AL.duplicate_ids))
    
#     if len(AL.duplicate_ids) == 82:
#         print(AL.duplicate_ids)
#         break
#     # print(key)

#     # print(len(al_gen_dict.keys()))
#     # print("")



# AL.duplicate_ids

# # AL_i = al_gen_dict[49]


# # .al_gen_dict[25]
# # out_data_dict_i["AL_dict"][0].al_gen_dict[25]
# # model = AL_i.model
# # model_i = model[
# #     (model["acquired"] == True) & \
# #     (model["duplicate"] == False)
# #     ]
# # # AL_i.indices_that_are_duplicates
# # model_i.sort_values("y_real").iloc[0:50].index.tolist()

# + {"jupyter": {"source_hidden": true}}
# layout["height"] = 37.795275591 * 7.05
# layout["height"] = 37.795275591 * 7.07
# layout["height"] = 37.795275591 * 7.055
# layout["height"] = 37.795275591 * 7.052
# layout["height"] = 37.795275591 * 7.048
# layout["height"] = 37.795275591 * 7.035
# layout["height"] = 37.795275591 * 7.041
# layout["height"] = 37.795275591 * 7.043  # too short
# layout["height"] = 37.795275591 * 7.045  # short
# layout["height"] = 37.795275591 * 7.047  # short
# layout["height"] = 37.795275591 * 7.052  # short
# layout["height"] = 37.795275591 * 7.09  # long
# layout["height"] = 37.795275591 * 7.08 # short
# layout["height"] = 37.795275591 * 7.085  # long
# layout["height"] = 37.795275591 * 7.082  # long
# layout["height"] = 37.795275591 * 7.07  # long
# layout["height"] = 37.795275591 * 7.06  # long
# layout["height"] = 37.795275591 * 7.04  # short
# layout["height"] = 37.795275591 * 7.045 # short
# layout["height"] = 37.795275591 * 7.049 # short
# layout["height"] = 37.795275591 * 7.051
# layout["height"] = 37.795275591 * 7.053 #short
# layout["height"] = 37.795275591 * 7.055 #l
# layout["height"] = 37.795275591 * 7.054
# layout["height"] = 37.795275591 * 7.053 #l
# layout["height"] = 37.795275591 * 7.052 #l
# layout["height"] = 37.795275591 * 7.051 # s
# layout["height"] = 37.795275591 * 7.056 # s
# layout["height"] = 37.795275591 * 7.06
# layout["height"] = 37.795275591 * 7.066
# layout["height"] = 37.795275591 * 7.07
# layout["height"] = 37.795275591 * 7.085
# layout["height"] = 37.795275591 * 7.092
# layout["height"] = 37.795275591 * 7.098
# layout["height"] = 37.795275591 * 7.2
# layout["height"] = 37.795275591 * 7.14

# layout["width"] = 37.9 * 5
# layout["width"] = 37.795275591 * 20
# layout["width"] = 37.795275591 * 10
# layout["width"] = 37.795275591 * 8
# layout["width"] = 37.795275591 * 6
# layout["width"] = 37.795275591 * 6.5

# + {"jupyter": {"source_hidden": true}}
# layout["height"] *= 2
# layout["width"] *= 4

# layout["paper_bgcolor"] = "white"
# layout["plot_bgcolor"] = "white"

# layout["showlegend"] = True

# fig.update_layout(layout)
# fig.show()

# + {"jupyter": {"source_hidden": true}}
# AL_i = out_data_dict_i["AL_dict"][0].al_gen_dict[25]
# # out_data_dict_i["AL_dict"][0].al_gen_dict[25]
# model = AL_i.model
# model_i = model[
#     (model["acquired"] == True) & \
#     (model["duplicate"] == False)
#     ]
# # AL_i.indices_that_are_duplicates
# model_i.sort_values("y_real").iloc[0:10].index.tolist()

# + {"jupyter": {"source_hidden": true}}
# # #############################################################################
# out_data_dict_i = process_data(
#     subdir="gp_ucb_True",
#     unique_scatter_props=dict(
#         name="gp_ucb w dupl",
#         marker=dict(color="red"),
#         ),
#     shared_scatter_props=shared_scatter_props,
#     ALPerf_account_duplicates=True,
#     top_ids_to_track=top_ids_to_track,
#     )

# trace_i = out_data_dict_i["trace"]
# data.append(trace_i)

# df_perf = out_data_dict_i["df_perf"]

# + {"jupyter": {"source_hidden": true}}
# df_perf

# df_m = df_perf


# columns_list = list(df_m.columns.levels[0])
# col = df_m.loc[: , columns_list[0]:columns_list[-1]]

# col


# y_mean = col.mean(axis=1)
# y_std = col.std(axis=1)

# y_mean
# y_std
