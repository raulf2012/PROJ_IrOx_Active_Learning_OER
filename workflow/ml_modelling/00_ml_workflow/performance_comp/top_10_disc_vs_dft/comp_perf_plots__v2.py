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

# #########################################################
# Local Import ############################################
from layout import layout
from inputs import stoich_i

# +
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow/191102_new_workflow"))
from al_data import al_data_files_dict

files_list_gp_ucb = al_data_files_dict[stoich_i]["files_list_gp_ucb"]
files_list_random = al_data_files_dict[stoich_i]["files_list_random"]

# +
# if stoich_i == "AB2":
#     files_list_gp_ucb = al_data_files_dict[stoich_i]["files_list_ab2_gp_ucb"]
#     files_list_random = al_data_files_dict[stoich_i]["files_list_ab2_random"]
# elif stoich_i == "AB3":
#     files_list_gp_ucb = al_data_files_dict[stoich_i]["files_list_ab3_gp_ucb"]
#     files_list_random = al_data_files_dict[stoich_i]["files_list_ab3_random"]
# -

files_list_random[0]

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


from inputs import top_ids_to_track_ab2, top_ids_to_track_ab3
# -

plot_guidlines = False

# + Collapsed="false" jupyter={"outputs_hidden": false}
if stoich_i == "AB2":
    top_ids_to_track = top_ids_to_track_ab2
elif stoich_i == "AB3":
    top_ids_to_track = top_ids_to_track_ab3
else:
    print("ISDJIFSDJI")


# + Collapsed="false"
def process_data(
    subdir=None,
    shared_scatter_props=None,
    unique_scatter_props=None,
    ALPerf_account_duplicates=True,
    top_ids_to_track=None,
    files_list=None,
    color2=None,
    ):
    """
    """
    # | - process_data
    out_data_dict = dict()

    # #############################################################################
    if files_list is not None:
        tmp = 42
    else:
        files_list = os.listdir(
            os.path.join(
                # dir_i,
                data_path_root,
                # "out_data",
                subdir))
        files_list = [i for i in files_list if "pickle" in i]
        files_list = [i for i in files_list if "AL_" in i]

    # print(files_list)

    data_dict = dict()
    for file_i in files_list:
        # #########################################################################
        num = file_i.split("_")[-1].split(".")[0]

        file_path_i =os.path.join(
            data_path_root, subdir, file_i)
        # COMBAK
        # with open(file_path_i, "rb") as fle:
        with open(file_i, "rb") as fle:
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

    # Adding 0 to trace
    df_ave.loc[0] = [0, 0]
    df_ave = df_ave.sort_index()

    dx = df_ave.index.values[-1] - df_ave.index.values[-2]
    # last_data_point_ind = df_ave.index.values[-1] + dx
    last_data_point_ind = df_ave.index.values[-1] + 50
    df_ave.loc[last_data_point_ind] = [10, 0]

    # df_ave.loc[260] = [10, 0]
    df_ave = df_ave.sort_index()


    out_data_dict["df_perf"] = df_m

    traces = []
    # #############################################################################
    trace = go.Scatter(
        x=df_ave.index.tolist(),
        y=df_ave["y_mean"],
        line=dict(
            width=1.,
            ),
        )
    trace.update(**shared_scatter_props)
    trace.update(**unique_scatter_props)
    traces.append(trace)
    # #########################################################################
    trace = go.Scatter(
        x=df_ave.index.tolist(),
        y=df_ave["y_mean"] + df_ave["y_std"],
        line=dict(
            width=0.5,
            ),
        )
    trace.update(**shared_scatter_props)
    trace.update(**unique_scatter_props)
    traces.append(trace)

    trace = go.Scatter(
        x=df_ave.index.tolist(),
        y=df_ave["y_mean"] - df_ave["y_std"],
        fill="tonexty",
        line=dict(
            width=0.5,
            # color="red",
            ),
        )
    trace.update(**shared_scatter_props)
    trace.update(**unique_scatter_props)
    traces.append(trace)
    # #########################################################################

    trace = go.Scatter(
        x=df_ave.index.tolist(),
        y=df_ave["y_mean"] - df_ave["y_std"],
        line=dict(
            width=0.5,
            color=color2,
            ),
        )
    trace.update(**shared_scatter_props)
    trace.update(**unique_scatter_props)
    traces.append(trace)

    trace = go.Scatter(
        x=df_ave.index.tolist(),
        y=df_ave["y_mean"] + df_ave["y_std"],
        line=dict(
            width=0.5,
            color=color2,
            ),
        )
    trace.update(**shared_scatter_props)
    trace.update(**unique_scatter_props)
    traces.append(trace)

    # out_data_dict["trace"] = trace
    out_data_dict["trace"] = traces
    out_data_dict["df_ave"] = df_ave

    return(out_data_dict)
    #__|

# + Collapsed="false" jupyter={"outputs_hidden": false}
data = []

# + Collapsed="false" active=""
#
#
#
#

# + [markdown] Collapsed="false"
# # Random | w/ Duplicates

# + Collapsed="false" jupyter={"outputs_hidden": false}
#############################################################################
color_i = "rgb(100,100,100,0.5)"
out_data_dict_i = process_data(
    # files_list=files_list,
    files_list=files_list_random,

    # from inputs import files_list_gp_ucb, files_list_random

    subdir="random_True",
    unique_scatter_props=dict(
        name="random w/ dupl",
        marker=dict(color=color_i),
        error_y=dict(
            # color=color_i,
            # color="red",
            ),
        ),
    shared_scatter_props=shared_scatter_props,
    ALPerf_account_duplicates=True,
    top_ids_to_track=top_ids_to_track,
    color2="rgb(100,100,100,1.)",
    )
trace_i = out_data_dict_i["trace"]
data.extend(trace_i)

df_perf_random = out_data_dict_i["df_perf"]
# df_perf
# -

print("num of runs random:", "\n", df_perf_random.shape[1])

# +
df_ave = out_data_dict_i["df_ave"]

x_interc0 = np.interp(
    num_disc,
    df_ave.y_mean.tolist(),
    df_ave.index.tolist(),
    )

# +
# assert False

# + [markdown] Collapsed="false"
# # Random | w/o Duplicates

# + Collapsed="false" jupyter={}
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

# + [markdown] Collapsed="false"
# # GP-UCB | w/ Duplicates

# + Collapsed="false" jupyter={"outputs_hidden": false}
# #############################################################################
# color_i = "red"
color_i = "rgba(0,100,255,0.5)"

out_data_dict_i = process_data(
    # files_list=files_list,
    files_list=files_list_gp_ucb,
    subdir="gp_ucb_True/01_attempt",
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
    # color_i = "rgba(0,100,255,0.5)"
    color2="rgba(0,100,255,1.)",
    )
trace_i = out_data_dict_i["trace"]
data.extend(trace_i)

df_perf_gpucb = out_data_dict_i["df_perf"]
# df_perf
# -

print("num of runs GP-UCB:", "\n", df_perf_gpucb.shape[1])

# +
df_ave = out_data_dict_i["df_ave"]

x_interc1 = np.interp(
    num_disc,
    df_ave.y_mean.tolist(),
    df_ave.index.tolist(),
    )

# + [markdown] Collapsed="false"
# # GP-UCB | w/o Duplicates

# + Collapsed="false" jupyter={}
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

# # Red guide-lines

# + Collapsed="false" jupyter={"outputs_hidden": false}
shared_shape_dict = dict(
    xref="x", yref="y",
    type="line",
    line=dict(
        color="red",
        width=1.5,
        dash="dot",
        ),
    )

shapes = tuple([

    go.layout.Shape(
        x0=x_interc1,  y0=-1,
        x1=x_interc1, y1=num_disc,
        **shared_shape_dict),

    go.layout.Shape(
        x0=0,  y0=num_disc,
        x1=x_interc1, y1=num_disc,
        **shared_shape_dict),


    go.layout.Shape(
        x0=x_interc0,  y0=-1,
        x1=x_interc0, y1=num_disc,
        **shared_shape_dict),

    go.layout.Shape(
        x0=0,  y0=num_disc,
        x1=x_interc0, y1=num_disc,
        **shared_shape_dict),

    ])

# + [markdown] Collapsed="false"
# # Plotting

# + Collapsed="false"
layout["height"] = 37.795275591 * 7.12
layout["width"] = 37.795275591 * 6.3

layout["paper_bgcolor"] = "rgba(0,0,0,0)"
layout["plot_bgcolor"] = "rgba(0,0,0,0)"

fig = go.Figure(data=data, layout=layout)

if stoich_i == "AB2":
    x_range = [-0.8, 470]
    y_range = [-0.3, 10.6]
elif stoich_i == "AB3":
    x_range = [-0.8, 250]
    y_range = [-0.3, 10.6]

if plot_guidlines:
    shapes = shapes
else:
    shapes = None

fig.layout.update(
    shapes=shapes,
    # xaxis=dict(range=[-0.8, 250]),
    # yaxis=dict(range=[-0.3, 10.6]),
    xaxis=dict(range=x_range),
    yaxis=dict(range=y_range),
    )

# fig = my_plotly_plot(
my_plotly_plot(
    figure=fig,
    plot_name=stoich_i + "_" + "al_performance",
    write_html=True,
    write_png=False,
    png_scale=6.0,
    write_pdf=True,
    write_svg=False,
    try_orca_write=True,
    )


fig.layout.update(paper_bgcolor="white")
# fig.show()

tmp = 42
# -

fig.show()

# +
# fig.layout.update(dict(
#     height=500,
#     width=600,
#     showlegend=True,
#     ))

# fig.show()

# + Collapsed="false" jupyter={"outputs_hidden": false}
# Pickling data ######################################################
import os; import pickle
directory = "out_data"
if not os.path.exists(directory):
    os.makedirs(directory)

# #####################################################################
with open(os.path.join(directory, stoich_i + "_" + "fig_al_perf.pickle"), "wb") as fle:
    pickle.dump(fig, fle)
# #####################################################################

# #####################################################################
with open(os.path.join(directory, stoich_i + "_" + "df_random.pickle"), "wb") as fle:
    pickle.dump(df_perf_random, fle)
with open(os.path.join(directory, stoich_i + "_" + "df_gbucb.pickle"), "wb") as fle:
    pickle.dump(df_perf_gpucb, fle)
# #####################################################################

# +
# df_perf_random

# + Collapsed="false" active=""
#
#
#
#
#

# + jupyter={}
# data_i = []
# df_i = df_perf_gpucb
# for i_cnt, col_i in enumerate(df_i.columns):
#     trace_i = go.Scatter(
#         x=df_i[col_i].index,
#         y=df_i[col_i].values,
#         # name=names_list[i_cnt],
#         name=col_i[0],
#         )
#     data_i.append(trace_i)
    
# fig = go.Figure(data=data_i)
# # fig.show()

# + jupyter={}
# data_i = []
# df_i = df_perf_random
# for i_cnt, col_i in enumerate(df_i.columns):
#     name_i = col_i[0]

#     trace_i = go.Scatter(
#         x=df_i[col_i].index,
#         y=df_i[col_i].values,
#         # name=col_i[0],
#         name=name_i,
#         )
#     data_i.append(trace_i)

# fig = go.Figure(data=data_i)
# # fig.show()

# + jupyter={}
# # df_ave

# df_ave.loc[0] = [0, 0]
# df_ave = df_ave.sort_index()

# df_ave


# df_perf
# pifehohu
# geheneva
# nisoponi

# if stoich_i == "AB3":
#     path_i = os.path.join(
#         os.environ["PROJ_irox"],
#         "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/00_abx_al_runs/out_data/AB3/gp_ucb_True",

#         # "01_attempt/AL_geheneva.pickle",
#         # "01_attempt/AL_pifehohu.pickle",
        
#         # NEW RUNS
#         # "TEST_AL_2_fugunefo.pickle",
#         "TEST_AL_2_seruladi.pickle",
#         )

# pre_path = os.path.join(
#     os.environ["PROJ_irox"],
#     "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/00_abx_al_runs/out_data/AB3/gp_ucb_True",
#     )

# files_list = [
#     pre_path + "/01_attempt/AL_geheneva.pickle",
#     pre_path + "/01_attempt/AL_nisoponi.pickle",
#     pre_path + "/01_attempt/AL_pifehohu.pickle",
#     pre_path + "/01_attempt/AL_suturomo.pickle",
#     pre_path + "/01_attempt/AL_vobifoko.pickle",
#     pre_path + "/TEST_AL_masahiti.pickle",

#     # pre_path + "/TEST_AL_2_devehowo.pickle",  # Not a good run for some reason
#     pre_path + "/TEST_AL_2_fugunefo.pickle",
#     pre_path + "/TEST_AL_2_hilerika.pickle",
#     pre_path + "/TEST_AL_2_pomogobu.pickle",
#     pre_path + "/TEST_AL_2_seruladi.pickle",

    
#     pre_path + "/TEST_AL_3_bikufupi.pickle",
#     pre_path + "/TEST_AL_3_dakubiku.pickle",
#     pre_path + "/TEST_AL_3_duloputo.pickle",
#     pre_path + "/TEST_AL_3_fahovara.pickle",
#     pre_path + "/TEST_AL_3_fulomoto.pickle",
#     pre_path + "/TEST_AL_3_laburike.pickle",
#     pre_path + "/TEST_AL_3_libapidi.pickle",
#     pre_path + "/TEST_AL_3_nikimido.pickle",
#     pre_path + "/TEST_AL_3_raluduhu.pickle",
#     pre_path + "/TEST_AL_3_supemono.pickle",

#     pre_path + "/TEST_AL_4beradeka.pickle",
#     pre_path + "/TEST_AL_4buruduwe.pickle",
#     pre_path + "/TEST_AL_4degekoku.pickle",
#     pre_path + "/TEST_AL_4deromeru.pickle",
#     pre_path + "/TEST_AL_4forafago.pickle",
#     pre_path + "/TEST_AL_4hefehepa.pickle",
#     pre_path + "/TEST_AL_4kaveboma.pickle",
#     pre_path + "/TEST_AL_4kihalage.pickle",
#     pre_path + "/TEST_AL_4lidirope.pickle",
#     pre_path + "/TEST_AL_4mebetige.pickle",
#     pre_path + "/TEST_AL_4megimodi.pickle",
#     pre_path + "/TEST_AL_4mohomato.pickle",
#     pre_path + "/TEST_AL_4moponuso.pickle",
#     pre_path + "/TEST_AL_4mukigapa.pickle",
#     pre_path + "/TEST_AL_4nekumuno.pickle",
#     pre_path + "/TEST_AL_4nipagula.pickle",
#     pre_path + "/TEST_AL_4nupetofi.pickle",
#     pre_path + "/TEST_AL_4rifibume.pickle",
#     pre_path + "/TEST_AL_4somageho.pickle",
#     pre_path + "/TEST_AL_4wolewoba.pickle",

#     pre_path + "/TEST_AL_5bofufada.pickle",
#     pre_path + "/TEST_AL_5delepaku.pickle",
#     pre_path + "/TEST_AL_5derohebi.pickle",
#     pre_path + "/TEST_AL_5dotesiga.pickle",
#     pre_path + "/TEST_AL_5dumokeru.pickle",
#     pre_path + "/TEST_AL_5fapehudo.pickle",
#     pre_path + "/TEST_AL_5fenumanu.pickle",
#     pre_path + "/TEST_AL_5gomememu.pickle",
#     pre_path + "/TEST_AL_5gulumobu.pickle",
#     pre_path + "/TEST_AL_5huwihime.pickle",
#     pre_path + "/TEST_AL_5kanototo.pickle",
#     pre_path + "/TEST_AL_5kibapeto.pickle",
#     pre_path + "/TEST_AL_5kifomimi.pickle",
#     pre_path + "/TEST_AL_5rukeraku.pickle",
#     pre_path + "/TEST_AL_5sipeteni.pickle",
#     pre_path + "/TEST_AL_5sokepefo.pickle",
#     pre_path + "/TEST_AL_5vikowagu.pickle",
#     pre_path + "/TEST_AL_5volibavo.pickle",
#     pre_path + "/TEST_AL_5wigidipu.pickle",
#     pre_path + "/TEST_AL_5wolowewu.pickle",
    
#     pre_path + "/TEST_AL_6_dodemuho.pickle",
#     pre_path + "/TEST_AL_6_fisopova.pickle",
#     pre_path + "/TEST_AL_6_gelabere.pickle",
#     pre_path + "/TEST_AL_6_gelenuni.pickle",
#     pre_path + "/TEST_AL_6_haligagu.pickle",
#     pre_path + "/TEST_AL_6_higusare.pickle",
#     pre_path + "/TEST_AL_6_kagesinu.pickle",
#     pre_path + "/TEST_AL_6_liwuderu.pickle",
#     pre_path + "/TEST_AL_6_lopukipu.pickle",
#     pre_path + "/TEST_AL_6_nubinada.pickle",
#     pre_path + "/TEST_AL_6_pitumimi.pickle",
#     pre_path + "/TEST_AL_6_pudapepi.pickle",
#     pre_path + "/TEST_AL_6_rividisi.pickle",
#     pre_path + "/TEST_AL_6_rovomama.pickle",
#     pre_path + "/TEST_AL_6_sovawafo.pickle",
#     pre_path + "/TEST_AL_6_takoliwo.pickle",
#     pre_path + "/TEST_AL_6_tehileme.pickle",
#     pre_path + "/TEST_AL_6_tibuduka.pickle",
#     pre_path + "/TEST_AL_6_warusiha.pickle",
#     pre_path + "/TEST_AL_6_wevuwofu.pickle",

#     ]

# pre_path = os.path.join(
#     os.environ["PROJ_irox"],
#     "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/00_abx_al_runs/out_data/AB3/random_True",
#     # "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/00_abx_al_runs/out_data/AB3/gp_ucb_True",
#     )

# files_list = [
#     pre_path + "/AL_duhakiro.pickle",
#     pre_path + "/AL_firegohi.pickle",
#     pre_path + "/AL_givohegu.pickle",
#     pre_path + "/AL_tiweluku.pickle",
#     pre_path + "/AL_vevenuwa.pickle",
#     pre_path + "/TEST_AL_2_fibataha.pickle",
#     pre_path + "/TEST_AL_2_kitagego.pickle",
#     pre_path + "/TEST_AL_3_bafigika.pickle",
#     pre_path + "/TEST_AL_3_bafomepo.pickle",
#     pre_path + "/TEST_AL_3_besurogi.pickle",
#     pre_path + "/TEST_AL_3_dotupibu.pickle",
#     pre_path + "/TEST_AL_3_fuliguso.pickle",
#     pre_path + "/TEST_AL_3_gosirinu.pickle",
#     pre_path + "/TEST_AL_3_gowutoga.pickle",
#     pre_path + "/TEST_AL_3_hiroguwe.pickle",
#     pre_path + "/TEST_AL_3_kavasaki.pickle",
#     pre_path + "/TEST_AL_3_mutefoti.pickle",
#     pre_path + "/TEST_AL_3_nivalula.pickle",
#     pre_path + "/TEST_AL_3_pasimuha.pickle",
#     pre_path + "/TEST_AL_3_sukiwalo.pickle",
#     pre_path + "/TEST_AL_3_temonofo.pickle",
#     pre_path + "/TEST_AL_3_vesudewa.pickle",
#     pre_path + "/TEST_AL_3_vuwugupi.pickle",
#     pre_path + "/TEST_AL_3_walipebi.pickle",
#     pre_path + "/TEST_AL_3_wetipotu.pickle",
#     pre_path + "/TEST_AL_3_wusabupa.pickle",
#     pre_path + "/TEST_AL_3_wutonovi.pickle",

#     pre_path + "/TEST_AL_4_benegeka.pickle",
#     pre_path + "/TEST_AL_4_dehebiko.pickle",
#     pre_path + "/TEST_AL_4_dinisefa.pickle",
#     pre_path + "/TEST_AL_4_fefefigi.pickle",
#     pre_path + "/TEST_AL_4_fefesama.pickle",
#     pre_path + "/TEST_AL_4_fivokito.pickle",
#     pre_path + "/TEST_AL_4_fuwasufi.pickle",
#     pre_path + "/TEST_AL_4_gekuporu.pickle",
#     pre_path + "/TEST_AL_4_gepapeba.pickle",
#     pre_path + "/TEST_AL_4_gerisiwe.pickle",
#     pre_path + "/TEST_AL_4_goderiwo.pickle",
#     pre_path + "/TEST_AL_4_hofavasi.pickle",
#     pre_path + "/TEST_AL_4_kavadosu.pickle",
#     pre_path + "/TEST_AL_4_kudadega.pickle",
#     pre_path + "/TEST_AL_4_masufika.pickle",
#     pre_path + "/TEST_AL_4_menireve.pickle",
#     pre_path + "/TEST_AL_4_metekovo.pickle",
#     pre_path + "/TEST_AL_4_mibunova.pickle",
#     pre_path + "/TEST_AL_4_milesumi.pickle",
#     pre_path + "/TEST_AL_4_nepubene.pickle",
#     pre_path + "/TEST_AL_4_nuvopeki.pickle",
#     pre_path + "/TEST_AL_4_petosuso.pickle",
#     pre_path + "/TEST_AL_4_pokikugi.pickle",
#     pre_path + "/TEST_AL_4_ragifipa.pickle",
#     pre_path + "/TEST_AL_4_rodadibi.pickle",
#     pre_path + "/TEST_AL_4_rovukuma.pickle",
#     pre_path + "/TEST_AL_4_sanaruri.pickle",
#     pre_path + "/TEST_AL_4_sapaveme.pickle",
#     pre_path + "/TEST_AL_4_sawuhewe.pickle",
#     pre_path + "/TEST_AL_4_sifisulu.pickle",
#     pre_path + "/TEST_AL_4_sifobuwe.pickle",
#     pre_path + "/TEST_AL_4_suhegope.pickle",
#     pre_path + "/TEST_AL_4_teruwufo.pickle",
#     pre_path + "/TEST_AL_4_tikeluvi.pickle",
#     pre_path + "/TEST_AL_4_vabesipo.pickle",
#     pre_path + "/TEST_AL_4_visisese.pickle",
#     pre_path + "/TEST_AL_4_wemawumu.pickle",
#     pre_path + "/TEST_AL_4_wibemafu.pickle",
#     pre_path + "/TEST_AL_4_wobumone.pickle",
#     pre_path + "/TEST_AL_4_wohaguha.pickle",

#     ]

# + jupyter={}
# df_perf_random


# # Simple Plotly Plot
# import plotly.graph_objs as go
# trace = go.Scatter(
#     x=df_ave.index.values,
#     y=df_ave.y_mean)
# data = [trace]
# fig = go.Figure(data=data)
# fig.show()

# + jupyter={}
# df_ave[
#     (df_ave.index.values < 300) & \
#     (df_ave.index.values > 200)
#     ]
