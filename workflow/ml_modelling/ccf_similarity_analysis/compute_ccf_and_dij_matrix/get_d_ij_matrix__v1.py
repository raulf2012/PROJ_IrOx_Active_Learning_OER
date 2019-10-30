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

import time
t0 = time.time()

# + {"jupyter": {"source_hidden": true}}
import os
import sys
import copy

import pickle

import numpy as np

import pandas as pd

from ase.visualize import view

# Plotly
import chart_studio.plotly as py
# import plotly.graph_objs as go
import plotly.graph_objects as go

from StructurePrototypeAnalysisPackage.ccf import (
    struc2ccf,
    cal_ccf_d,
    cal_inter_atomic_d,
    d2ccf,
    weight_f,
    pearson_cc,
    gaussian_f,
    element_tag,
    cell_range,
    count_atoms_dict,
    )

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import (
    static_irox_structures_path,
    bulk_dft_data_path,
    )
# -

# # Script Inputs

# +
r_cut_off = 10.

r_vector = np.arange(1, 10, 0.02)

verbose = True
# -

# ## Read Data

# +
# with open("out_data/df_ccf_test.pickle", "rb") as fle:
with open("out_data/df_ccf.pickle", "rb") as fle:
    df_ccf = pickle.load(fle)
# df_ccf = df_ccf[0:40]

# #############################################################################
with open("out_data/df_d_ij_all.pickle", "rb") as fle:
# with open("out_data/df_d_ij_all_temp.pickle", "rb") as fle:
    df_d_comp_prev = pickle.load(fle)

# +
print("df_ccf.shape:", df_ccf.shape)

print("df_d_comp_prev.shape:", df_d_comp_prev.shape)
# -

# # Process data (Create D_ij matrix)

# + {"jupyter": {"outputs_hidden": true}}
df = df_ccf

result = np.zeros((len(df), len(df)))
for i_cnt, (row_name_i, row_i) in enumerate(df.iterrows()):
    print(str(i_cnt).zfill(4), 75 * "#")
    
    row_in_prev_data = False

    for j_cnt, (row_name_j, row_j) in enumerate(df.iterrows()):
        if i_cnt == j_cnt:
            continue


        index_in_cols = row_name_j in df_d_comp_prev.columns
        index_in_rows = row_name_i in df_d_comp_prev.index
        if index_in_cols and index_in_rows:
            if True:
#             try:
                d_ij = df_d_comp_prev.loc[row_name_i, row_name_j]

                if verbose:
                    row_in_prev_data = True
#                     print("Parsed d_ij from previous data")
        else:
            if verbose:
                tmp = 42
                # print("Computing d_ij from scratch")

            d_ij = cal_ccf_d(
                row_i["ccf"],
                row_j["ccf"])
            # d_ij = 42.

        result[i_cnt][j_cnt] = d_ij

    if row_in_prev_data:
        print("This row is the previous data (GOOD)")
df_d_comp = pd.DataFrame(
    result,
    index=df.index,
    columns=df.index)


# Pickling data ######################################################
import os; import pickle
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)


# with open(os.path.join(directory, "df_d_ij_all_temp.pickle"), "wb") as fle:
with open(os.path.join(directory, "df_d_ij_all.pickle"), "wb") as fle:
    pickle.dump(df_d_comp, fle)
# #####################################################################
# -

# path_i = os.path.join("out_data", "df_d_ij_all_temp.pickle")
path_i = os.path.join("out_data", "df_d_ij_all.pickle")
with open(path_i, "rb") as fle:
    df_d_comp = pickle.load(fle)

# + {"active": ""}
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

time_elapsed = time.time() - t0
print(
    "Notebook took ",
    round(time_elapsed / 60, 4),
    "minutes to finish")

# + {"active": ""}
#
#
#
# -

'vfxezkmsv2' in df_d_comp.index

df_d_comp.describe()

"bsv4nex29l" in df_d_comp.index

# + {"jupyter": {"source_hidden": true}, "toc-hr-collapsed": true}
# max_d = df_d_comp.max().max()
# min_d = df_d_comp.min().min()
# data_range_z = (max_d - min_d)
# print("data_range_z: ", data_range_z)

# # d_thresh = 0.01
# # d_thresh = 0.05
# d_thresh = 0.075
# d_thresh_stan = d_thresh * (1 / data_range_z)

# # #############################################################################
# # #############################################################################
# trace_i = go.Heatmap(
#     z=df_d_comp.values,
#     x=df_d_comp.index.tolist(),
#     y=df_d_comp.index.tolist(),
#     colorscale=[
#         [0.0, "black"],
#         [0.000001, "red"],
#         [d_thresh_stan, "pink"],
#         [d_thresh_stan + 0.000001, "rgb(220,220,220)"],
#         [1.0, "white"],
# #         [1.2, "blue"],
#         ]

#     )
# data = [trace_i]

# layout = go.Layout(width=1000, height=1000)
# fig = go.Figure(data=data, layout=layout)
# fig

## Calculate d_ij

# df = df_ccf

# result = np.zeros((len(df), len(df)))
# for i_cnt, (row_name_i, row_i) in enumerate(df.iterrows()):
#     print(str(i_cnt).zfill(4), 75 * "#")
#     for j_cnt, (row_name_j, row_j) in enumerate(df.iterrows()):
#         if i_cnt == j_cnt:
#             continue
# #         if j_cnt < i_cnt:
# #             continue

#         d_ij = cal_ccf_d(
#             row_i["ccf"],
#             row_j["ccf"])

#         result[i_cnt][j_cnt] = d_ij

# df_d_comp = pd.DataFrame(
#     result,
#     index=df.index,
#     columns=df.index,
#     )

# # df_d_comp

# # Pickling data ######################################################
# import os; import pickle
# directory = "out_data"
# if not os.path.exists(directory): os.makedirs(directory)
# with open(os.path.join(directory, "df_d_ij_static_irox_structures.pickle"), "wb") as fle:
#     pickle.dump(df_d_comp, fle)
# # #####################################################################

# path_i = os.path.join("out_data", "df_d_ij_static_irox_structures.pickle")
# with open(path_i, "rb") as fle:
#     df_d_comp = pickle.load(fle)

# def get_ccf_df(ccf, r_cut_off, r_vector):
#     """
#     """
#     df_i = pd.DataFrame()
#     for key, value in ccf.items():
#         df_i[key] = value

#     df_i = df_i.set_index(r_vector)
    
#     return(df_i)

# with open("out_data/df_dft_irox_structures_ccf.pickle", "rb") as fle:
#     df_ccf = pickle.load(fle)
# # df_ccf = df_ccf[0:80]

# Process Static (Pre-optimized) IrOx Structures

## Read Data

# # with open("out_data/df_static_irox_structures_ccf.pickle", "rb") as fle:
# with open("out_data/df_ccf.pickle", "rb") as fle:
#     df_ccf = pickle.load(fle)

# + {"jupyter": {"source_hidden": true}}
# # row_name_i
# # row_name_j
# # row_name_j_tmp = 'zuzwxhvuxe_tmp'

# index_in_cols = row_name_j in df_d_comp_prev.columns
# index_in_rows = row_name_j in df_d_comp_prev.index
# if index_in_cols and index_in_rows:
#     try:
#         d_ij = df_d_comp_prev.loc[row_name_i, row_name_j]
#     except:
#         tmp = 42
#         print("Couldn't find entry in previous D_ij matrix")

# + {"jupyter": {"source_hidden": true}}
# df_d_comp_prev.shape

# index_in_cols
# index_in_rows

# # df_d_comp_prev.loc['95c29e9f6h']

# '95c29e9f6h' in df_d_comp_prev.index
