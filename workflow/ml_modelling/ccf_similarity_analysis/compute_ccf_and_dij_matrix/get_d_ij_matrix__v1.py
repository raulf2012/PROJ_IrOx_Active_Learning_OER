# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
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

# # Computing Structural Similarity of All Structures
# ---
#
# All pair-wise similarities are computed d_ij

read_from_PROJ_DATA = False
read_from_PROJ_DATA = True

# + [markdown] {"Collapsed": "false"}
# # Import Modules

# + {"Collapsed": "false", "jupyter": {}}
import os
print(os.getcwd())
import sys

import pickle

import numpy as np
import pandas as pd

# + {"Collapsed": "false"}
import time
t0 = time.time()

# + [markdown] {"Collapsed": "false"}
# # Script Inputs

# + {"Collapsed": "false"}
r_cut_off = 10.

r_vector = np.arange(1, 10, 0.02)

verbose = True

# + [markdown] {"Collapsed": "false"}
# ## Read Data

# + {"Collapsed": "false"}
with open("out_data/df_ccf.pickle", "rb") as fle:
    df_ccf = pickle.load(fle)

# #########################################################
if read_from_PROJ_DATA:
    path_i = os.path.join(
        os.environ["PROJ_DATA"], "04_IrOx_surfaces_OER",
        "PROJECT_COMPUTED_OUT_DATA/PROJ_IrOx_Active_Learning_OER",
        "workflow/ml_modelling/ccf_similarity_analysis/compute_ccf_and_dij_matrix",
        "out_data/df_d_ij_all.pickle")
    with open(path_i, "rb") as fle:
        df_d_comp_prev = pickle.load(fle)

else:
    try:
        with open("out_data/df_d_ij_all.pickle", "rb") as fle:
            df_d_comp_prev = pickle.load(fle)
    except:
        df_d_comp_prev = pd.DataFrame()


# #########################################################
print("df_ccf.shape:", df_ccf.shape)
print("df_d_comp_prev.shape:", df_d_comp_prev.shape)

# + [markdown] {"Collapsed": "false"}
# # Process data (Create D_ij matrix)

# + {"Collapsed": "false"}
df = df_ccf

result = np.zeros((len(df), len(df)))
for i_cnt, (row_name_i, row_i) in enumerate(df.iterrows()):
    if verbose:
        print(str(i_cnt).zfill(4), 75 * "#")
    
    row_in_prev_data = False
    for j_cnt, (row_name_j, row_j) in enumerate(df.iterrows()):
        if i_cnt == j_cnt:
            continue

        index_in_cols = row_name_j in df_d_comp_prev.columns
        index_in_rows = row_name_i in df_d_comp_prev.index
        if index_in_cols and index_in_rows:
            d_ij = df_d_comp_prev.loc[row_name_i, row_name_j]
            row_in_prev_data = True

            if verbose:
                tmp = 42
                # print("Parsed d_ij from previous data")
        else:
            if verbose: print("Computing d_ij from scratch")
            d_ij = cal_ccf_d(row_i["ccf"], row_j["ccf"])

        # d_ij = 1.11
        result[i_cnt][j_cnt] = d_ij

    if row_in_prev_data:
        if verbose: print("This row is the previous data (GOOD)")

# + {"Collapsed": "false"}
# #########################################################
df_d_comp = pd.DataFrame(
    result,
    index=df.index,
    columns=df.index)


# Pickling data ######################################################
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)
with open(os.path.join(directory, "df_d_ij_all.pickle"), "wb") as fle:
    pickle.dump(df_d_comp, fle)
# #####################################################################

# + {"Collapsed": "false"}
# path_i = os.path.join("out_data", "df_d_ij_all_temp.pickle")
path_i = os.path.join("out_data", "df_d_ij_all.pickle")
with open(path_i, "rb") as fle:
    df_d_comp = pickle.load(fle)

# + {"Collapsed": "false"}
time_elapsed = time.time() - t0
print(
    "Notebook took ",
    round(time_elapsed / 60, 4),
    "minutes to finish")
# -

print(20 * "# # ")
print("All done!")
assert False

# + {"Collapsed": "false", "active": ""}
#
#
#
#
#
#
#
#
#

# + {"Collapsed": "false", "jupyter": {"source_hidden": true}}
# import copy
# import pandas as pd

# from ase.visualize import view

# # Plotly
# import chart_studio.plotly as py
# import plotly.graph_objects as go

# # from StructurePrototypeAnalysisPackage.ccf import (
# from spap.ccf import (
#     struc2ccf,
#     cal_ccf_d,
#     cal_inter_atomic_d,
#     d2ccf,
#     weight_f,
#     pearson_cc,
#     gaussian_f,
#     element_tag,
#     cell_range,
#     count_atoms_dict,
#     )

# sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
# from proj_data_irox import (
#     static_irox_structures_path,
#     bulk_dft_data_path,
#     )
