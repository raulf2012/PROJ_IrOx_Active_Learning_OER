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

# # Import Modules

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
import plotly.graph_objs as go


# #############################################################################
from ccf import(
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
    bulk_dft_data_path,
    unique_ids_path,
    prototypes_data_path,
    static_irox_structures_path,
    oqmd_irox_data_path,
    voronoi_features_data_path,
    )
# -

# # Script Inputs

# +
r_cut_off = 10.

r_vector = np.arange(1, 10, 0.02)
# -

# # Read Data

# +
with open(static_irox_structures_path, "rb") as fle:
    df_structures = pickle.load(fle)

atoms_0 = df_structures.iloc[0]["atoms"]

# atoms_1 = df_structures.iloc[50]["atoms"]

atoms_1 = copy.deepcopy(atoms_0)
atoms_1.rattle(
    stdev=0.2,
#     stdev=0.01,
    )
# -

assert False


def get_ccf_df(ccf, r_cut_off, r_vector):
    """
    """
    df_i = pd.DataFrame()
    for key, value in ccf.items():
        df_i[key] = value

    df_i = df_i.set_index(r_vector)
    
    return(df_i)


# +
ccf_0 = struc2ccf(atoms_0, r_cut_off, r_vector)
ccf_1 = struc2ccf(atoms_1, r_cut_off, r_vector)

df_0 = get_ccf_df(ccf_0, r_cut_off, r_vector)
df_1 = get_ccf_df(ccf_1, r_cut_off, r_vector)
# -

cal_ccf_d(ccf_0, ccf_1)

# +
data = []
df_i = df_0
for col_i in list(df_i):
    trace_i = go.Scatter(
        x=list(df_i.index),
        y=df_i[col_i],
        mode="lines+markers",

        line=dict(
            color='red',
            width=0.8,
            ),

        )
    data.append(trace_i)

df_i = df_1
for col_i in list(df_i):
    trace_i = go.Scatter(
        x=list(df_i.index),
        y=df_i[col_i],
        mode="lines+markers",

        line=dict(
            color='black',
            width=0.8,
            ),

        )
    data.append(trace_i)

fig = go.Figure(data=data)
fig.show()

# + {"active": ""}
#
#
#

# + {"jupyter": {"source_hidden": true}}
# ccf_0 = struc2ccf(atoms_0, r_cut_off, r_vector)
# ccf_50 = struc2ccf(atoms_50, r_cut_off, r_vector)


# cal_ccf_d(ccf_0, ccf_50)

# # i_a_d =cal_inter_atomic_d(atoms_i, 5.)
# # d2ccf(i_a_d, r_cut_off, r_vector)

# ccf_0_list = []
# ccf_1_list = []

# r_vector_array = np.arange(1, 10, 2.)
# for r_vector in r_vector_array:
#     print(r_vector)

#     ccf_0 = struc2ccf(atoms_0, r_cut_off, r_vector)
#     ccf_1 = struc2ccf(atoms_1, r_cut_off, r_vector)

#     ccf_0_list.append(ccf_0)
#     ccf_1_list.append(ccf_1)

# df_0 = pd.DataFrame(
#     ccf_0_list
#     )

# ccf_0 = struc2ccf(atoms_0, r_cut_off, r_vector)
# ccf_1 = struc2ccf(atoms_1, r_cut_off, r_vector)
