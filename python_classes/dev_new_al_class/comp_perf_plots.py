# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python [conda env:research-new]
#     language: python
#     name: conda-env-research-new-py
# ---

# + jupyter={}
import os
import sys

import pickle

import chart_studio.plotly as py
import plotly.graph_objs as go

# #############################################################################
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "python_classes/active_learning"))
from active_learning import (
    ALBulkOpt,
    # ALGeneration,
    RegressionModel,
    FingerPrints,
    CandidateSpace,
    )

sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "python_classes"))
from ccf_similarity.ccf import CCF
# -

# # Script Inputs

# +
perc_of_structs = 50

shared_scatter_props = dict(
    mode="markers+lines",
    marker=dict(
        symbol="circle",
        size=8,
        line=dict(
            color='black',
            width=2
            )
        ),

    )


data_path_root = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow",
    "dev_new_al_class/out_data")

# +
from al_analysis import ALPerformance

data = []

# +
path_i = os.path.join(data_path_root, "ub-ucb_acq__bin_5.pickle")
with open(path_i, "rb") as fle:
    AL = pickle.load(fle)

# #############################################################################
ALPerf = ALPerformance(
    ALBulkOpt=AL,
    verbose=False)
ALPerf.num_sys_discovered(
    perc_of_structs=perc_of_structs,
    account_duplicates=True)

# #############################################################################
df = ALPerf.num_sys_discovered_df
trace = go.Scatter(
    x=df["num_dft"],
    y=df["num_ids_discovered"],
    marker=dict(color="pink"),
    )
trace.update(**shared_scatter_props)
data.append(trace)

# +
path_i = os.path.join(data_path_root, "random_acq__bin_5.pickle")
with open(path_i, "rb") as fle:
    AL = pickle.load(fle)

# #############################################################################
ALPerf = ALPerformance(
    ALBulkOpt=AL,
    verbose=False)
ALPerf.num_sys_discovered(
    perc_of_structs=perc_of_structs,
    account_duplicates=True)

# #############################################################################
df = ALPerf.num_sys_discovered_df
trace = go.Scatter(
    x=df["num_dft"],
    y=df["num_ids_discovered"],
    name="random w/ duplicates",
    marker=dict(color="black"),
    )
trace.update(**shared_scatter_props)
data.append(trace)

# +
path_i = os.path.join(data_path_root, "random_acq__bin_5.pickle")
with open(path_i, "rb") as fle:
    AL = pickle.load(fle)

# #############################################################################
ALPerf = ALPerformance(
    ALBulkOpt=AL,
    verbose=False)
ALPerf.num_sys_discovered(
    perc_of_structs=perc_of_structs,
    account_duplicates=False)

# #############################################################################
df = ALPerf.num_sys_discovered_df
trace = go.Scatter(
    x=df["num_dft"],
    y=df["num_ids_discovered"],
    name="random w/o duplicates",
    marker=dict(color="grey"),
    )
trace.update(**shared_scatter_props)
data.append(trace)
# -

# # Create figure

fig = go.Figure(data=data)
fig.show()
