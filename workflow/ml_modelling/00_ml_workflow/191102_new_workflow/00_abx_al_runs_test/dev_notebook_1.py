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

# + {"jupyter": {"source_hidden": true}}
# %%capture

import os
import sys

import copy

import time
import pickle

import numpy as np
import pandas as pd

import gpflow

from sklearn.decomposition import PCA

# #############################################################################
from catlearn.regression.gaussian_process import GaussianProcess
from catlearn.preprocess.clean_data import (
    clean_infinite,
    clean_variance,
    clean_skewness)
from catlearn.preprocess.scaling import standardize

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

from al_analysis import ALAnalysis, ALAnimation

# #############################################################################
from IPython.display import display
# -

# # Script Inputs

from inputs import (
    stoich_i,
    verbose,
    gp_settings,
    name_i,
    )

# # Read Data

# + {"jupyter": {"source_hidden": true}}
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling"))
from ml_methods import get_data_for_al

out_dict = get_data_for_al(
    stoich=stoich_i, verbose=False,
    drop_too_many_atoms=True,
#     drop_too_many_atoms=False,
    )

df_bulk_dft = out_dict["df_bulk_dft"]
df_bulk_dft = df_bulk_dft[df_bulk_dft["source"] == "raul"]

# df_bulk_dft = df_bulk_dft[["atoms", "energy_pa"]]
df_bulk_dft = df_bulk_dft[["atoms", "dH"]]
df_bulk_dft.columns.values[1] = "y_real"

df_features_pre = out_dict["df_features_pre"]
df_features_post = out_dict["df_features_post"]

df_ids = out_dict["df_ids"]


df_static_irox = out_dict["df_static_irox"]
df_dij = out_dict["df_dij"]
# -

# # Filter to candidates w/ DFT energy

# + {"jupyter": {"source_hidden": true}}
ids_w_dft = df_bulk_dft.index

# TEMP | Reduce size of candidate space
# np.random.seed(8)
# ids_w_dft = np.sort(np.random.choice(np.sort(ids_w_dft), size=200))
ids_w_dft = list(set(ids_w_dft))

df_bulk_dft = df_bulk_dft.loc[ids_w_dft]

df_features_pre = df_features_pre.loc[ids_w_dft]
df_features_post = df_features_post.loc[ids_w_dft]

# + {"active": ""}
#
#
#
#
# -

# # CCF Class

# + {"jupyter": {"source_hidden": true}}
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "python_classes"))
from ccf_similarity.ccf import CCF

d_thresh = 0.02
CCF = CCF(
    df_dij=df_dij,
    d_thresh=d_thresh)

# + {"jupyter": {"source_hidden": true}}
RM = RegressionModel(
    opt_hyperparameters=True,
    gp_settings_dict=gp_settings,
    verbose=verbose,
    )

FP = FingerPrints(
    df_features_pre,
    df_features_post=df_features_post,
    pca_mode="num_comp",  # 'num_comp' or 'perc'
    pca_comp=10,
    pca_perc=None,
    verbose=verbose,
    )

CS = CandidateSpace(
    Y_data=df_bulk_dft,
    Y_key="y_real",
    FingerPrints=FP,
    )

# +
AL = ALBulkOpt(
    CandidateSpace=CS,
    RegressionModel=RM,
    DuplicateFinder=CCF,  # Optional
    # num_seed_calcs=11,
    num_seed_calcs=5,
    acquisition_bin=10,
    stop_mode="num_generations",
#     stop_mode=None,
    stop_num_generations=3,
    name=name_i,
    verbose=verbose,
    acquisition_method="gp_ucb",
    )

run_al = True
if run_al:
    AL.run_AL()
    AL.duplicate_system_history_analysis()
    AL.__save_state__()
# -

assert False

# + {"active": ""}
#
#
#

# + {"jupyter": {"source_hidden": true}}
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow",
    "dev_new_al_class/out_data",

    "TEST__acq_5_all.pickle"
    # "TEST__acq_10.pickle"
    )

with open(path_i, "rb") as fle:
    AL = pickle.load(fle)

ALAnim = ALAnimation(
    ALBulkOpt=AL,
    marker_color_dict=id_color_dict,
    verbose=True)

if False:
    ALAnim.create_animation(
        # duration_long=1000 * 0.5,
        # duration_short=800 * 0.5,
        duration_long=1000 * 4,
        duration_short=800 * 4,
        serial_parallel="parallel",  # 'serial' or 'parallel'
        # marker_color_dict=id_color_dict,
        )

# + {"jupyter": {"source_hidden": true}}
from al_analysis import ALPerformance

ALPerf = ALPerformance(
    ALBulkOpt=AL,
    verbose=False,
    )
ALPerf.num_sys_discovered()

df = ALPerf.num_sys_discovered_df

# df

# + {"jupyter": {"source_hidden": true}}
import chart_studio.plotly as py
import plotly.graph_objs as go
import os


trace = go.Scatter(
    x=df["num_dft"],
    y=df["num_ids_discovered"],
    mode="markers",
    marker=dict(
        symbol="circle",
        color="grey",
        size=14,
        line=dict(
            color='black',
            width=2
            )
        ),
    )

data = [trace]

fig = go.Figure(data=data)
fig.show()

# + {"active": ""}
#
#
#
# -

# # SANDBOX

# + {"jupyter": {"source_hidden": true}}
# #############################################################################
self = AL

al_gen_dict = AL.al_gen_dict
duplicate_ids = AL.duplicate_ids

duplicate_swap_dict = AL.duplicate_swap_dict
# #############################################################################

# + {"jupyter": {"source_hidden": true}}
CandidateSpace = AL.CandidateSpace
FingerPrints = CandidateSpace.FingerPrints
df_pre = FingerPrints.df_pre
all_indices = df_pre.index.tolist()

# + {"jupyter": {"source_hidden": true}}
# #############################################################################
gen_i = 0

AL_i = al_gen_dict[gen_i]
self = AL_i

model = self.model
prev_duplicate_ids = self.prev_duplicate_ids
indices_that_are_duplicates = self.indices_that_are_duplicates
# #############################################################################
