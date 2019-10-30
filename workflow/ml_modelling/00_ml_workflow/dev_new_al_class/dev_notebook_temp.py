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
    ALGeneration,
    RegressionModel,
    FingerPrints,
    CandidateSpace,
    )

from al_analysis import ALAnalysis, ALAnimation
# -

# # Script Inputs

# +
stoich_i = "AB2"
verbose = False
num_gen_stop = 3

gp_settings = {
    "noise": 0.02542,
    "sigma_l": 1.0049,
    "sigma_f": 5.19,
    "alpha": 0.018,
    }
# -

# # Read Data

# +
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling"))
from ml_methods import get_data_for_al

out_dict = get_data_for_al(
    stoich="AB2", verbose=False,
    drop_too_many_atoms=True,
#     drop_too_many_atoms=False,
    )

df_bulk_dft = out_dict["df_bulk_dft"]
df_bulk_dft = df_bulk_dft[df_bulk_dft["source"] == "raul"]
df_bulk_dft = df_bulk_dft[["atoms", "energy_pa"]]
df_bulk_dft.columns.values[1] = "y_real"

df_features_pre = out_dict["df_features_pre"]
df_features_post = out_dict["df_features_post"]

df_ids = out_dict["df_ids"]


df_static_irox = out_dict["df_static_irox"]
df_dij = out_dict["df_dij"]

# +
ids_w_dft = df_bulk_dft.index

# TEMP | Reduce size of candidate space
np.random.seed(8)
# ids_w_dft = np.sort(np.random.choice(np.sort(ids_w_dft), size=30))
ids_w_dft = list(set(ids_w_dft))

df_bulk_dft = df_bulk_dft.loc[ids_w_dft]

df_features_pre = df_features_pre.loc[ids_w_dft]
df_features_post = df_features_post.loc[ids_w_dft]

# +
# ['6w6sbkvy6j' 'ckmlcq7hne' 'xfnavfzrmh' '7wva6g9284' 'm29j648g6i' 'xe7p8tc5z5']
# ['6w6sbkvy6j' '7wva6g9284' 'ckmlcq7hne' 'm29j648g6i' 'xe7p8tc5z5' 'xfnavfzrmh']

# + {"active": ""}
#
#
#
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

# +
RM = RegressionModel(
    opt_hyperparameters=True,
    gp_settings_dict=gp_settings,
    verbose=verbose,
    )

FP = FingerPrints(
    df_features_pre,
    df_features_post=df_features_post,
    pca_mode="num_comp",  # 'num_comp' or 'perc'
    pca_comp=11,
    pca_perc=None,
    verbose=verbose,
    )

CS = CandidateSpace(
    Y_data=df_bulk_dft,
    Y_key="y_real",
    FingerPrints=FP,
    )


name_i = "AL_" + stoich_i + "_" + str(num_gen_stop).zfill(2)
print("name:", name_i)
AL = ALBulkOpt(
    CandidateSpace=CS,
    RegressionModel=RM,
    DuplicateFinder=CCF,  # Optional
    num_seed_calcs=11,
    acquisition_bin=1,
    stop_mode="num_generations",
#     stop_mode=None,
    stop_num_generations=num_gen_stop,
    name="TEST_small",
    verbose=verbose,
    )

run_al = False
if run_al:
    AL.run_AL()

    AL.__save_state__()

# +
# assert False

# + {"active": ""}
#
#
#

# +
# #############################################################################
import pickle; import os
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow",
    "dev_new_al_class/out_data",
#     "TEST.pickle")
    "TEST_small.pickle")

with open(path_i, "rb") as fle:
    AL = pickle.load(fle)

ALAnim = ALAnimation(ALBulkOpt=AL, verbose=True)
# -

ALAnim.create_animation_parall(
    duration_long=1000 * 0.5,
    duration_short=800 * 0.5,
    )

# +
ALBulkOpt = ALAnim.ALBulkOpt
get_trace_j = ALAnim.get_trace_j


from pathos.multiprocessing import ProcessingPool
from multiprocessing import Pool, freeze_support
from functools import partial

inputs = [i.model for i in ALBulkOpt.al_gen_dict.values()]

results = Pool().map(partial(
    get_trace_j,  # METHOD

    # KWARGS
    prediction_key="y",
    uncertainty_key="err",
    plot_dft_instead_of_pred=True,
    trace_all_dft=True,
    trace_horiz_lines=True,
    marker_size=8,

    ), inputs)



# results
# -

assert False

# # TESTING | TEMP

# ## Testing ALBulkOpt

# + {"jupyter": {"source_hidden": true}}
self = AL

# #############################################################################
CandidateSpace = self.CandidateSpace
acquisition_bin = self.acquisition_bin
al_gen = self.al_gen
al_gen_dict = self.al_gen_dict
completed_ids = self.completed_ids
get_seed_ids = self.get_seed_ids
mode = self.mode
num_seed_calcs = self.num_seed_calcs
run_AL = self.run_AL
seed_ids = self.seed_ids
verbose = self.verbose
# #############################################################################
# -

# ## Testing ALGeneration

# +
AL_i = AL.al_gen_dict[3]
self = AL_i

# #############################################################################
completed_ids = self.completed_ids
CandidateSpace = self.CandidateSpace
model = self.model
verbose = self.verbose
# df_train = self.df_train
# df_test = self.df_test
verbose = self.verbose
acquisition_bin = self.acquisition_bin
RegressionModel = self.RegressionModel
# #############################################################################
# -

# ## Testing CandidateSpace

# + {"active": ""}
#
#
#
