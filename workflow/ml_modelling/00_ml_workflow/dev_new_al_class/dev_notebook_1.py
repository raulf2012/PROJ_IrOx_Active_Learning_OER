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
import os
import sys

import pickle

import numpy as np
import pandas as pd

# + {"jupyter": {"source_hidden": true}}
import numpy as np
import pandas as pd

import gpflow

from sklearn.decomposition import PCA


from catlearn.regression.gaussian_process import GaussianProcess
from catlearn.preprocess.clean_data import (
    clean_infinite,
    clean_variance,
    clean_skewness)
from catlearn.preprocess.scaling import standardize

# + {"jupyter": {"source_hidden": true}}
sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import (ids_to_discard__too_many_atoms_path)
# -

# # Script Inputs

# +
verbose = False

stoich_i = "AB2"

num_gen_stop = 4
# -

# # Read Data

# +
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling"))
from ml_methods import get_data_for_al

out_dict = get_data_for_al(stoich="AB2", verbose=False)

df_bulk_dft = out_dict["df_bulk_dft"]
df_bulk_dft = df_bulk_dft[df_bulk_dft["source"] == "raul"]
df_bulk_dft = df_bulk_dft[["atoms", "energy_pa"]]
df_bulk_dft.columns.values[1] = "y_real"

df_features_pre = out_dict["df_features_pre"]
df_features_post = out_dict["df_features_post"]

# +
# #############################################################################
with open(ids_to_discard__too_many_atoms_path, "rb") as fle:
    ids_to_drop__too_many_atoms = pickle.load(fle)
    ids_to_drop__too_many_atoms = \
        [i for i in ids_to_drop__too_many_atoms if i in df_features_pre.index]

# #############################################################################
df_features_pre = df_features_pre.drop(
    labels=ids_to_drop__too_many_atoms,
    axis=0)

# #############################################################################
df_features_post = df_features_post.loc[
    [i for i in df_features_post.index if i in df_features_pre.index]]

# #############################################################################
df_bulk_dft = df_bulk_dft.loc[
    [i for i in df_bulk_dft.index if i in df_features_pre.index]]
# -

sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "python_classes/active_learning"))
from active_learning import (
    ALBulkOpt,
    ALGeneration,
    RegressionModel,
    FingerPrints,
    CandidateSpace)

# + {"jupyter": {"source_hidden": true}}
gp_settings = {
    "noise": 0.02542,
    "sigma_l": 1.0049,
    "sigma_f": 5.19,
    "alpha": 0.018,
    }

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

AL = ALBulkOpt(
    CandidateSpace=CS,
    RegressionModel=RM,
    num_seed_calcs=11,
    acquisition_bin=40,
    stop_mode="num_generations",
    stop_num_generations=num_gen_stop,
    name=name_i,
    verbose=verbose,
    )

# +
AL.run_AL()

AL.__save_state__()
# -

# Pickling data ###############################################################
import os; import pickle
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)
with open(os.path.join(directory, "temp_" + name_i), "wb") as fle:
    pickle.dump(AL, fle)
# #############################################################################

temp_data = AL.al_gen_dict[0]
with open(os.path.join(directory, "temp_single_gen.pickle"), "wb") as fle:
    pickle.dump(temp_data, fle)

# +
with open(os.path.join(directory, "df_bulk_dft.pickle"), "wb") as fle:
    pickle.dump(df_bulk_dft, fle)

with open(os.path.join(directory, "df_features_post.pickle"), "wb") as fle:
    pickle.dump(df_features_post, fle)

with open(os.path.join(directory, "df_features_pre.pickle"), "wb") as fle:
    pickle.dump(df_features_pre, fle)

# + {"active": ""}
#
#
#
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

# + {"jupyter": {"source_hidden": true}}
AL_i = AL.al_gen_dict[3]
self = AL_i

# #############################################################################
completed_ids = self.completed_ids
CandidateSpace = self.CandidateSpace
verbose = self.verbose
df_train = self.df_train
df_test = self.df_test
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

# + {"jupyter": {"source_hidden": true}}
# acquisition_method="gp_ucb"  # 'gp_ucb' or 'random'
# # #####################################################################
# acquisition_bin = self.acquisition_bin
# model = self.model
# # #####################################################################
# if acquisition_method == "gp_ucb":
#     acquisition_ids_ordered = self.acquisition_gp_ucb(model, kappa=1.)
# elif acquisition_method == "random":
#     acquisition_ids_ordered = self.acquisition_random(model)


# # Model ordered based on acquisition function
# model_tmp = model.loc[acquisition_ids_ordered]

# # Remove rows for which DFT data is not available
# model_data_avail = model_tmp[~model_tmp["y_real"].isna()]

# # Remove rows which have already been acquired
# # Only acquire what hasn't already been acquired
# model_data_avail = model_data_avail[model_data_avail["acquired"] == False]

# new_acquis_ids = model_data_avail.index[0:acquisition_bin].tolist()

# +
# # Pickling data ###############################################################
# import os; import pickle
# directory = "out_data"
# if not os.path.exists(directory): os.makedirs(directory)
# with open(os.path.join(directory, "AL_ab2.pickle"), "wb") as fle:
#     pickle.dump(AL, fle)
# # #############################################################################
