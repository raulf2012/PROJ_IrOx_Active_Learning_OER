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

# %load_ext autoreload
# %autoreload 2

# + jupyter={}
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

sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "python_classes"))
from ccf_similarity.ccf import CCF

from pathos.multiprocessing import ProcessingPool

# +
from pathos.multiprocessing import ProcessingPool
class Bar:
    def foo(self, name):
        return(len(str(name)))
    def boo(self, things):
        for thing in things:
            self.sum += self.foo(thing)
        return(self.sum)
    sum = 0

b = Bar()
results = ProcessingPool().map(b.boo, [[12,3,456],[8,9,10],['a','b','cde']])

results
# -

assert False

# +
# #############################################################################
import pickle; import os
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow",
    "dev_new_al_class/out_data",
    "TEST_small.pickle")
    # "AL_AB2_05.pickle")
    # "AL_AB2_06.pickle")

with open(path_i, "rb") as fle:
    AL = pickle.load(fle)

# #############################################################################
self = AL

seed_ids = AL.seed_ids
index_acq_gen_dict = AL.index_acq_gen_dict
# #############################################################################
AL_i = AL.al_gen_dict[2]
self = AL_i

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


model[model["acquired"] == True]
# -

AL_i.__run_duplicate_analysis__()

self = AL_i

# +
# #####################################################################
acquisition_bin = self.acquisition_bin
model = self.model
DuplicateFinder = self.DuplicateFinder
index_acq_gen_dict = self.index_acq_gen_dict
# #####################################################################


# #####################################################################
# #####################################################################
# | - Apply 'gen_acquired' to model df
def method(row_i, index_acq_gen_dict):
    index_i = row_i.name
    gen_i = index_acq_gen_dict.get(index_i, np.nan)
    return(gen_i)

model["gen_acquired"] = model.apply(
    method, axis=1,
    args=(index_acq_gen_dict, ))
#__|


# #####################################################################
# #####################################################################
# #####################################################################
model_acq = model[model["acquired"] == True]

# Only consider duplicates in the set of structures that have been computed
filter_ids = model_acq.index.tolist()

simil_dict_master = dict()
for index_i in model_acq.index.tolist():
    simil_dict = DuplicateFinder.i_all_similar(
        index_i, filter_ids=filter_ids)

    simil_dict_master[index_i] = simil_dict

keys_to_delete = []
for key, val in simil_dict_master.items():
    if val == dict() or val is None:
        keys_to_delete.append(key)

for key in keys_to_delete:
    del simil_dict_master[key]




if len(simil_dict_master.keys()) == 0:
    self.indices_that_are_duplicates = []

else:
    # #####################################################################
    # #####################################################################
    # #####################################################################
    keys = list(simil_dict_master.keys())

    tmp_list = [np.array(list(i.keys())) for i in simil_dict_master.values()]
    all_ids_from_duplicate_analysis = keys + list(np.hstack(tmp_list))
    all_ids_from_duplicate_analysis = list(set(all_ids_from_duplicate_analysis))

    # #####################################################################
    # #####################################################################
    # #####################################################################

    # Tracks ids that have already been identified as duplicates
    # Don't consider further, already being removed/treated
    indices_that_are_duplicates = []
    for key, val in simil_dict_master.items():
        
        if key in indices_that_are_duplicates:
            continue

        ids_of_duplicates = [key] + list(val.keys())

        ids_of_duplicates = \
            [i for i in ids_of_duplicates if i not in indices_that_are_duplicates]

        # Skip loop if no duplicate ids are present
        if len(ids_of_duplicates) <= 1:
            continue

        df_tmp = model.loc[ids_of_duplicates].sort_values("gen_acquired")
        assert df_tmp.shape[0] > 1, "Only one row in df_tmp"
        self.TEMP__df_tmp = df_tmp

        earlist_gen = df_tmp.iloc[0]["gen_acquired"]

        # Check that there is only 1 row from previous generations
        # If this is working, then all duplicates are removed as they occur,
        # so there shouldn't be any left overs
        earliest_acq_row = df_tmp.iloc[0]
        generations_acquired = df_tmp["gen_acquired"].tolist()

        if len(list(set(generations_acquired))) == 1:
            print("All duplicates acquired at the same gen | OK")
        else:
            mess = "There shouldn't be more than one duplicate from previous generations"
            num_early_gens = generations_acquired.count(earliest_acq_row["gen_acquired"])
            # assert num_early_gens == 1, mess


        # Are there multiple early gen rows to choose from?
        # Should only happen if multiple are acquired at once
        multiple_early_gens_present = False
        if len(list(set(generations_acquired))) == 1:
            print("multiple_early_gens_present")
            print("TEMP")
            multiple_early_gens_present = True
            # break


        selected_row = \
            df_tmp[df_tmp["gen_acquired"] == earlist_gen].sort_values("y_real").iloc[0]

        # | - OLD | Trying to replace value for lowest energy duplicate
        # lowest_y_row = df_tmp.sort_values("y_real").iloc[0]
        # TEMP
        # lowest_y_row = df_tmp.sort_values("y_real").iloc[1]
        # if earliest_acq_row.name != lowest_y_row.name:
        #     print(earliest_acq_row.name, lowest_y_row.name)
        # #     model.loc[lowest_y_row.name]
        # #     model.rename(
        # #         index={
        # #             lowest_y_row.name: earliest_acq_row.name + "_TEMP",
        # #             earliest_acq_row.name: lowest_y_row.name,
        # #             }, inplace=True)
        # #     model.rename(
        # #         index={
        # #             earliest_acq_row.name + "_TEMP": earliest_acq_row.name,
        # #             }, inplace=True)
        #__|

        indices_that_are_duplicates_i = df_tmp.index.tolist()
        indices_that_are_duplicates_i.remove(selected_row.name)

        indices_that_are_duplicates.extend(indices_that_are_duplicates_i)


    indices_that_are_duplicates = list(set(indices_that_are_duplicates))
    self.indices_that_are_duplicates = indices_that_are_duplicates

    [i for i in all_ids_from_duplicate_analysis if i not in indices_that_are_duplicates]

# -

df_tmp

# +
ids_of_duplicates

# generations_acquired

ids_of_duplicates = [key] + list(val.keys())

# ids_of_duplicates = \
#     [i for i in ids_of_duplicates if i not in indices_that_are_duplicates]

# +
val.keys()

print(key)

print(val)
# -

indices_that_are_duplicates

assert False



# + active=""
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

# +
# for key, AL_i in AL.al_gen_dict.items():
#     duplicates = AL_i.indices_that_are_duplicates

# model["duplicates"] = [True if i in duplicates else False for i in model.index.tolist()]

# model

# +
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling"))
from ml_methods import get_data_for_al, get_ml_dataframes

out_dict = get_data_for_al(
    stoich="AB2",
    verbose=False,
    drop_too_many_atoms=True)

# df_features_pre = out_dict["df_features_pre"]
df_static_irox = out_dict["df_static_irox"]
# df_bulk_dft = out_dict["df_bulk_dft"]
df_dij = out_dict["df_dij"]

all_indices = df_static_irox.index


sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "python_classes"))
from ccf_similarity.ccf import CCF

d_thresh = 0.02
CCF = CCF(
    df_dij=df_dij,
    d_thresh=d_thresh)

# +
model_acq = model[model["acquired"] == True]

# Only consider duplicates in the set of structures that have been computed
filter_ids = model_acq.index.tolist()

simil_dict_master = dict()
for index_i in model_acq.index.tolist():
    simil_dict = CCF.i_all_similar(
        index_i, filter_ids=filter_ids)

    simil_dict_master[index_i] = simil_dict

keys_to_delete = []
for key, val in simil_dict_master.items():
    if val == dict() or val is None:
        keys_to_delete.append(key)

for key in keys_to_delete:
    del simil_dict_master[key]

# simil_dict_master

# +
keys = list(simil_dict_master.keys())

tmp_list = [np.array(list(i.keys())) for i in simil_dict_master.values()]
all_ids_from_duplicate_analysis = keys + list(np.hstack(tmp_list))

all_ids_from_duplicate_analysis = list(set(all_ids_from_duplicate_analysis))

all_ids_from_duplicate_analysis


# +
def method(row_i, index_acq_gen_dict):
    index_i = row_i.name
    gen_i = index_acq_gen_dict.get(index_i, np.nan)
    return(gen_i)

model["gen_acquired"] = model.apply(
    method, axis=1,
    args=(index_acq_gen_dict, ))

# model[~model["gen_acquired"].isna()]

# +
# model[model["acquired"] == True].loc[[
#     "xhbabrx4zq",
#     "zszinjv3zf",
#     "727lmkmq74",
#     "zgntxjxrvj",
#     ]]

# +
indices_that_are_duplicates = []
for key, val in simil_dict_master.items():
    ids_of_duplicates = [key] + list(val.keys())

    ids_of_duplicates = \
        [i for i in ids_of_duplicates if i not in indices_that_are_duplicates]

    # Skip loop if no duplicate ids are present
    if len(ids_of_duplicates) == 0:
        continue

    df_tmp = model.loc[ids_of_duplicates].sort_values("gen_acquired")
    earlist_gen = df_tmp.iloc[0]["gen_acquired"]

    # Check that there is only 1 row from previous generations
    # If this is working, then all duplicates are removed as they occur,
    # so there shouldn't be any left overs
    earliest_acq_row = df_tmp.iloc[0]
    generations_acquired = df_tmp["gen_acquired"].tolist()
    if len(list(set(generations_acquired))) == 1:
        print("All duplicates acquired at the same gen | OK")
    else:
        mess = "There shouldn't be more than one duplicate from previous generations"
        num_early_gens = generations_acquired.count(earliest_acq_row["gen_acquired"])
        # assert num_early_gens == 1, mess


    # Are there multiple early gen rows to choose from?
    # Should only happen if multiple are acquired at once
    multiple_early_gens_present = False
    if len(list(set(generations_acquired))) == 1:
        print("multiple_early_gens_present")
        multiple_early_gens_present = True

#         break

    selected_row = \
        df_tmp[df_tmp["gen_acquired"] == earlist_gen].sort_values("y_real").iloc[0]

    # lowest_y_row = df_tmp.sort_values("y_real").iloc[0]
    # TEMP
    # lowest_y_row = df_tmp.sort_values("y_real").iloc[1]
    # if earliest_acq_row.name != lowest_y_row.name:
    #     print(earliest_acq_row.name, lowest_y_row.name)
    # #     model.loc[lowest_y_row.name]
    # #     model.rename(
    # #         index={
    # #             lowest_y_row.name: earliest_acq_row.name + "_TEMP",
    # #             earliest_acq_row.name: lowest_y_row.name,
    # #             }, inplace=True)
    # #     model.rename(
    # #         index={
    # #             earliest_acq_row.name + "_TEMP": earliest_acq_row.name,
    # #             }, inplace=True)

    
    indices_that_are_duplicates_i = df_tmp.index.tolist()
    indices_that_are_duplicates_i.remove(selected_row.name)

    indices_that_are_duplicates.extend(indices_that_are_duplicates_i)
    
indices_that_are_duplicates = list(set(indices_that_are_duplicates))
# -

[i for i in all_ids_from_duplicate_analysis if i not in indices_that_are_duplicates]

assert False

# df_tmp.index.tolist()
# print(df_tmp.index.tolist())
import copy
ids_to_drop = copy.deepcopy(ids_of_duplicates)
ids_to_drop.remove(earliest_acq_row.name)

# +
print(lowest_y_row.name)

print(earliest_acq_row.name)

# +
ids_to_drop = model.index.intersection(ids_to_drop).tolist()
# ids_of_duplicates

model.drop(labels=ids_to_drop, axis=0, inplace=True)
# -

model.shape

# +
# lowest_y_row.name
# earliest_acq_row.name

# model.loc['xhbabrx4zq']
model.loc['zszinjv3zf']

# 'zszinjv3zf' in 
# model.index.tolist()

# # model.drop?

# + active=""
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
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

assert False

# +
# index_list = []
# new_acquisition_dict = {}
# for gen_i, AL_i in AL.al_gen_dict.items():
#     index_acq_gen_dict_i = dict()
#     for index_j in AL_i.new_acquisition:
#         index_list.append(index_j)
#         index_acq_gen_dict_i[index_j] = int(gen_i)

#     new_acquisition_dict.update(index_acq_gen_dict_i)

# mess = "Seems like an id was acquired in more than 1 generation?"
# assert len(index_list) == len(set(index_list)), mess

# +
import os
import sys

import pickle

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import (
    bulk_dft_data_path,
    ids_to_discard__too_many_atoms_path,
    unique_ids_path,
    df_dij_path)

# +
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling"))
from ml_methods import get_data_for_al, get_ml_dataframes

out_dict = get_data_for_al(
    stoich="AB2",
    verbose=False,
    drop_too_many_atoms=True)

print(out_dict.keys())

df_features_pre = out_dict["df_features_pre"]
df_static_irox = out_dict["df_static_irox"]
df_bulk_dft = out_dict["df_bulk_dft"]
df_dij = out_dict["df_dij"]

all_indices = df_static_irox.index

sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "python_classes"))
from ccf_similarity.ccf import CCF

d_thresh = 0.02
CCF = CCF(
    df_dij=df_dij,
    d_thresh=d_thresh,
    )

# +
# CCF.df_dij.loc["v5ckmsnqxi"]
# CCF.df_dij.loc["z3ngxhbrz4"]

# +
index_i = "mkbrzh8kv5"

filter_ids = all_indices

simil_dict_master = dict()
for index_i in all_indices:
    simil_dict = CCF.i_all_similar(index_i, filter_ids=filter_ids)
    simil_dict_master[index_i] = simil_dict

# + active=""
#
#
#
#
#
#
#
#

# + jupyter={}
# v5ckmsnqxi

# + jupyter={}
# CCF.df_dij.loc[]

# + jupyter={}
# tmp = [i for i in all_indices if i in df_dij.index.tolist()]
# len(tmp)
# len(all_indices)

# + jupyter={}
# self = CCF
# d_thresh = self.d_thresh
# df_dij = self.df_dij

# index_i = "mkbrzh8kv5"
# # index_j = "folatese_05"

# # filter_ids = all_indices
# filter_ids = None
# # #####################################################################


# if filter_ids is not None:
#     filter_ids_inter = df_dij.index.intersection(filter_ids)
#     df_dij = df_dij.loc[filter_ids_inter, filter_ids_inter]

# row_i = df_dij.loc[index_i]
# row_i = row_i.drop(labels=index_i)

# out_dict = row_i[row_i < d_thresh].to_dict()

# out_dict

# + jupyter={}
# self = CCF
# d_thresh = self.d_thresh
# df_dij = self.df_dij


# index_i = "budabebu_36"
# index_j = "folatese_05"

# dij = df_dij.loc[index_i, index_j]

# similar = False
# if dij < d_thresh:
#     similar = True
# elif dij > d_thresh:
#     similar = False
# else:
#     assert False, "AHHHHHHHH!!!"

# # return(similar)

# CCF.i_j_similar()

# + active=""
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

assert False

# + jupyter={}
df_static_irox["static_id"]

# + jupyter={}
df_dij = out_dict["df_dij"]


ids_static = df_dij.index.intersection(df_static_irox["static_id"])
ids_completed_post_dft = df_dij.index.intersection(df_features_pre.index)

ids_dij = ids_static.tolist() + ids_completed_post_dft.tolist()

df_dij = df_dij.loc[ids_dij, ids_dij]

# + jupyter={}
# "bsv4nex29l" in df_dij.index

# for index in df_features_pre.index.tolist():
#     if index not in df_dij.index.tolist():
#         print(index)

# + jupyter={}
# df_static_irox.head()
# df_static_irox = df_static_irox.loc[
#     df_static_irox.index.intersection(
#         df_features_pre.index
#         ).unique()
#     ]

# + active=""
#
#
#
#
#
#
#
#
#

# + jupyter={}
assert False

# + jupyter={}
out_dict = get_ml_dataframes(
    names=[
        # "bulk_dft_data_path",
        # "unique_ids_path",
        # "prototypes_data_path",
        "static_irox_structures_path",
        # "static_irox_structures_kirsten_path",
        # "oqmd_irox_data_path",
        # "df_features_pre_opt_path",
        # "df_features_pre_opt_kirsten_path",
        # "df_features_post_opt_path",
        # "df_features_path",
        # "df_features_cleaned_path",
        # "df_features_cleaned_pca_path",
        # "oer_bulk_structures_path",
        # "df_ccf_path",
        "df_dij_path",
        # "ids_to_discard__too_many_atoms_path",
        ],

    )

out_dict.keys()


df_static_irox = out_dict["static_irox_structures"]
df_static_irox[
    (df_static_irox["stoich"] == stoich) & \
    (df_static_irox["source"] == "chris")
    ]

# + jupyter={}
df_dij_path_tmp = df_dij_path[0:-18] + "df_d_ij_all_temp.pickle"
with open(df_dij_path_tmp, "rb") as fle:
    df_dij_dft = pickle.load(fle)
    print("df_dij_dft.shape:", df_dij_dft.shape)

# + jupyter={}
rows_equal_cols = all(df_dij_dft.index == df_dij_dft.columns)



# + jupyter={}
assert False
# -

# # Figuring out the issue with the GP (giving error)

# + jupyter={}
import pandas as pd

# + jupyter={}
# #############################################################################
import pickle; import os
path_i = os.path.join(
    os.environ["HOME"],
    "__temp__",
    "TEMP.pickle")
with open(path_i, "rb") as fle:
    data = pickle.load(fle)
# #############################################################################

# + jupyter={}
data.keys()

K = data["K"]
kernel_list = data["kernel_list"]
train_matrix = data["train_matrix"]
theta = data["theta"]

# + jupyter={}
kernel_list

# + jupyter={}
from scipy.linalg import cho_solve, cho_factor
L, lower = cho_factor(K, overwrite_a=False, lower=True, check_finite=True)

# + jupyter={}
pd.DataFrame(K)
