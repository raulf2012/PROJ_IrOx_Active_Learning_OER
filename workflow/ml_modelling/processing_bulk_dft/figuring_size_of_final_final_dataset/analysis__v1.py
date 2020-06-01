# ---
# jupyter:
#   jupytext:
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

# +
# | - Import Modules
import os
print(os.getcwd())
import sys

import pickle

import numpy as np
import pandas as pd
#__|

# +
# /mnt/f/Dropbox/01_norskov/00_git_repos/PROJ_IrOx_Active_Learning_OER/workflow/ml_modelling

sys.path.insert(0, "/mnt/f/Dropbox/01_norskov/00_git_repos/PROJ_IrOx_Active_Learning_OER/workflow/ml_modelling")

from ml_methods import get_ml_dataframes, get_data_for_al

# +
# # %%capture

al_data_dict = get_data_for_al(
    stoich='AB3',
    verbose=False,
    drop_too_many_atoms=True,
    )

# +
al_data_dict.keys()
df_bulk_dft_0 = al_data_dict["df_bulk_dft"]

df_bulk_dft_0.shape

# +
# al_data_dict.keys()

df_static_irox = al_data_dict["df_static_irox"]

df_static_irox
# -

assert False

stoich = "AB2"
verbose = True
drop_too_many_atoms = True
# drop_too_many_atoms = False

# +
# def get_data_for_al(
# stoich="AB2",
# verbose=True,
# drop_too_many_atoms=True,
# ):
"""
"""

# | - get_data_for_al


sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import (ids_to_discard__too_many_atoms_path)

# | - Get all necessary dfs
df_dict = get_ml_dataframes(
    names=[
        "bulk_dft_data_path",
        "unique_ids_path",
        # "prototypes_data_path",
        "static_irox_structures_path",
        # "static_irox_structures_kirsten_path",
        # "oqmd_irox_data_path",
        "df_features_pre_opt_path",
        "df_features_pre_opt_kirsten_path",
        "df_features_post_opt_path",
        # "oer_bulk_structures_path",
        # "df_ccf_path",
        "df_dij_path",
        # "ids_to_discard__too_many_atoms_path",
        ],
    )

df_ids = df_dict.get("unique_ids", None)
df_bulk_dft = df_dict.get("bulk_dft_data", None)
df_bulk_dft = df_bulk_dft[df_bulk_dft.source == "raul"]

df_features_pre = df_dict.get("df_features_pre_opt", None)
# df_features_pre = df_dict.get("df_features_pre_opt_kirsten", None)
df_features_post = df_dict.get("df_features_post_opt", None)

df_dij = df_dict.get("df_dij", None)


df_static_irox = df_dict.get("static_irox_structures", None)
#__|

# +
# | - Filter ids to user specifications
df_ids = df_ids[
    (df_ids["stoich"] == stoich) & \
    (df_ids["source"] != "oqmd") & \
    (df_ids["source"] != "raul_oer") & \
    [True for i in range(len(df_ids))]]
ids = df_ids["unique_ids"].values

#__|

# #########################################################################

# | - DFT dataframe
df_i = df_bulk_dft

# print("isidfjisdjifjsidjf8yu2894h90832uy4908tyu98023wht0982quj098gtfujw3e")
# print(df_i.index.shape)
# print(df_i.index.unique().shape)

# Common ids between user ids and df
common_ids = list(set(df_i.index) & set(ids))

ids_not_in__df_i = [i for i in ids if i not in common_ids]

df_i = df_i.loc[common_ids]

if verbose:
    print("len(ids):", len(ids))
    print("len(common_ids)", len(common_ids))
    print("len(ids_not_in__bulk_dft_data):", len(ids_not_in__df_i))
    print("\n", "df_i.shape: ", df_i.shape, sep="")

df_bulk_dft = df_i
#__|
# -

with open(ids_to_discard__too_many_atoms_path, "rb") as fle:
    ids_to_drop__too_many_atoms = pickle.load(fle)



# +
ids_too_big = np.intersect1d(
    ids,
    ids_to_drop__too_many_atoms,
    )
len(ids_too_big)

# tmp

# +
ids = [i for i in ids if i not in ids_too_big]

len(ids)

# df_bulk_dft = df_bulk_dft.loc[ids]
# df_bulk_dft.shape

# +
tmp = np.intersect1d(
    df_bulk_dft.index.values,
    ids,
    )

len(tmp)

# +
np.setdiff1d(
    tmp,
    df_bulk_dft_0.index.values,
    )


# tmp0 = np.intersect1d(
#     tmp, 
#     df_bulk_dft_0.index.values,
#     )
# len(tmp0)

# +
# df_bulk_dft.loc["6dzhcimdxs"]

# +
# tmp

# + active=""
#
#
#
#
#
#
# -

for i in df_bulk_dft.index.unique().tolist():
    if i not in tmp:
        print(i)

np.setdiff1d(df_bulk_dft.index.unique().tolist(), tmp)

"xfbhx5nyvd" in tmp

# +
# atoms = df_bulk_dft.loc["xfbhx5nyvd"]["atoms"]

# atoms.get_number_of_atoms()
# -

df_bulk_dft.shape

assert False

# +
# 10 systems for AB3 not computed
# 210 systems for AB2 not computed

# +
109

10 + 210 - 10 - 2
# -

df_bulk_dft.shape

# +
# | - Features post-DFT
df_i = df_features_post

# Common ids between user ids and df
common_ids = list(set(df_i.index) & set(ids))

ids_not_in__df_i = [i for i in ids if i not in common_ids]

df_i = df_i.loc[common_ids]

if verbose:
    print("len(ids):", len(ids))
    print("len(common_ids)", len(common_ids))
    print("len(ids_not_in__bulk_dft_data):", len(ids_not_in__df_i))
    print("\n", "df_i.shape: ", df_i.shape, sep="")

df_features_post = df_i

# Only use post-DFT features from my data set
df_features_post = \
    df_features_post[df_features_post["data"]["source"] == "raul"]

df_features_post = df_features_post["voronoi"]
#__|
# -

df_bulk_dft.shape

# +
# | - Dropping certain rows
all_ids = list(set(
    df_bulk_dft.index.tolist() + \
    df_features_pre.index.tolist() + \
    df_features_post.index.tolist() ))


ids_to_drop = []

# #########################################################################
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/processing_bulk_dft/static_prototypes_structures/out_data",
    "ids_to_discard__proto_dupl.pickle")
with open(path_i, "rb") as fle:
    ids_to_discard__proto_dupl = pickle.load(fle)
    ids_to_drop.extend(ids_to_discard__proto_dupl)
# #########################################################################

if drop_too_many_atoms:
    # #####################################################################
    with open(ids_to_discard__too_many_atoms_path, "rb") as fle:
        ids_to_drop__too_many_atoms = pickle.load(fle)
        ids_to_drop.extend(ids_to_drop__too_many_atoms)

        # ids_to_drop = ids_to_drop__too_many_atoms
        # ids_to_drop = [i for i in ids_to_drop if i in all_ids]
    # #####################################################################


ids_to_drop = [i for i in ids_to_drop if i in all_ids]
print("len(ids_to_drop)", len(ids_to_drop))

df_features_pre = df_features_pre.drop(
    labels=ids_to_drop, axis=0)

tmp = df_bulk_dft.index.intersection(
    df_features_pre.index
    ).unique()
print("dksljkfjsijfijsdijfi", len(df_bulk_dft.index) - len(tmp))

df_bulk_dft = df_bulk_dft.loc[
    df_bulk_dft.index.intersection(
        df_features_pre.index
        ).unique()
    ]

#__|

# +
# len(ids_to_drop__too_many_atoms)

ids_to_drop__too_many_atoms

# ids.values

ids_to_drop_abx = np.intersect1d(
    ids_to_drop__too_many_atoms,
    ids.values,
    )

len(ids_to_drop_abx)
# -

tmp = np.intersect1d(
    ids_not_in__df_i,
    # ids_to_drop_abx,
    ids_to_drop__too_many_atoms,
    )
len(tmp)

697 - 131

df_bulk_dft.shape

458 + 29 

# +
# 458 structures for AB2
# 243 structures for AB2
# -

458 + 243

# +





# + jupyter={}
# out_dict = dict()

# out_dict["df_features_post"] = df_features_post
# out_dict["df_features_pre"] = df_features_pre
# out_dict["df_bulk_dft"] = df_bulk_dft

# # TEMP
# out_dict["df_ids"] = df_ids
# out_dict["df_dij"] = df_dij
# out_dict["df_static_irox"] = df_static_irox


# # return(out_dict)
# #__|

# + jupyter={}
# # | - Featurs pre-DFT
# df_i = df_features_pre

# # Common ids between user ids and df
# common_ids = list(set(df_i.index) & set(ids))

# ids_not_in__df_i = [i for i in ids if i not in common_ids]

# df_i = df_i.loc[common_ids]

# if verbose:
#     print("len(ids):", len(ids))
#     print("len(common_ids)", len(common_ids))
#     print("len(ids_not_in__bulk_dft_data):", len(ids_not_in__df_i))
#     print("\n", "df_i.shape: ", df_i.shape, sep="")

# df_features_pre = df_i

# df_features_pre = df_features_pre["voronoi"]
# #__|
