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

# # New ML Active Learning Workflow
# ---

# + active=""
# Following system is failing with voro fingerprinting:
# z39g648rnl
# -

# # Import Modules

# +
import os
import sys
import copy

import pickle

import pandas as pd

from sklearn.decomposition import PCA

from protosearch.ml_modelling.fingerprint import (
    FingerPrint,
    VoronoiFingerprint
    )

# #############################################################################
pd.set_option('display.max_rows', None)
os.getcwd()



# #############################################################################
sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import (
    bulk_dft_data_path,
    unique_ids_path,
    prototypes_data_path,
    static_irox_structures_path,
    static_irox_structures_kirsten_path,
    df_features_pre_opt_path,
    df_features_pre_opt_kirsten_path,
    df_features_post_opt_path,
    )
# -

directory = "out_data"
if not os.path.exists(directory):
    os.makedirs(directory)

# # Read Data

# +
with open(bulk_dft_data_path, "rb") as fle:
    df_bulk_dft = pickle.load(fle)

with open(static_irox_structures_path, "rb") as fle:
    df_struct = pickle.load(fle)

with open(static_irox_structures_kirsten_path, "rb") as fle:
    df_struct_kirsten = pickle.load(fle)

df_ids = pd.read_csv(unique_ids_path)
# -

"6dzhcimdxs" in df_bulk_dft.index

# +
# assert False
# -

# # Remove rows with missing atoms objects

# Removing missing data
df_bulk_dft = df_bulk_dft[df_bulk_dft["atoms"].notnull()]
df_struct = df_struct[df_struct["atoms"].notnull()]

df_struct.columns

# # Processing pre-opt Fingerprints

# +
# FP_struct = FingerPrint(**{
#     "feature_methods": ["voronoi"],
#     "input_data": df_struct,
# #     "input_data": df_combined,
#     "input_index": ["atoms"]})

# FP_struct.generate_fingerprints()
# df_features_pre_opt = FP_struct.fingerprints

# # #############################################################################
# with open(df_features_pre_opt_path, "wb") as fle:
#     pickle.dump(df_features_pre_opt, fle)

# +
# FP_struct = FingerPrint(**{
#     "feature_methods": ["voronoi"],
#     "input_data": df_struct_kirsten,
# #     "input_data": df_combined,
#     "input_index": ["atoms"]})

# FP_struct.generate_fingerprints()
# df_features_pre_opt = FP_struct.fingerprints

# # #############################################################################
# with open(df_features_pre_opt_kirsten_path, "wb") as fle:
#     pickle.dump(df_features_pre_opt, fle)

# +
with open(df_features_pre_opt_path, "rb") as fle:
    df_features_pre_opt = pickle.load(fle)

with open(df_features_pre_opt_kirsten_path, "rb") as fle:
    df_features_pre_opt_kirsten = pickle.load(fle)
# -

with open(df_features_post_opt_path, "rb") as fle:
    df_features_post_opt = pickle.load(fle)

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

# +
df_features_post_opt["data", "INDEX_OLD"] = df_features_post_opt.index
df_features_post_opt["data", "INDEX_NEW"] = df_features_post_opt["data", "INDEX_OLD"] + "_" + df_features_post_opt["data"]["source"]

df_features_post_opt = df_features_post_opt.set_index(df_features_post_opt["data", "INDEX_NEW"])
df_features_post_opt = df_features_post_opt.drop(labels=[["data", "INDEX_NEW"]], axis=1)

print("df_features_post_opt.shape:", df_features_post_opt.shape)
# -

df_bulk_dft["INDEX_NEW"] = df_bulk_dft.index + "_" + df_bulk_dft["source"]
df_bulk_dft["INDEX_OLD"] = df_bulk_dft.index
df_bulk_dft = df_bulk_dft.set_index("INDEX_NEW")

ids_to_process = [i for i in df_bulk_dft.index if i not in df_features_post_opt.index]
df_bulk_dft_not_processed = df_bulk_dft.loc[ids_to_process]

print("df_bulk_dft_not_processed.shape:", df_bulk_dft_not_processed.shape)
# df_bulk_dft_not_processed.head()

# +
# df_bulk_dft[df_bulk_dft["stoich"] == "AB2"].sort_values("energy_pa")

# +
FP_struct = FingerPrint(**{
    "feature_methods": ["voronoi"],
    "input_data": df_bulk_dft_not_processed,
    "input_index": ["atoms"]})

FP_struct.generate_fingerprints()
df_features_post_opt_new = FP_struct.fingerprints

# Add the 'source' column to features dataframe since there are duplicate ids
# due to the fact that Chris and I ran the same structures
df_features_post_opt_new["data", "source"] = df_bulk_dft_not_processed["source"]
df_features_post_opt_new["data", "INDEX_OLD"] = df_bulk_dft_not_processed["INDEX_OLD"]

# +
df_features_post_opt_comb = pd.concat([
    df_features_post_opt,
    df_features_post_opt_new])

nan_mask = df_features_post_opt_comb["voronoi"].isnull().any(axis="columns")

df_features_post_opt_comb_cpy = copy.deepcopy(df_features_post_opt_comb)

df_features_post_opt_comb = df_features_post_opt_comb.loc[~nan_mask]
df_nan_in_voro = df_features_post_opt_comb_cpy.loc[nan_mask]

print(
    "Does `df_features_post_opt` have any NaN values in it: ",
    "\n -->",
    df_features_post_opt_comb["voronoi"].isnull().any(axis="columns").any())

# +
old_indices = df_features_post_opt_comb["data", "INDEX_OLD"]

df_features_post_opt_comb = df_features_post_opt_comb.set_index(old_indices)
# df_features_post_opt_comb

# +
index_renamed = df_features_post_opt_comb.index.rename("id_unique")
df_features_post_opt_comb = df_features_post_opt_comb.set_index(index_renamed)
df_features_post_opt_comb = df_features_post_opt_comb.drop(("data", "INDEX_OLD"), axis=1)

index_renamed = df_bulk_dft.index.rename("id_unique")
df_bulk_dft = df_bulk_dft.set_index(index_renamed)
df_bulk_dft = df_bulk_dft.drop(("INDEX_OLD"), axis=1)

# +
# #############################################################################
# with open(df_features_post_opt_path, "wb") as fle:
#     pickle.dump(df_features_post_opt_comb, fle)

# + active=""
#
#
#
# -

assert False

df_features_

# + active=""
#
#
#
#

# + jupyter={}
ab2_ids = df_ids[df_ids["stoich"] == "AB2"]["unique_ids"]

tmp = df_features_pre_opt.loc[ab2_ids].describe().loc["std"].tolist()
print("df_features_pre_opt.shape:", df_features_pre_opt.shape)

len([i for i in tmp if i < 0.00000000001])

print(df_features_pre_opt.loc[:, df_features_pre_opt.loc[ab2_ids].describe().loc["std"] < 0.00001].shape)
tmpa = df_features_pre_opt.loc[:, df_features_pre_opt.loc[ab2_ids].describe().loc["std"] < 0.00001]

tmpa = [i[1] for i in tmpa.columns]

# + jupyter={}
ab2_ids = df_ids[df_ids["stoich"] == "AB3"]["unique_ids"]

tmp = df_features_pre_opt.loc[ab2_ids].describe().loc["std"].tolist()
print("df_features_pre_opt.shape:", df_features_pre_opt.shape)

len([i for i in tmp if i < 0.00000000001])

print(df_features_pre_opt.loc[:, df_features_pre_opt.loc[ab2_ids].describe().loc["std"] < 0.00001].shape)
tmpb = df_features_pre_opt.loc[:, df_features_pre_opt.loc[ab2_ids].describe().loc["std"] < 0.00001]

tmpb = [i[1] for i in tmpb.columns]

# + jupyter={}

tmpa

[i for i in tmpa if i not in tmpb]

len(set(tmpa) & set(tmpb))

len(tmpa)

# + jupyter={}
271 - 170

# 125 - 102

# + active=""
# 125 columns of no info for IrO2 and IrO3
# 102/101 columns of no info for IrO2/3
