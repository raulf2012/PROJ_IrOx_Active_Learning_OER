# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_irox] *
#     language: python
#     name: conda-env-PROJ_irox-py
# ---

# # New ML Active Learning Workflow
# ---

# + active=""
# Following system is failing with voro fingerprinting:
# z39g648rnl
# -

read_from_PROJ_DATA = False
read_from_PROJ_DATA = True

# # Import Modules

# +
import os
print(os.getcwd())
import sys

import copy
import pickle

import pandas as pd

from sklearn.decomposition import PCA


# #############################################################################
from protosearch.ml_modelling.fingerprint import (
    FingerPrint,
    VoronoiFingerprint
    )


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

# with open(static_irox_structures_kirsten_path, "rb") as fle:
#     df_struct_kirsten = pickle.load(fle)

df_ids = pd.read_csv(unique_ids_path)
# -

# # Remove rows with missing atoms objects

# +
# Removing missing data
df_bulk_dft = df_bulk_dft[df_bulk_dft["atoms"].notnull()]
df_struct = df_struct[df_struct["atoms"].notnull()]

df_bulk_dft["INDEX_NEW"] = df_bulk_dft.index + "_" + df_bulk_dft["source"]
df_bulk_dft["INDEX_OLD"] = df_bulk_dft.index
df_bulk_dft = df_bulk_dft.set_index("INDEX_NEW")

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
# -

# # Read Current df_features_pre

if read_from_PROJ_DATA:
    path_i = os.path.join(
        os.environ["PROJ_DATA"], "04_IrOx_surfaces_OER",
        "PROJECT_COMPUTED_OUT_DATA/PROJ_IrOx_Active_Learning_OER",
        "workflow/ml_modelling/voronoi_featurize",
        "out_data/df_features_pre_opt.pickle")
    with open(path_i, "rb") as fle:
        df_features_pre_opt__before = pickle.load(fle)
else:
    try:
        with open(df_features_pre_opt_path, "rb") as fle:
            df_features_pre_opt__before = pickle.load(fle)
    except:
        df_features_pre_opt__before = pd.DataFrame()

# +
ids_to_process = [i for i in df_struct.index if i not in df_features_pre_opt__before.index]

df_struct_to_process = df_struct.loc[ids_to_process]
# -

# # Processing pre-opt Fingerprints

if df_struct_to_process.shape[0] == 0:
    df_features_pre_opt__after = pd.DataFrame()
else:
    FP_struct = FingerPrint(**{
        "feature_methods": ["voronoi"],
        "input_data": df_struct_to_process,
        "input_index": ["atoms"]})

    FP_struct.generate_fingerprints()
    df_features_pre_opt__after = FP_struct.fingerprints

# +
df_features_pre_opt = pd.concat([
    df_features_pre_opt__after,
    df_features_pre_opt__before,
    ])

# #############################################################################
with open(df_features_pre_opt_path, "wb") as fle:
    pickle.dump(df_features_pre_opt, fle)

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

# +
# assert False
# -

# # Getting Features for Post-DFT Structures

if read_from_PROJ_DATA:
    path_i = os.path.join(
        os.environ["PROJ_DATA"], "04_IrOx_surfaces_OER",
        "PROJECT_COMPUTED_OUT_DATA/PROJ_IrOx_Active_Learning_OER",
        "workflow/ml_modelling/voronoi_featurize",
        "out_data/df_features_post_opt.pickle")
    with open(path_i, "rb") as fle:
        df_features_post_opt = pickle.load(fle)

        df_features_post_opt["data", "INDEX_OLD"] = df_features_post_opt.index
        df_features_post_opt["data", "INDEX_NEW"] = df_features_post_opt["data", "INDEX_OLD"] + "_" + df_features_post_opt["data"]["source"]

        df_features_post_opt = df_features_post_opt.set_index(df_features_post_opt["data", "INDEX_NEW"])
        df_features_post_opt = df_features_post_opt.drop(labels=[["data", "INDEX_NEW"]], axis=1)
else:

    try:
        with open(df_features_post_opt_path, "rb") as fle:
            df_features_post_opt = pickle.load(fle)

        df_features_post_opt["data", "INDEX_OLD"] = df_features_post_opt.index
        df_features_post_opt["data", "INDEX_NEW"] = df_features_post_opt["data", "INDEX_OLD"] + "_" + df_features_post_opt["data"]["source"]

        df_features_post_opt = df_features_post_opt.set_index(df_features_post_opt["data", "INDEX_NEW"])
        df_features_post_opt = df_features_post_opt.drop(labels=[["data", "INDEX_NEW"]], axis=1)

        print("df_features_post_opt.shape:", df_features_post_opt.shape)

    except:
        df_features_post_opt = pd.DataFrame()

# +
# df_features_post_opt

# df_bulk_dft.index

# +
ids_to_process = [i for i in df_bulk_dft.index if i not in df_features_post_opt.index]
df_bulk_dft_not_processed = df_bulk_dft.loc[ids_to_process]

print("df_bulk_dft_not_processed.shape:", df_bulk_dft_not_processed.shape)
# df_bulk_dft_not_processed.head()
# -

if df_bulk_dft_not_processed.shape[0] == 0:
    df_features_post_opt_new = pd.DataFrame()
else:
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
# -

#############################################################################
with open(df_features_post_opt_path, "wb") as fle:
    pickle.dump(df_features_post_opt_comb, fle)

print(20 * "# # ")
print("All done!")
assert False

# + active=""
#
#
#
#
#

# + jupyter={"source_hidden": true}
# df_bulk_dft[df_bulk_dft["stoich"] == "AB2"].sort_values("energy_pa")

# # df_features_pre_opt

# # df_features_pre_opt.shape

# # df_features_pre_opt.loc["6dzhcimdxs"]

# "6dzhcimdxs" in df_bulk_dft.index

# assert False

# df_struct.index.unique()


# assert False

# + jupyter={"source_hidden": true}
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

# with open(df_features_pre_opt_kirsten_path, "rb") as fle:
#     df_features_pre_opt_kirsten = pickle.load(fle)

# + jupyter={"source_hidden": true}
# df_features_pre_opt

# df_features_pre_opt__before.head()

# + jupyter={"source_hidden": true}
# sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "workflow/ml_modelling"))
# from ml_methods import get_ml_dataframes
# ml_data_dict = get_ml_dataframes(names=["df_features_pre_opt_path"])
# df_features_pre_opt = ml_data_dict["df_features_pre_opt"]

# + jupyter={"source_hidden": true}
# TEMP
# df_struct = df_struct.sample(n=23)

# + jupyter={"source_hidden": true}
# with open(df_features_pre_opt_path, "rb") as fle:
#     df_features_pre_opt = pickle.load(fle)
