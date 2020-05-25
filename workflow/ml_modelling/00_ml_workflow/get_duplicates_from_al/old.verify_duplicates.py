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

# + [markdown] Collapsed="false"
# # Import Modules

# + Collapsed="false"
import os
import sys                                 

import pickle

import pandas as pd

# #########################################################
# Local Imports
sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import df_dij_path

from al_data import main_AB2_run, main_AB3_run
# -

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "workflow/ml_modelling"))
from ml_methods import get_ml_dataframes

# # Script Inputs

# +
stoich_i = "AB2"

if stoich_i == "AB2":
    AL_data_path = main_AB2_run
elif stoich_i == "AB3":
    AL_data_path = main_AB3_run

# + [markdown] Collapsed="false"
# # Import Data

# +
with open(AL_data_path, "rb") as fle:
    AL_i = pickle.load(fle)

# #########################################################
DF_dict = get_ml_dataframes(
    names=[
        'bulk_dft_data_path',
        'unique_ids_path',
        'prototypes_data_path',
        'static_irox_structures_path',
        'static_irox_structures_kirsten_path',
        'oqmd_irox_data_path',
        'df_features_pre_opt_path',
        'df_features_pre_opt_kirsten_path',
        'df_features_post_opt_path',
        'oer_bulk_structures_path',
        'df_ccf_path',
        'df_dij_path',
        'ids_to_discard__too_many_atoms_path',
        'df_prototype_dft_path',
        'df_prototype_static_path',
        ]
    )

df_bulk_dft = DF_dict["bulk_dft_data"]
# unique_ids = DF_dict["unique_ids"]
# prototypes_data = DF_dict["prototypes_data"]
# static_irox_structures = DF_dict["static_irox_structures"]
# static_irox_structures_kirsten = DF_dict["static_irox_structures_kirsten"]
# oqmd_irox_data = DF_dict["oqmd_irox_data"]
# df_features_pre_opt = DF_dict["df_features_pre_opt"]
# df_features_pre_opt_kirsten = DF_dict["df_features_pre_opt_kirsten"]
# df_features_post_opt = DF_dict["df_features_post_opt"]
# oer_bulk_structures = DF_dict["oer_bulk_structures"]
# df_ccf = DF_dict["df_ccf"]
df_dij = DF_dict["df_dij"]
# ids_to_discard__too_many_atoms = DF_dict["ids_to_discard__too_many_atoms"]
# df_dft_final_final = DF_dict["df_dft_final_final"]

df_prototype_dft_path = DF_dict["df_prototype_dft"]
df_prototype_static_path = DF_dict["df_prototype_static"]

# +
from misc_modules.pandas_methods import drop_columns

df_prototype_dft_path = drop_columns(df=df_prototype_dft_path, columns=[
    # 'p_name',
    'structure_name',
    # 'spacegroup',
    'wyckoffs',
    'species',
    ],
    keep_or_drop="drop")

df_bulk_dft = drop_columns(df=df_bulk_dft, columns=[
    "atoms",
    "form_e_chris",
    "id",
    "id_old",
    "path",
    "volume",
    # "",
    ], keep_or_drop="drop")

# #########################################################
shared_indices = df_bulk_dft.index.intersection(df_prototype_dft_path.index)
df_info = pd.concat([
    df_prototype_dft_path.loc[shared_indices],
    df_bulk_dft[df_bulk_dft.source == "raul"].loc[shared_indices],
    ],
    axis=1)

df_info = df_info.sort_values("dH")
df_info = df_info.drop_duplicates()

# + Collapsed="false" active=""
#
#
#

# + Collapsed="false"
last_gen_key = list(AL_i.al_gen_dict.keys())[-1]
AL_gen_f = AL_i.al_gen_dict[last_gen_key]

model = AL_gen_f.model
model_notdupl = model[model.duplicate == False]

indices_not_duplicates = model_notdupl.index

# + Collapsed="false"
# Subset of df_dij that contains entries of final AL generation that aren't duplicates
df_dij_sub = df_dij.loc[indices_not_duplicates, indices_not_duplicates]
# df_dij_sub[df_dij_sub < 0.01]


for i_cnt, row_i in df_dij_sub.iterrows():
    id_i = row_i.name

    # duplicate_ids = row_i[row_i < 0.02].index.tolist()
    duplicate_ids = row_i[row_i < 0.03].index.tolist()
    if id_i in duplicate_ids:
        duplicate_ids.remove(id_i)

    if len(duplicate_ids) > 0:
        print(id_i, ":", duplicate_ids)

# +
df_dij.loc["8p8evt9pcg", "xw9y6rbkxr"]

tmp_list = model_notdupl.sort_values("y_real").iloc[0:15].index.tolist()

df_dij.loc[tmp_list, tmp_list]

# + active=""
#
#
#
#

# +
# 9yz2mt8hbh

# +
model_notdupl.sort_values("y_real")

# id_i = "64cg6j9any"
# id_i = "n36axdbw65"
# id_i = "clc2b1mavs"
id_i = "ck638t75z3"
# id_i = "mkbj6e6e9p"
# id_i = "b49kx4c19q"
# id_i = "85z4msnl6o"
# id_i = "bpc2nk6qz1"


# df_dij[df_dij.loc[id_i] < 0.01]
# df_dij.loc[id_i][df_dij.loc[id_i] < 0.2]
# df_dij.loc[id_i] < 0.4
# df_dij.loc[id_i]

# df_i = df_dij.loc[id_i][df_dij.loc[id_i] < 0.05]
df_i = df_dij.loc[id_i][df_dij.loc[id_i] < 0.08]

df_i.index.tolist()

# +
# df_prototype_dft_path.loc[df_i.index]
# df_bulk_dft.loc[df_i.index].sort_values("dH")

df_info.loc[df_i.index].sort_values("dH")

# +


# df_info.shape
# -

# model_notdupl.sort_values("y_real").iloc[0:10].index.tolist()
model_notdupl.sort_values("y_real").iloc[0:10]

# + active=""
#
#
#

# + jupyter={}
# [i for i in df_dij.index if "byna" in i]

# + Collapsed="false" jupyter={}
# model_notdupl.sort_values("y_real")

# model.sort_values("y_real").iloc[0:15]

# + Collapsed="false" jupyter={}
# /mnt/f/Dropbox/01_norskov/00_git_repos/PROJ_IrOx_Active_Learning_OER/workflow/ml_modelling/00_ml_workflow/191102_new_workflow/00_abx_al_runs/out_data/AB2/gp_ucb_True/TEST_AL_muvasubu.pickle

# # #############################################################################
# # AL_data_path += "/AB3/gp_ucb_True"
# # AL_data_path += "/AL_geheneva.pickle"

# # #####################################
# AL_data_path += "/AB2/gp_ucb_True"
# # #####################################
# # AL_data_path += "/AL_piritapo.pickle"
# # AL_data_path += "/TEST_AL_wakuhifa.pickle"
# AL_data_path += "/TEST_AL_muvasubu.pickle"

# + Collapsed="false" jupyter={}
# # AL_data_path = "/home/raulf2012/Dropbox/01_norskov/00_git_repos/PROJ_IrOx_Active_Learning_OER/workflow/ml_modelling/00_ml_workflow/191102_new_workflow/00_abx_al_runs/out_data"
# AL_data_path = os.path.join(
#     os.environ["PROJ_irox"],
#     "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/00_abx_al_runs/out_data",
#     )

# # #####################################
# AL_data_path += "/AB3/gp_ucb_True/01_attempt"
# # #####################################
# AL_data_path += "/AL_geheneva.pickle"

# # AL_geheneva.pickle
