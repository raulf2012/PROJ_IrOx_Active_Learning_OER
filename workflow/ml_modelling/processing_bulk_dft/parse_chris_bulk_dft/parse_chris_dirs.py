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

# # Parsing Chris's DFT Data on NERSC
# ---
#
#
# Author(s): Raul A. Flores

# # Notes
# ---

# + {"active": ""}
# ID IrO3-146 isn't availabe in the parsed folders
#
# -171.82382373
# Ir8O24
#
# # #############################################################################
# # Have these columns by the end
# ['atoms',
#  'energy',
#  'energy_pa',
#  'force_max',
#  'force_sum',
#  'form_e_chris',
#  'path',
#  'stoich',
#  'id_old',
#  'source']
# -

# # Import Modules

# +
import os
import sys

sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "chris_prototypes_structures/oqmd_iro3",
    ))

sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "data",
    ))

import pickle

# import numpy as np
import pandas as pd
from ase import io
# -

# # Script Inputs

root_dir = os.path.join(
    "/mnt/f/GDrive/norskov_research_storage",
    "00_projects/PROJ_irox_2/chris_nersc_files")

# # Parse ID List from Files from Chris

# +
file_path_i = os.path.join(
    os.environ["PROJ_irox"],
    "data/ml_irox_data/iro2_training_data.csv")
train_data_iro2 = pd.read_csv(file_path_i)
train_data_iro2.set_index("id", inplace=True)

file_path_i = os.path.join(
    os.environ["PROJ_irox"],
    "data/ml_irox_data/iro3_training_data.csv")
train_data_iro3 = pd.read_csv(file_path_i)
train_data_iro3.set_index("id", inplace=True)

train_data_dict = {
    "iro2": train_data_iro2,
    "iro3": train_data_iro3,
    }
# -

# # Parse NERSC DFT Data
#
# Comment out to read pickled data file instead (saves time)

# + {"jupyter": {"source_hidden": true}}
# master_data_list = []
# id_list_nersc = []
# for subdir, dirs, files in os.walk(root_dir):
#     if "gas_references" in subdir:
#         continue
#     if "IrO2/Old_ML_calcs" in subdir:
#         continue
#     if "__old__" in subdir:
#         continue
#     if "volume" in subdir:
#         continue

#     if "OUTCAR" in files:
#         dir_i = subdir
#         id_i = int(subdir.split("/")[-1].split("_")[0])

#         try:
#             atoms_i = io.read(os.path.join(dir_i, "OUTCAR"))
#         except:
#             atoms_i = None

#         replace_path_snip = os.path.join(
#             "/mnt/f/GDrive/norskov_research_storage"
#             "00_projects/PROJ_irox_2/chris_nersc_files/"
#             )
#         path_short_i = dir_i.replace(replace_path_snip, "")

#         dict_i = {
#             "id_old": id_i,
#             "atoms": atoms_i,
#             "path": path_short_i}

#         master_data_list.append(dict_i)


# # Save Data ###################################################################
# directory = "out_data"
# if not os.path.exists(directory):
#     os.makedirs(directory)

# with open(os.path.join(directory, "parse_data.pickle"), "wb") as fle:
#     pickle.dump(master_data_list, fle)

# +
with open("out_data/parse_data.pickle", "rb") as fle:
    master_data_list = pickle.load(fle)

df = pd.DataFrame(master_data_list)


# -

# # Process dataframe

# +
# #############################################################################
# #############################################################################
def method(row_i):
    """
    """
    if "IrO2" in row_i["path"]:
        sys_i = "AB2"
    elif "IrO3" in row_i["path"]:
        sys_i = "AB3"
    else:
        sys_i = None
    return(sys_i)

df["stoich"] = df.apply(
    method,
    axis=1)

# #############################################################################
# #############################################################################
def method(row_i):
    """
    """
    ignore_list = [
        "IrO3/winnersIrO3",
        "IrO3/full_relax",
        "IrO3/full_relax_ML1",
        "IrO3/full_relax_ML2",
        "IrO3/full_relax_ML3",
        "IrO3/full_relax_ML4",
        "IrO3/single_point",
        "IrO3/volume_relax",
        "IrO3/volume_relax_ML1",
        "IrO3/volume_relax_ML2",
        "IrO3/volume_relax_ML3",
        "IrO3/volume_relax_ML4",
        ]
    ignore = False
    for ignore_seg_i in ignore_list:
        if ignore_seg_i in row_i["path"]:
            ignore = True

    return(ignore)

df["ignore_tag"] = df.apply(
    method,
    axis=1)

# #############################################################################
# #############################################################################
def method(row_i):
    """
    """
    if "volume" in row_i["path"]:
        out = True
    else:
        out = False
    return(out)

df["volume_tag"] = df.apply(
    method,
    axis=1)

# #############################################################################
# #############################################################################


df = df[df["ignore_tag"] == False]
df = df[df["volume_tag"] == False]


good_bye_list = [
    "ignore_tag",
    "volume_tag",
    ]

# df_dft_calcs.drop(good_bye_list, axis=1, inplace=True)
df.drop(good_bye_list, axis=1, inplace=True)
# -

# # IrO2

# +
# df_iro2 = df[df["stoich"] == "IrO2"]
df_iro2 = df[df["stoich"] == "AB2"]

master_data = []
for id_i, row_i in train_data_dict["iro2"].iterrows():

    if row_i["source"] != "chris":
        continue

    form_e_chris_i = row_i["form_e_chris"]

    df_i = df_iro2[df_iro2["id_old"] == id_i]

    if len(df_i) == 0:
        print(id_i, " | There are no rows for this id!!!!")

    df_0 = df_i[df_i["path"].str.contains("final_opt_new1-3")]
    df_1 = df_i[df_i["path"].str.contains("final_relax")]

    row_i = None
    if len(df_0) > 0:
        if len(df_0) > 1:
            print("NOOOOOOOOOO!!!!!!!!!")
        row_j = df_0.iloc[0]

    else:
        if len(df_1) > 0:
            if len(df_1) > 1:
                print("NOOOOOOOOOO!!!!!!!!!")
            row_j = df_1.iloc[0]
        else:
            tmp = 42

    master_data.append(row_j)

df_iro2_unique = pd.concat(master_data, axis=1, sort=True).transpose()
# df_iro2_unique.set_index("id", inplace=True)
# -

# # IrO3

# +
df_iro3 = df[df["stoich"] == "AB3"]

master_data = []
for id_i, row_i in train_data_dict["iro3"].iterrows():

    if row_i["source"] != "chris":
        continue

    form_e_chris_i = row_i["form_e_chris"]

    df_i = df_iro3[df_iro3["id_old"] == id_i]

    if len(df_i) == 0:
        print(id_i, " | There are no rows for this id!!!!")
        row_j = pd.Series({"id_old": int(id_i)})
    else:
        df_0 = df_i[df_i["path"].str.contains("final_opt_new1-3_sorted")]
        if len(df_0) > 0:
            if len(df_0) > 1:
                print("NOOOOOOOOOO!!!!!!!!!")
            row_j = df_0.iloc[0]
        else:
            row_j = df_i.iloc[0]

    master_data.append(row_j)

df_iro3_unique = pd.concat(master_data, axis=1, sort=True).transpose()
# df_iro3_unique = df_iro3_unique.astype({"id": int})
# df_iro3_unique.set_index("id", inplace=True)
# -

# # Combining IrO2 and IrO3 dataframes

# +
df_iro2_dft = df_iro2_unique
df_iro3_dft = df_iro3_unique

df_dft_calcs = pd.concat([
    df_iro2_dft,
    df_iro3_dft,
    ])
# -

# # Save data to pickle

directory = "out_data"
if not os.path.exists(directory):
    os.makedirs(directory)
with open("out_data/df_dft_calcs.pickle", "wb") as fle:
    pickle.dump(df_dft_calcs, fle)

df_dft_calcs[df_dft_calcs["id_old"] == 192].iloc[0]["atoms"].get_potential_energy()

new_id_list = [192, 195, 118, 136, 160, 202, 72, 126, 111, 106, 27]
df_dft_calcs[df_dft_calcs["id_old"].isin(new_id_list)].loc[76]["path"]
