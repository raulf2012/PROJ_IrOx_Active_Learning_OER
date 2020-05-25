# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python [conda env:research-new]
#     language: python
#     name: conda-env-research-new-py
# ---

# # Import Modules

# + jupyter={}
import os
import sys
import pickle

import pandas as pd

# #############################################################################
sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import (
    bulk_dft_data_path, unique_ids_path,
    static_irox_structures_path)

# #############################################################################
from IPython.display import display

# +
with open(bulk_dft_data_path, "rb") as fle:
    df_bulk_dft = pickle.load(fle)

with open(static_irox_structures_path, "rb") as fle:
    df_static = pickle.load(fle)
    
df_ids = pd.read_csv(unique_ids_path)


# +
def method(row_i):
    atoms = row_i["atoms"]
    num_atoms = atoms.get_number_of_atoms()
    return(num_atoms)

df_static["num_atoms"] = df_static.apply(method, axis=1)

# +
row_i = df_static.iloc[0]

atoms = row_i["atoms"]
num_atoms = atoms.get_number_of_atoms()
# -

ab3_ids = df_ids[df_ids["stoich"] == "AB3"]["unique_ids"].tolist()

# +
print("df_static.shape:", df_static.shape)

print("df_bulk_dft.shape:", df_bulk_dft.shape)
# -

# # IrO3 Structures

df_static_ab3 = df_static[
    (df_static["stoich"] == "AB3") & \
    (df_static["source"] == "chris") & \
    [True for i in range(len(df_static))]
    ]

# +
df_bulk_dft_ab3 = df_bulk_dft[
    (df_bulk_dft["stoich"] == "AB3") & \
    (df_bulk_dft["source"] == "raul") & \
    [True for i in range(len(df_bulk_dft))]
    ]

# df_bulk_dft_ab3

# +
df_bulk_dft_ab3.index

ab3_ids_not_finished = [i for i in df_static_ab3.index if i not in df_bulk_dft_ab3.index]

df_i = df_static_ab3.loc[ab3_ids_not_finished]

print("Structures with greater than 100 atoms")
df_tmp = df_i[df_i["num_atoms"] >= 100]

display(df_tmp)

# +
df_i = df_i.drop(df_tmp.index)

# df_i

# +
# #############################################################################
path_i = os.path.join(
    os.environ["PROJ_DATA"],
    "04_IrOx_surfaces_OER/ml_bulk_irox_dft/iro3",
    "df_dict.pickle")
with open(path_i, "rb") as fle:
    df_dict = pickle.load(fle)
df_iro3_sherlock = df_dict["df"]
df_new_jobs_sherlock = df_dict["df_new_jobs"]
df_iro3_sherlock["cluster"] = "sherlock"

# #############################################################################
# Parsing Sherlock IrO3 DFT Data
path_i = os.path.join(
    os.environ["PROJ_DATA"],
    "04_IrOx_surfaces_OER/ml_bulk_irox_dft/iro3",
    "df_dict_sher.pickle")
with open(path_i, "rb") as fle:
    df_dict = pickle.load(fle)
df_iro3_nersc = df_dict["df"]
df_new_jobs_nersc = df_dict["df_new_jobs"]
df_iro3_nersc["cluster"] = "nersc"

# #############################################################################
df_iro3 = pd.concat([
    df_iro3_nersc,
    df_iro3_sherlock 
    ], axis=0, sort=True)

# +
ids_running_on_nersc = [
    46,
    70,
    83,
    144,
    152,
    173,
    206,
    220,
    228,
    ]

ids_running_on_sher = [
    2,
    22,
    51,
    101,
    250,
    ]

ids_running = ids_running_on_nersc + ids_running_on_sher

len(ids_running)

# + active=""
# 236
# -

df_i = df_i[~df_i["id_old"].isin(ids_running)]

# +
# df_iro3

# for id_i in df_i["id_old"].tolist():
#     df_tmp = df_iro3[df_iro3["id"] == str(id_i).zfill(3)]

#     print(id_i)
#     display(df_tmp.sort_values("revision"))
#     print(40 * "-")
    
    
completed_ids = df_iro3[
    (df_iro3["completed"] == True) & \
    (df_iro3["isif"] == 2)
    ]["id"].tolist()

not_done_ids = df_iro3[~df_iro3["id"].isin(completed_ids)]["id"].unique().tolist()
for ind_i, id_i in enumerate(not_done_ids):
    df_tmp = df_iro3[df_iro3["id"] == str(id_i).zfill(3)]

    print(ind_i, " | ", id_i)
    display(df_tmp.sort_values("revision"))
    print(40 * "-")


# + active=""
#
#
#
#

# + active=""
#
#
#
#
# -



# +
# type(df_iro3.loc[22].iloc[0]["completed"])

# df_iro3["completed"].isnull()
