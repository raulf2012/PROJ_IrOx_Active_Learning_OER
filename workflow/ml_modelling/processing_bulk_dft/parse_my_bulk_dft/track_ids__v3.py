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

# +
import os
import sys
import pickle

import numpy as np
import pandas as pd
pd.set_option('display.max_rows', None)

# #############################################################################
sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import (
    bulk_dft_data_path, unique_ids_path,
    static_irox_structures_path)

# #############################################################################
from IPython.display import display

# +
from out_data.inputs_nersc_iro2 import ignore_ids as ignore_ids_nersc_iro2
# from out_data.inputs_nersc_iro3 import ignore_ids as ignore_ids_nersc_iro3
ignore_ids_nersc_iro3 = []

from out_data.inputs_sher_iro2 import ignore_ids as ignore_ids_sher_iro2
from out_data.inputs_sher_iro3 import ignore_ids as ignore_ids_sher_iro3

from out_data.inputs_slac_iro2 import ignore_ids as ignore_ids_slac_iro2
# from out_data.inputs_slac_iro3 import ignore_ids as ignore_ids_slac_iro3
ignore_ids_slac_iro3 = []

ignore_ids_iro2 = ignore_ids_nersc_iro2 + ignore_ids_sher_iro2 + ignore_ids_slac_iro2
ignore_ids_iro3 = ignore_ids_nersc_iro3 + ignore_ids_sher_iro3 + ignore_ids_slac_iro3

print("len(ignore_ids_iro2):", len(ignore_ids_iro2))
print("len(ignore_ids_iro3):", len(ignore_ids_iro3))

ignore_ids_dict = {
    "AB2": ignore_ids_iro2,
    "AB3": ignore_ids_iro3}


# Pickling data ######################################################
import os; import pickle
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)
with open(os.path.join(directory, "ignore_ids.pickle"), "wb") as fle:
    pickle.dump(ignore_ids_dict, fle)
# #####################################################################

# +
from out_data.ids_to_run_nersc_iro2 import ids_to_run as ids_to_run_nersc_iro2
from out_data.ids_to_run_nersc_iro3 import ids_to_run as ids_to_run_nersc_iro3

from out_data.ids_to_run_slac_iro2 import ids_to_run as ids_to_run_slac_iro2
# from out_data.ids_to_run_slac_iro3 import ids_to_run as ids_to_run_slac_iro3
ids_to_run_slac_iro3 = []

from out_data.ids_to_run_sher_iro2 import ids_to_run as ids_to_run_sher_iro2
from out_data.ids_to_run_sher_iro3 import ids_to_run as ids_to_run_sher_iro3

ids_to_run_iro2 = ids_to_run_nersc_iro2 + ids_to_run_sher_iro2 + ids_to_run_slac_iro2
ids_to_run_iro3 = ids_to_run_nersc_iro3 + ids_to_run_sher_iro3 + ids_to_run_slac_iro3

ids_to_run_iro3 = list(set(ids_to_run_iro3))
ids_to_run_iro2 = list(set(ids_to_run_iro2))

ids_to_run_sher = ids_to_run_sher_iro2 + ids_to_run_sher_iro3
ids_to_run_nersc = ids_to_run_nersc_iro2 + ids_to_run_nersc_iro3
ids_to_run_slac = ids_to_run_slac_iro2 + ids_to_run_slac_iro3

ids_to_run_sher = list(set(ids_to_run_sher))
ids_to_run_nersc = list(set(ids_to_run_nersc))
ids_to_run_slac = list(set(ids_to_run_slac))
# -

# # Script Inputs

stoich_i = "AB2"

# # Read Data

# +
with open(bulk_dft_data_path, "rb") as fle:
    df_bulk_dft = pickle.load(fle)
    df_bulk_dft = df_bulk_dft[(df_bulk_dft["source"] == "raul")]

with open(static_irox_structures_path, "rb") as fle:
    df_static = pickle.load(fle)
    df_static = df_static[(df_static["source"] == "chris")]

df_ids = pd.read_csv(unique_ids_path)

with open("out_data/df_new_jobs.pickle", "rb") as fle:
    df_new_jobs = pickle.load(fle)

with open("out_data/df_irox_long.pickle", "rb") as fle:
    df_irox_long = pickle.load(fle)
# -

unique_ids_in_process = df_new_jobs[df_new_jobs["stoich"] == stoich_i]["id"]
unique_ids_in_process = unique_ids_in_process.unique().tolist()
unique_ids_in_process = [i for i in unique_ids_in_process if i.isdigit()]


# +
def method(row_i):
    atoms = row_i["atoms"]
    num_atoms = atoms.get_number_of_atoms()
    return(num_atoms)

df_static["num_atoms"] = df_static.apply(method, axis=1)
df_static = df_static.drop(["atoms", "path", "source",], axis=1)


# + active=""
#
#
#
#
#

# +
def method(row_i):
    new_column_values_dict = {}

    stoich_i = row_i["stoich"]

    id_old = row_i["id_old"]
    id_old_str = str(id_old).zfill(3)

    # #########################################################################
    if row_i.name in df_bulk_dft.index:
        new_column_values_dict["done"] = True
    else:
        new_column_values_dict["done"] = False

    # #########################################################################
    if id_old_str in unique_ids_in_process:
        new_column_values_dict["being_processed"] = True
    else:
        new_column_values_dict["being_processed"] = False

    # #########################################################################
    if stoich_i == "AB2":
        ignore_ids_i = ignore_ids_iro2
    elif stoich_i == "AB3":
        ignore_ids_i = ignore_ids_iro3
    # id_old = str(row_i["id_old"]).zfill(3)

    if id_old_str in ignore_ids_i:
        new_column_values_dict["ignored"] = True
    else:
        new_column_values_dict["ignored"] = False

    # #########################################################################
    # stoich_i = row_i["stoich"] == "AB2"
    if stoich_i == "AB2":
        ids_to_run_i = ids_to_run_iro2
    elif stoich_i == "AB3":
        ids_to_run_i = ids_to_run_iro3

    if id_old in ids_to_run_i:
        new_column_values_dict["ids_to_run"] = True
    else:
        new_column_values_dict["ids_to_run"] = False


    # #########################################################################
    if stoich_i == "AB2":
        if id_old in ids_to_run_sher_iro2:
            new_column_values_dict["cluster"] = "sherlock"
        elif id_old in ids_to_run_nersc_iro2:
            new_column_values_dict["cluster"] = "nersc"
        elif id_old in ids_to_run_slac_iro2:
            new_column_values_dict["cluster"] = "slac"
        else:
            new_column_values_dict["cluster"] = np.nan

    elif stoich_i == "AB3":
        if id_old in ids_to_run_sher_iro3:
            new_column_values_dict["cluster"] = "sherlock"
        elif id_old in ids_to_run_nersc_iro3:
            new_column_values_dict["cluster"] = "nersc"
        elif id_old in ids_to_run_slac_iro3:
            new_column_values_dict["cluster"] = "slac"
        else:
            new_column_values_dict["cluster"] = np.nan


    # #########################################################################
    df_new_jobs_i = df_new_jobs[
        (df_new_jobs["stoich"] == stoich_i) & \
        (df_new_jobs["id"] == str(id_old).zfill(3))
        ]

    if len(df_new_jobs_i) == 1:
        row_j = df_new_jobs_i.iloc[0]

        action_j = row_j["action"]

        bool_0 = action_j == 'Job is busy, will skip'
        bool_1 = action_j == 'Time out or failed | Restarting isif 3 calc'
        bool_2 = action_j == 'Time out or failed | Restarting isif 7 calc'
        bool_3 = action_j == 'Job done | ISIF 3 done | --> isif 2'
        
        if bool_0 or bool_1 or bool_2 or bool_3:
            new_column_values_dict["running"] = True
        else:
            new_column_values_dict["running"] = False

        bool_0 = action_j == "Error, need manual attention"
        bool_1 = action_j == "Couldn't figure out what to do"
        if bool_0 or bool_1:
            new_column_values_dict["error"] = True

    elif len(df_new_jobs_i) == 0:
        new_column_values_dict["running"] = np.nan
    elif len(df_new_jobs_i) > 1:
        # print(df_new_jobs_i)
        # display(df_new_jobs_i)
        new_column_values_dict["running"] = "TEMP, more than 1 sys"

    # #########################################################################
    for key, value in new_column_values_dict.items():
        row_i[key] = value
    return(row_i)

df_i = df_static
df_static = df_i.apply(method, axis=1)

# + active=""
#
#
#
#
#
#
# -

# 697 total AB2  structures in df_static
#
# ------------------------------------
# 87 systems with num_atoms >= 100
#
# 132 systems with num_atoms >= 75

# +
697 # Total AB2
565 # AB2 | < 75 atoms
105 # AB2 | < 75 atoms | Not  done
# 32 # AB2 | < 75 atoms | Not  done | Not ignored
18 # AB2 | < 75 atoms | Not  done | Not ignored


# #####################################
14  # Currently running
7  # That have identified error

# +
df = df_static.drop("static_id", axis=1)
df = df[
    (df["stoich"] == stoich_i) & \
    (df["num_atoms"] < 75) & \
#     (df["num_atoms"] > 75) & \
#     (df["num_atoms"] < 100) & \
#     (df["being_processed"] == False) & \
    (df["done"] == False) & \
    (df["ignored"] == False) & \
#     (df["running"] == True) & \
    [True for i in range(len(df))]
    ]
print("df.shape:", df.index.unique().shape)
# display(df)



from running_ids import sherlock_ids_running, slac_ids_running


slac_ids_running

print(len([i for i in df["id_old"].tolist() if str(i).zfill(3) in sherlock_ids_running]))
print(len([i for i in df["id_old"].tolist() if str(i).zfill(3) in slac_ids_running]))

# sherlock_ids_running
# df["id_old"].tolist()

tmp_list_0 = np.sort(
[str(i).zfill(3) for i in df["id_old"].tolist() if str(i).zfill(3) not in sherlock_ids_running + slac_ids_running]
)
print(tmp_list_0)

# +
df_tmp = df_new_jobs[
    (df_new_jobs["id"].isin(tmp_list_0)) & \
    (df_new_jobs["stoich"] == "AB2")
    ].sort_values("id")
display(df_tmp); print("")

df = df_static
df_static[
    (df["id_old"].isin(tmp_list_0)) & \
    (df["stoich"] == "AB2") & \
    (df["error"] != True)
    ].sort_values("id_old")
# -



# 60,
# 89,
# 118,
155,
193,
# 236,
# 250,
275,
303,
# 483,
# 649

# +
repeated_ids_list = []
for i in ids_master_list:
    if ids_master_list.count(i) > 1:
        repeated_ids_list.append(i)


set(repeated_ids_list)
