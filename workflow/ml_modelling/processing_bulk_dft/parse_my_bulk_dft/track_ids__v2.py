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

# +
df_static[df_static["id_old"] == 182]
# df_static

df_bulk_dft[df_bulk_dft["id_old"] == 182]

# +
# df_irox_long[df_irox_long["id"] == "250"].sort_values("revision")

df_irox_long[
    (df_irox_long["stoich"] == "AB2") & \
    (df_irox_long["id"] == "250")
    ].sort_values("revision")

df_new_jobs[
    (df_new_jobs["stoich"] == "AB2") & \
    (df_new_jobs["id"] == "250")
    ]
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
row_i =df_static.iloc[0]

id_old = str(row_i["id_old"]).zfill(3)

if id_old in unique_ids_in_process:
    tmp = 42

# +
df_bulk_dft[
    (df_bulk_dft["stoich"] == "AB2") & \
    (df_bulk_dft["source"] == "raul")
    # () & \
    ].index.unique()

#.unique()

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

697 - \
(
    132 + \
    469 + 
    0
    )

565 - (469 + 53)

df_bulk_dft_iro2_ids_list = df_bulk_dft[
    (df_bulk_dft["stoich"] == "AB2") & \
    (df_bulk_dft["source"] == "raul")
    # () & \
    ].index.unique().tolist()

df_static.shape

# + jupyter={}
# df = df_static.drop("static_id", axis=1)
# df = df[
#     (df["running"] == "TEMP, more than 1 sys") & \
# #     (df["stoich"] == stoich_i) & \
# #     (df["num_atoms"] < 75) & \
# #     (df["num_atoms"] > 75) & \
# #     (df["num_atoms"] < 100) & \
# #     (df["being_processed"] == False) & \
# #     (df["done"] == False) & \
# #     (df["ignored"] == False) & \
# #     (df["running"] == True) & \
#     [True for i in range(len(df))]
#     ]
# print("df.shape:", df.index.unique().shape)
# # display(df)

# ids_in_2_clusters = df["id_old"].tolist()

# ids_in_2_clusters


# for id_i in ids_in_2_clusters:
#     tmp = 42

#     id_i = str(id_i).zfill(3)

#     print(
#         "id_i:", id_i,
#         "|",
#         df_irox_long[df_irox_long["id"] == id_i]["source"].unique().tolist()    
#         )

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

np.sort(
[i for i in df["id_old"].tolist() if str(i).zfill(3) not in sherlock_ids_running + slac_ids_running]
)

# -

len(sherlock_ids_running + slac_ids_running)

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
# df_irox_long.head()

# idf_irox_long[df_irox_long["id"] == "155"]
# -

df_static.head()
df_static[df_static["id_old"] == 155]

# +
list_tmp = [i for i in df["id_old"].tolist() if str(i).zfill(3) not in list(sherlock_ids_running + slac_ids_running)]

len(list_tmp)

list_tmp

df_new_jobs[
    (df_new_jobs["stoich"] == "AB2") & \
    (df_new_jobs["id"].isin([str(i).zfill(3) for i in list_tmp]))
    # (df_new_jobs["id"].isin(["250", "275"]))
    # (df_new_jobs["id"] == "250")
    ].sort_values("id")


# df_new_jobs[df_new_jobs["id"] == "250"]

# +
df_static.head()

# df_static[df_static["id_old"] == "182"]
# -

[str(i).zfill(3) for i in list_tmp]

# +
# df_static.head()

# df_static[df_static["id_old"] == 250]

# [i for i in df.index.unique().tolist() if i not  in df_bulk_dft_iro2_ids_list]

# df_bulk_dft.loc['vwxrnun48g']

# 'vwxrnun48g' in df_bulk_dft_iro2_ids_list

# + active=""
#
#
#
#
#
#

# +
# df_ab2_sher = df_irox_long[df_irox_long["source"] == "sherlock"]

# +
ids_master_list = []

df_ab2_i = df_irox_long[
    (df_irox_long["source"] == "sherlock") & \
    (df_irox_long["stoich"] == "AB2")
    ]

df_ab2_i_complete = df_ab2_i[
    (df_ab2_i["isif"] == 2.) & \
    (df_ab2_i["completed"] == True)
    ]
print(
    df_ab2_i_complete["id"].shape,
    df_ab2_i_complete["id"].unique().shape,
    )
ids_i = df_ab2_i_complete["id"].unique().tolist()
ids_master_list.extend(ids_i)

print("sherlock:", df_ab2_i_complete.shape)



# #############################################################################
df_ab2_i = df_irox_long[
    (df_irox_long["source"] == "slac") & \
    (df_irox_long["stoich"] == "AB2")
    ]

df_ab2_i_complete = df_ab2_i[
    (df_ab2_i["isif"] == 2.) & \
    (df_ab2_i["completed"] == True)
    ]
print(
    df_ab2_i_complete["id"].shape,
    df_ab2_i_complete["id"].unique().shape,
    )
ids_i = df_ab2_i_complete["id"].unique().tolist()
ids_master_list.extend(ids_i)

print("slac:", df_ab2_i_complete.shape)

# #############################################################################
df_ab2_i = df_irox_long[
    (df_irox_long["source"] == "nersc") & \
    (df_irox_long["stoich"] == "AB2")
    ]

df_ab2_i_complete = df_ab2_i[
    (df_ab2_i["isif"] == 2.) & \
    (df_ab2_i["completed"] == True)
    ]
print(
    df_ab2_i_complete["id"].shape,
    df_ab2_i_complete["id"].unique().shape,
    )
ids_i = df_ab2_i_complete["id"].unique().tolist()
ids_master_list.extend(ids_i)

print("nersc:", df_ab2_i_complete.shape)
# -



# +
repeated_ids_list = []
for i in ids_master_list:
    # print(i)
    if ids_master_list.count(i) > 1:
#         print(i)
        repeated_ids_list.append(i)


set(repeated_ids_list)

# +
# # df_irox_long.head()

# df_tmp = df_irox_long[
#     (df_irox_long["id"] == "303") & \
#     (df_irox_long["stoich"] == "AB2")
#     ].sort_values(["source", "revision"])

# print(df_tmp.loc[745]["atoms"].get_potential_energy())
# print(df_tmp.loc[1475]["atoms"].get_potential_energy())

# df_tmp

# +
len(ids_master_list)

len(set(ids_master_list))

ids_master_list[0:4]

ids_master_list = [int(i) for i in ids_master_list]

# +
# df_bulk_dft[df_bulk_dft["id_old"] == 338]

# df_bulk_dft

# +
ids_list_from_main_script = df_bulk_dft[
    (df_bulk_dft["stoich"] == "AB2") & \
    (df_bulk_dft["source"] == "raul")
    ]["id_old"].tolist()


ids_not_in_main_list = [i for i in ids_master_list if i not in ids_list_from_main_script]


ids_not_in_main_list
# -

import pickle; import os
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/processing_bulk_dft/parse_my_bulk_dft/out_data",
    "parse_data.pickle"
    )
with open(path_i, "rb") as fle:
    job_dirs, completed_ids = pickle.load(fle)

# + jupyter={}
# print(len(job_dirs))
# print(len(completed_ids))

# # [i for i in ids_0 if i not in job_dirs]

# + jupyter={}
# df_static[df_static["num_atoms"] >= 100]

# + jupyter={}
# # (464,)

# unique_ids_in_process = df_new_jobs[df_new_jobs["stoich"] == "AB2"]["id"]
# unique_ids_in_process = unique_ids_in_process.unique().tolist()

# [i for i in unique_ids_in_process if i.isdigit()]

# # unique_ids_in_process[0].isdigit()

# + jupyter={}
# all_ab2_ids = df_static[
#     (df_static["stoich"] == "AB2") & \
#     [True for i in range(len(df_static))]
#     ].index.tolist()


# ab2_completed_ids = df_static[
#     (df_static["done"] == True) & \
# #     (df_static["ignored"] == True) & \
# #     (df_static["num_atoms"] >= 100) & \
# #     (df_static["ids_to_run"] == False) & \
# #     (df_static["running"] == True) & \
# #     (df_static["error"] == True) & \
#     (df_static["stoich"] == "AB2") & \
#     [True for i in range(len(df_static))]
#     ].index.tolist()

# ab2_errored_ids = df_static[
# #     (df_static["done"] == True) & \
# #     (df_static["ignored"] == True) & \
# #     (df_static["num_atoms"] >= 100) & \
# #     (df_static["ids_to_run"] == False) & \
# #     (df_static["running"] == True) & \
#     (df_static["error"] == True) & \
#     (df_static["stoich"] == "AB2") & \
#     [True for i in range(len(df_static))]
#     ].index.tolist()


# ab2_large_atoms_ids = df_static[
# #     (df_static["done"] == True) & \
# #     (df_static["ignored"] == True) & \
#     (df_static["num_atoms"] >= 100) & \
# #     (df_static["ids_to_run"] == False) & \
# #     (df_static["running"] == True) & \
# #     (df_static["error"] == True) & \
#     (df_static["stoich"] == "AB2") & \
#     [True for i in range(len(df_static))]
#     ].index.tolist()


# ab2_ids_accounted_for = ab2_completed_ids + ab2_errored_ids + ab2_large_atoms_ids


# ab2_ids_not_acc = [i for i in all_ab2_ids if i not in ab2_ids_accounted_for]
# # [i for i all_ab3_ids]

# # df_new_jobs

# df_static.loc[ab2_ids_not_acc]

# + jupyter={}
# all_ab3_ids = df_static[
#     (df_static["stoich"] == "AB3") & \
#     [True for i in range(len(df_static))]
#     ].index.tolist()


# ab3_completed_ids = df_static[
#     (df_static["done"] == True) & \
# #     (df_static["ignored"] == True) & \
# #     (df_static["num_atoms"] >= 100) & \
# #     (df_static["ids_to_run"] == False) & \
# #     (df_static["running"] == True) & \
# #     (df_static["error"] == True) & \
#     (df_static["stoich"] == "AB3") & \
#     [True for i in range(len(df_static))]
#     ].index.tolist()

# ab3_errored_ids = df_static[
# #     (df_static["done"] == True) & \
# #     (df_static["ignored"] == True) & \
# #     (df_static["num_atoms"] >= 100) & \
# #     (df_static["ids_to_run"] == False) & \
# #     (df_static["running"] == True) & \
#     (df_static["error"] == True) & \
#     (df_static["stoich"] == "AB3") & \
#     [True for i in range(len(df_static))]
#     ].index.tolist()


# ab3_large_atoms_ids = df_static[
# #     (df_static["done"] == True) & \
# #     (df_static["ignored"] == True) & \
#     (df_static["num_atoms"] >= 100) & \
# #     (df_static["ids_to_run"] == False) & \
# #     (df_static["running"] == True) & \
# #     (df_static["error"] == True) & \
#     (df_static["stoich"] == "AB3") & \
#     [True for i in range(len(df_static))]
#     ].index.tolist()


# ab3_ids_accounted_for = ab3_completed_ids + ab3_errored_ids + ab3_large_atoms_ids


# ab3_ids_not_acc = [i for i in all_ab3_ids if i not in ab3_ids_accounted_for]
# # [i for i all_ab3_ids]

# # df_new_jobs

# df_static.loc[ab3_ids_not_acc]

# + jupyter={}
# df_new_jobs["action"].unique().tolist()

# [
#  'ALL DONE! | ISIF 2',
#  'Ignoring this id',
#  'Job is busy, will skip',
#  'Error, need manual attention',
#  "Couldn't figure out what to do",
#  'Time out or failed | Restarting isif 3 calc',
#  'Time out or failed | Restarting isif 7 calc',
#  'Job done | ISIF 3 done | --> isif 2',
# ]

# [
#  'ALL DONE! | ISIF 2',
#  'Ignoring this id',
#  'Job is busy, will skip',
#  'Error, need manual attention',
#  "Couldn't figure out what to do",
#  'Time out or failed | Restarting isif 3 calc',
#  'Time out or failed | Restarting isif 7 calc',
#  'Job done | ISIF 3 done | --> isif 2',
# ]

# "Error, need manual attention"
# "Couldn't figure out what to do"

# + jupyter={}
# row_i = df_static.iloc[1]

# id_old = row_i["id_old"]

# id_old in ids_to_run_nersc

# row_i
# stoich_i = row_i["stoich"]
# id_old = row_i["id_old"]


# df_new_jobs_i = df_new_jobs[
#     (df_new_jobs["stoich"] == stoich_i) & \
#     (df_new_jobs["id"] == str(id_old).zfill(3))
#     ]

# assert len(df_new_jobs_i) == 1, "SIFISDF"

# row_j = df_new_jobs_i.iloc[0]

# action_j = row_j["action"]

# bool_0 = action_j == 'Job is busy, will skip'
# bool_1 = action_j == 'Time out or failed | Restarting isif 3 calc'
# bool_2 = 'Time out or failed | Restarting isif 7 calc'
# bool_3 = 'Job done | ISIF 3 done | --> isif 2'
# if bool_0 or bool_1 or bool_2 or bool_3:
#     tmp = 42

# + jupyter={}
# ignore_ids_iro2

# [i for i in set(ignore_ids_iro2) if i not in ignore_ids_iro2]




# for i in ignore_ids_iro2:
    
#     if ignore_ids_iro2.count(i) > 1:
#         print(i)
# "054"

# "054" in ignore_ids_sher_iro2

# "054" in ignore_ids_nersc_iro2


# ignore_ids_sher_iro2.count("054")

# + jupyter={}
# 236 + 96 + 115

# + jupyter={}
# ids_0 = df_ab2_i["id"].unique()

# ids_0.shape
# # df_ab2_i["id"].unique()
