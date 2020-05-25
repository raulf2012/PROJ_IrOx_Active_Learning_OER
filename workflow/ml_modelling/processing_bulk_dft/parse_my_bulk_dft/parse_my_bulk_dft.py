# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_irox] *
#     language: python
#     name: conda-env-PROJ_irox-py
# ---

# # Parsing my bulk DFT data
# ---
#
# 130 calculations with the previous criteria by which jobs are included
# 103 after requiring the job to be done 100%

# # Import Modules

# +
import os
import sys

import pickle

from ase import io
import pandas as pd

from misc_modules.pandas_methods import drop_columns
# -

# # Notebook Prep

# Creating 'out_data' dir
directory = "out_data"
if not os.path.exists(directory):
    os.makedirs(directory)

# # IrO2 Bulk Data

# +
            # #############################################################################
# Parsing Sherlock IrO2 DFT Data
path_i = os.path.join(
    os.environ["PROJ_DATA"],
    "04_IrOx_surfaces_OER/ml_bulk_irox_dft/iro2",
    "df_dict_nersc.pickle")
with open(path_i, "rb") as fle:
    df_dict = pickle.load(fle)
    df_new_jobs_nersc_iro2 = df_dict["df_new_jobs"]
    df_new_jobs_nersc_iro2["source"] = "nersc"

    df_iro2_nersc = df_dict["df"]
    df_iro2_nersc["source"] = "nersc"

# #############################################################################
# Parsing Sherlock IrO2 DFT Data
path_i = os.path.join(
    os.environ["PROJ_DATA"],
    "04_IrOx_surfaces_OER/ml_bulk_irox_dft/iro2",
    "df_dict_sher.pickle")
with open(path_i, "rb") as fle:
    df_dict = pickle.load(fle)
    df_new_jobs_sher_iro2 = df_dict["df_new_jobs"]
    df_new_jobs_sher_iro2["source"] = "sherlock"

    df_iro2_sherlock = df_dict["df"]
    df_iro2_sherlock["source"] = "sherlock"

# #############################################################################
# Parsing SLAC IrO2 DFT Data
path_i = os.path.join(
    os.environ["PROJ_DATA"],
    "04_IrOx_surfaces_OER/ml_bulk_irox_dft/iro2",
    "df_dict_slac.pickle")
with open(path_i, "rb") as fle:
    df_dict = pickle.load(fle)

    df_new_jobs_slac_iro2 = df_dict["df_new_jobs"]
    df_new_jobs_slac_iro2["source"] = "slac"

    df_iro2_slac = df_dict["df"]
    df_iro2_slac["source"] = "slac"

# #############################################################################
df_new_jobs_iro2 = pd.concat([
    df_new_jobs_nersc_iro2,
    df_new_jobs_sher_iro2,
    df_new_jobs_slac_iro2,
    ])
df_new_jobs_iro2["stoich"] = "AB2"


# #############################################################################
df_iro2_long = pd.concat([
    df_iro2_nersc,
    df_iro2_sherlock,
    df_iro2_slac,
    ], axis=0, sort=True)

df_iro2_long["stoich"] = "AB2"

# +
# df_iro2_long.head()

df_iro2_long[df_iro2_long["id"] == "338"]

# +
data_list = []
grouped = df_iro2_long.groupby(["id"])
for name, group in grouped:

    # df_succ = group[group["job_state"] == "SUCCEEDED"]
    df_succ = group[group["completed"] == True]
    isif_2_done = 2 in df_succ["isif"].tolist()

    if len(df_succ) > 0 and isif_2_done:
        latest_succ_rev = df_succ.sort_values("revision").iloc[-1]
        data_list.append(latest_succ_rev)
df_iro2 = pd.DataFrame(data_list)

# Droping all unnecessary columns
df_iro2 = drop_columns(df=df_iro2, columns=["atoms", "path"], keep_or_drop="keep")

# Adding stoich column
df_iro2["stoich"] = "AB2"

# + active=""
# # #############################################################################
# # #############################################################################
# # #############################################################################
# # #############################################################################
# -

# # IrO3 Bulk Data

# +
# #############################################################################
# Parsing NERSC IrO3 DFT Data
path_i = os.path.join(
    os.environ["PROJ_DATA"],
    "04_IrOx_surfaces_OER/ml_bulk_irox_dft/iro3",
    "df_dict_nersc.pickle")
    # "df_dict.pickle")
with open(path_i, "rb") as fle:
    df_dict = pickle.load(fle)

    df_new_jobs_nersc_iro3 = df_dict["df_new_jobs"]
    df_new_jobs_nersc_iro3["source"] = "nersc"

    df_iro3_nersc = df_dict["df"]
    df_iro3_nersc["source"] = "nersc"

# #############################################################################
# Parsing Sherlock IrO3 DFT Data
path_i = os.path.join(
    os.environ["PROJ_DATA"],
    "04_IrOx_surfaces_OER/ml_bulk_irox_dft/iro3",
    "df_dict_sher.pickle")
with open(path_i, "rb") as fle:
    df_dict = pickle.load(fle)

    df_new_jobs_sher_iro3 = df_dict["df_new_jobs"]
    df_new_jobs_sher_iro3["source"] = "sherlock"

    df_iro3_sher = df_dict["df"]
    df_iro3_sher["source"] = "sherlock"

# #############################################################################
df_new_jobs_iro3 = pd.concat([
    df_new_jobs_nersc_iro3,
    df_new_jobs_sher_iro3,
    # df_new_jobs_slac_iro2,
    ])
df_new_jobs_iro3["stoich"] = "AB3"

# #############################################################################
df_iro3_long = pd.concat([
    df_iro3_nersc,
    df_iro3_sher,
    ], axis=0, sort=True)

df_iro3_long["stoich"] = "AB3"

# +
## Processing IrO3 Dataframe
data_list = []
grouped = df_iro3_long.groupby(["pre_path"])
for name, group in grouped:

    # df_succ = group[group["job_state"] == "SUCCEEDED"]
    df_succ = group[group["completed"] == True]

    isif_2_done = 2 in df_succ["isif"].tolist()

    if len(df_succ) > 0 and isif_2_done:
        latest_succ_rev = df_succ.sort_values("revision").iloc[-1]
        data_list.append(latest_succ_rev)
df_iro3 = pd.DataFrame(data_list)

# Droping all unnecessary columns
df_iro3 = drop_columns(df=df_iro3, columns=["atoms", "path"], keep_or_drop="keep")

# Adding stoich column
df_iro3["stoich"] = "AB3"

# + active=""
# # #############################################################################
# # #############################################################################
# # #############################################################################
# # #############################################################################
# -

# # Combining dataframes

# +
# #############################################################################
frames = [df_iro2, df_iro3]
df_m = pd.concat(frames)
print("df_m.shape:", df_m.shape)

# #############################################################################
frames = [df_new_jobs_iro2, df_new_jobs_iro3]
df_new_jobs = pd.concat(frames)
print("df_new_jobs.shape:", df_new_jobs.shape)

# #############################################################################
frames = [df_iro2_long, df_iro3_long]
df_irox_long = pd.concat(frames)
print("df_irox_long.shape:", df_irox_long.shape)


# -

# # Parsing ID from path

# +
# row_i = df_m.iloc[0]

def method(row_i):
    """
    """
    path_i = row_i["path"]
    folder_lists__isdigit = [i for i in path_i.split("/") if i.isdigit()]

    mess_i = "Must have only one folder in path that is numeric"
    assert len(folder_lists__isdigit) == 1, mess_i

    id_i =int(folder_lists__isdigit[0])

    return(id_i)

df_m["id_old"] = df_m.apply(
    method,
    axis=1)
# -

# # Writing Data

# +
with open(os.path.join(directory, "df_bulk_raul_iro2.pickle"), "wb") as fle:
    pickle.dump(df_iro2, fle)

with open(os.path.join(directory, "df_bulk_raul_iro3.pickle"), "wb") as fle:
    pickle.dump(df_iro3, fle)

with open(os.path.join(directory, "df_bulk_raul_irox.pickle"), "wb") as fle:
    pickle.dump(df_m, fle)

# #############################################################################
with open(os.path.join(directory, "df_new_jobs.pickle"), "wb") as fle:
    pickle.dump(df_new_jobs, fle)

with open(os.path.join(directory, "df_irox_long.pickle"), "wb") as fle:
    pickle.dump(df_irox_long, fle)
# -

print("AB2:", df_m[df_m["stoich"] == "AB2"].shape)
print("AB3:", df_m[df_m["stoich"] == "AB3"].shape)

assert False

# + {"active": ""}
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

# # Keeping track of competed jobs (Manually)

# + active=""
# # IrO2
# # ########################
#
# (368, 4)
# (442, 4)
# (445, 4)
# (446, 4)
# (448, 4)
# (450, 4)
# (454, 4)
# (460, 4)
# (463, 4)
# (468, 4)
# (470, 4)
# (475, 4)
# (478, 4)
# (480, 4)
# (482, 4)
# (484, 4)
# (485, 4)
# (487, 4)

# + active=""
# # IrO3
# # ########################
#
# (238, 4)
# (241, 4)
# (244, 4)
# (245, 4)
# (246, 4)
# (248, 4)
# (249, 4)

# + {"active": ""}
# # IrO2
# # ########################
#
# (99, 4)
# (100, 4)
# (102, 4)
# (106, 4)
# (108, 4)
# (131, 4)
# (227, 4)
# (230, 4)
# (231, 4)
# (233, 4)
# (235, 4)
# (236, 4)
# (237, 4)
# (239, 4)
# (242, 4)
# (244, 4)
# (245, 4)
# (247, 4)
# (250, 4)
# (253, 4)
# (254, 4)
# (256, 4)
# (258, 4)
# (266, 4)
# (272, 4)
# (273, 4)
# (276, 4)
# (279, 4)
# (283, 4)
# (285, 4)
# (287, 4)
# (289, 4)
# (292, 4)
# (295, 4)
# (297, 4)
# (305, 4)
# (310, 4)
# (313, 4)
# (317, 4)
# (322, 4)
# (326, 4)
# (330, 4)
# (335, 4)
# (340, 4)
# (354, 4)
# (357, 4)
# (359, 4)
# (367, 4)
#
# # IrO3
# # ########################
#
# (37, 4)
# (47, 4)
# (55, 4)
# (60, 4)
# (62, 4)
# (87, 4)
# (116, 4)
# (119, 4)
# (120, 4)
# (121, 4)
# (123, 4)
# (124, 4)
# (127, 4)
# (128, 4)
# (132, 4)
# (135, 4)
# (146, 4)
# (150, 4)
# (155, 4)
# (158, 4)
# (161, 4)
# (177, 4)
# (188, 4)

# + active=""
#
#
#

# +
df_iro2_slac
pd.set_option('display.max_rows', None)
# 13 | errored
# 58 | finished
df_tmp = df_new_jobs_slac_iro2

df_tmp = df_tmp[df_tmp["action"] != "ALL DONE! | ISIF 2"]

df_tmp = df_tmp[df_tmp["action"] != "Job is busy, will skip"]
print("df_tmp.shape:", df_tmp.shape)
display(df_tmp)

print("All jobs:", df_new_jobs_slac_iro2.shape[0], "\n")
print("Busy jobs:", df_new_jobs_slac_iro2[df_new_jobs_slac_iro2["action"] == "Job is busy, will skip"].shape[0])
print("Finished jobs:", df_new_jobs_slac_iro2[df_new_jobs_slac_iro2["action"] == "ALL DONE! | ISIF 2"].shape[0])
# -

# df_iro2_slac
#
# df_new_jobs_slac_iro2[df_new_jobs_slac_iro2["id"] == "590"]

# +
[i for i in df_new_jobs_sher_iro2["pre_path"].tolist() if "iro3" in i]

[i for i in df_iro2_sherlock["path"].tolist() if "iro3" in i]

# +
# df_m[df_m["stoich"] == "AB2"]
# -

df_m[df_m["id_old"] == 182]
