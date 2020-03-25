# | - Import Modules
import os
import sys

import pandas as pd

from ase import io

sys.path.insert(0, "..")
from methods import (
    parse_job_err,
    parse_finished_file,
    parse_job_state,
    is_job_submitted,
    get_isif_from_incar,
    get_number_of_ionic_steps,
    )
#__|

# | - Script Inputs
paths_to_parse_list = [
    "/scratch/users/flores12/PROJ_irox_ml_oer/beta_iro3_phase_calc/01_attempt",
    ]
#__|


# | - Parsing files
for path_i in paths_to_parse_list:
    print(path_i)
    print("DIFJIDJFISJDIFJISDJIFJSDi")
    # parse_job_err(path_i)


    data_list = []
    job_dirs = []
    for subdir, dirs, files in os.walk(path_i):
        print("TEMP: ", subdir)

        # if subdir in done_pre_paths:
        #     print("skipping this dir, already done")
        #     continue

        last_dir = subdir.split("/")[-1]
        cond_0 = last_dir[0] == "_"
        cond_1 = last_dir[1:].isdigit()
        if cond_0 and cond_1:
            print("Parsing dirs:", subdir)

            revision_i = int(last_dir[1:])

            job_pre_path_i = "/".join(subdir.split("/")[:-1])
            id_i = job_pre_path_i.split("/")[-1]

            try:
                atoms_path_i = os.path.join(subdir, "OUTCAR")
                atoms_i = io.read(atoms_path_i)
            except:
                atoms_i = None

            out_dict = dict(
                path=subdir,
                pre_path=job_pre_path_i,
                id=id_i,
                revision=revision_i,
                atoms=atoms_i)
            data_list.append(out_dict)

            job_dirs.append(subdir)

    df = pd.DataFrame(data_list)
#__|

# | - Apply job parse methods
def method(row_i):
    status_dict = parse_job_err(row_i["path"])
    for key, value in status_dict.items():
        row_i[key] = value
    return(row_i)
df = df.apply(method, axis=1)

# #############################################################################
def method(row_i):
    status_dict = parse_finished_file(row_i["path"])
    for key, value in status_dict.items():
        row_i[key] = value
    return(row_i)
df = df.apply(method, axis=1)

# #############################################################################
def method(row_i):
    status_dict = parse_job_state(row_i["path"])
    for key, value in status_dict.items():
        row_i[key] = value
    return(row_i)
df = df.apply(method, axis=1)

# #############################################################################
def method(row_i):
    status_dict = is_job_submitted(row_i["path"])
    for key, value in status_dict.items():
        row_i[key] = value
    return(row_i)
df = df.apply(method, axis=1)

# #############################################################################
def method(row_i):
    status_dict = get_isif_from_incar(row_i["path"])
    for key, value in status_dict.items():
        row_i[key] = value
    return(row_i)
df = df.apply(method, axis=1)

# #############################################################################
def method(row_i):
    status_dict = get_number_of_ionic_steps(row_i["path"])
    for key, value in status_dict.items():
        row_i[key] = value
    return(row_i)
df = df.apply(method, axis=1)
#__|



# | - Write and upload dataframe to Dropbox


# Write pickle file
# Pickling data ###############################################################
import os; import pickle
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)
with open(os.path.join(directory, "df_manual.pickle"), "wb") as fle:
    pickle.dump(object, fle)
# #############################################################################

# db_path = input_dict.get(
#     "proj_data_save_dir",
#     db_path_bkp)

proj_data_save_dir = os.path.join(
    "01_norskov/PROJECT_DATA",
    "04_IrOx_surfaces_OER/ml_bulk_irox_dft/manual_calcs/iro3")


# out_file_name = input_dict.get(
#     "out_file_name",
#     "df_dict.pickle")

out_file_name = "df_manual.pickle"

db_path = os.path.join(
    proj_data_save_dir,
    # out_file_name
    )


bash_comm = "rclone copy out_data/df_manual.pickle raul_dropbox:" + db_path
# bash_comm = "rclone copyto " + filename_i + " raul_dropbox:" + db_path
print(bash_comm)
os.system(bash_comm)
#__|
