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
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Import Modules

# +
# | - Import Modules
import os
import pickle

import pandas as pd
pd.set_option('display.max_rows', None)

from ase import io

import shutil

from methods import parse_job_err
from methods import parse_finished_file
from methods import parse_job_state
from methods import is_job_submitted
from methods import get_isif_from_incar
from methods import read_write_CONTCAR
from methods import get_number_of_ionic_steps


from methods import set_up__submit__new_job
#__|
# -

# # Script Inputs

try:
    from inputs import input_dict
#     root_dir = input_dict["root_dir"]
    root_dir = input_dict.get("root_dir", ".")
    ignore_ids = input_dict.get("ignore_ids", [])
except:
    print("Couldn't import inputs script")
    input_dict = {}
    root_dir = "."
    ignore_ids = []

# +
cwd = os.getcwd()

# Set to True when you're actually ready
change_file_sys = False

# root_dir = "__test__/job_folders_2"
# root_dir = "."
# root_dir = "__test__/iro2_calcs_test"

# + {"active": ""}
# # Job rules
#
# ## If the job is running or pending or configuring, ignore for now
# ## If the job timed_out or has a filed job_state, then resubmit
# ## If the job has succeeded job_state, then check if isif is 7 or 3,
#   ## if 7 resubmit with 3, if 3 then resubmit 1 more time with 3, then finally submit one more time with isif 2
# -

# # Parsing dir to find job folder

# +
done_pre_paths = []
if False:
    path_i = "out_data/df_dict.pickle"
    with open(path_i, "rb") as fle:
        data = pickle.load(fle)
        df_tmp = data["df"]
        done_pre_paths = df_tmp["pre_path"].tolist()

# except:
#     print("THIS FAILED!!")
#     done_pre_paths = []
# -

print("done_pre_paths:", done_pre_paths)

# +
data_list = []
job_dirs = []
for subdir, dirs, files in os.walk(root_dir):
    print("TEMP: ", subdir)

    if subdir in done_pre_paths:
        print("skipping this dir, already done")
        continue

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
# -

# # Reducing jobs to parse

# +
# # path_i = "/home/raulf2012/Dropbox/01_norskov/PROJECT_DATA/04_IrOx_surfaces_OER/ml_bulk_irox_dft/iro3/df_dict.pickle"

# try:
#     path_i = "out_data/df_dict.pickle"
#     with open(path_i, "rb") as fle:
#         data = pickle.load(fle)

#     df_new_jobs = data["df_new_jobs"]
# #     df_tmp = data["df"]

#     done_str = "ALL DONE! | ISIF 2"
#     completed_jobs_pre_path_list = df_new_jobs[
#         df_new_jobs["action"] == done_str]["pre_path"].tolist()

#     df = df[~df["pre_path"].isin(completed_jobs_pre_path_list)]
# except:
#     print("Couldn't reduce `pre_path` to save time")
# -

# # Parsing dirs to get job state info

# +
# | - Parsing dirs to get job state info
# #############################################################################
print("Parsing dirs to get job state info")
print(80 * "&")
print(80 * "&")
print(80 * "&")
print(80 * "&")

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

# +
import os
import pickle

directory = "out_data"
if not os.path.exists(directory):
    os.makedirs(directory)

with open(os.path.join(directory, "df.pickle"), "wb") as fle:
    pickle.dump(df, fle)
# -

# # TEST ---------------

# +
print(80 * "*"); print(80 * "*")

# unique_pre_paths = df["pre_path"].unique()

data_list = []
grouped = df.groupby(["pre_path"])
for name, group in grouped:
#     group = df[df["pre_path"] == "__test__/job_folders/000"]

    data_dict_i = {}

    group = group.sort_values("revision", ascending=False)


    # Number of completed isif = 3 jobs
    # Should be run maybe 2-3 times to fully converge
    df_tmp = group[group["isif"] == 3]
    num_completed_isif_3 = len(df_tmp[df_tmp["completed"] == True])

    data_dict_i["num_completed_isif_3"] = num_completed_isif_3

    latest_rev = group.iloc[0]

    data_dict_i["num_revisions"] = latest_rev["revision"]

    print("Main", latest_rev["path"])




    # #########################################################################
    # #########################################################################
    # #########################################################################





    path = latest_rev["path"]
    job_state = latest_rev["job_state"]
    timed_out = latest_rev["timed_out"]
    isif = latest_rev["isif"]
    completed = latest_rev["completed"]
    pre_path = latest_rev["pre_path"]
    id_i = latest_rev["id"]

    failed = latest_rev["error"]
    error_type = latest_rev["error_type"]

    data_dict_i["pre_path"] = pre_path
    data_dict_i["id"] = id_i

    skip_job = False

    new_job_file_dict = dict()

    # Ignoring manually selected jobs
    if id_i in ignore_ids:
        data_dict_i["action"] = "Ignoring this id"
        skip_job = True

    else:

        cond_0 = job_state == "RUNNING"
        cond_1 = job_state == "PENDING"
        cond_2 = job_state == "CONFIGURING"
        if cond_0 or cond_1 or cond_2:
            # | - If job is either running, pending or being configured
            # (whatever that means), just move on for now
            mess_i = "Job is busy, will skip"
            data_dict_i["action"] = mess_i

            skip_job = True
            pass
            #__|

        elif job_state == "SUCCEEDED" or completed:
            # | - SUCCEEDED
            # Picking the model.py script to use
            mess_i = "Job done"
            data_dict_i["action"] = mess_i

            read_write_CONTCAR(path, new_job_file_dict)

            if isif == 7:
                model_file_path = os.path.join(
                    os.environ["PROJ_irox"],
                    "run_nersc_vasp/ml_bulk_opt",
                    "bulk_opt_last.py")
                new_job_file_dict[model_file_path] = "model.py"

                data_dict_i["action"] += " | ISIF 7 done, --> ISIF 3"

            elif isif == 3:
                data_dict_i["action"] += " | ISIF 3 done"

                if num_completed_isif_3 < 3:
                    model_file_path = os.path.join(
                        os.environ["PROJ_irox"],
                        "run_nersc_vasp/ml_bulk_opt",
                        "bulk_opt_last.py")
                    new_job_file_dict[model_file_path] = "model.py"
                    data_dict_i["action"] += " | Rerunning isif 3"
                else:
                    model_file_path = os.path.join(
                        os.environ["PROJ_irox"],
                        "run_nersc_vasp/ml_bulk_opt",
                        "bulk_opt_init_final_isif_2.py")
                    new_job_file_dict[model_file_path] = "model.py"
                    data_dict_i["action"] += " | --> isif 2"


            elif isif == 2:
                skip_job = True
                data_dict_i["action"] = "ALL DONE! | ISIF 2"

                pass
            #__|

        elif timed_out or failed:
            # | - Timed out of failed
            print("timed out or failed")

            if error_type == "Error in SGRCON (symm error)":
                data_dict_i["action"] = "Error, need manual attention"
                skip_job = True
            else:
                data_dict_i["action"] = "Time out or failed"

                read_write_CONTCAR(path, new_job_file_dict)


                # Picking the model.py script to use
                if isif == 7:
                    model_file_path = os.path.join(
                        os.environ["PROJ_irox"],
                        "run_nersc_vasp/ml_bulk_opt",
                        "bulk_opt_init.py")
                    new_job_file_dict[model_file_path] = "model.py"
                    data_dict_i["action"] += " | Restarting isif 7 calc"

                elif isif == 3:
                    model_file_path = os.path.join(
                        os.environ["PROJ_irox"],
                        "run_nersc_vasp/ml_bulk_opt",
                        "bulk_opt_last.py")
                    new_job_file_dict[model_file_path] = "model.py"
                    data_dict_i["action"] += " | Restarting isif 3 calc"

                elif isif == 2:
                    model_file_path = os.path.join(
                        os.environ["PROJ_irox"],
                        "run_nersc_vasp/ml_bulk_opt",
                        "bulk_opt_init_final_isif_2.py")
                    new_job_file_dict[model_file_path] = "model.py"
                    data_dict_i["action"] += " | Restarting isif 2 calc"
            #__|

        else:
            skip_job = True

            mess_i = "Couldn't figure out what to do"
            data_dict_i["action"] = mess_i
            print(mess_i)


    if not skip_job and change_file_sys:
        set_up__submit__new_job(
            latest_rev,
            new_job_file_dict,
            run_calc=True)

    data_list.append(data_dict_i)

df_new_jobs = pd.DataFrame(data_list)

# +
df_not_done = df_new_jobs[df_new_jobs["action"] != "ALL DONE! | ISIF 2"]

df_timed_out = df_new_jobs[df_new_jobs["action"] == "Time out or failed | Restarting isif 3 calc"]

df_not_sure = df_new_jobs[df_new_jobs["action"] == "Couldn't figure out what to do"]
# -

print("")
print(80 * "#")
print("df_not_sure | df_timed_out | df_not_done | df_new_jobs")
print(80 * "#")

# # Saving data and uploading to Dropbox

# +
print("")

directory = "out_data"
if not os.path.exists(directory):
    os.makedirs(directory)

df_dict = {
    "df": df,
    "df_new_jobs": df_new_jobs,
    }

import datetime
datetime_i = datetime.datetime.now().strftime("%y%m%d__%H_%M")

filename_i = "out_data/" + "df_dict.pickle"
with open(filename_i, "wb") as fle:
    pickle.dump(df_dict, fle)

filename_i = "out_data/" + datetime_i + "_df_dict.pickle"
with open(filename_i, "wb") as fle:
    pickle.dump(df_dict, fle)



if os.environ["USER"] == "flores12":

    db_path_bkp = os.path.join(
        "01_norskov/00_git_repos/PROJ_IrOx_Active_Learning_OER",
        "run_nersc_vasp/ml_bulk_opt/run_all_bulks/out_data/")

    db_path = input_dict.get(
        "proj_data_save_dir",
        db_path_bkp)

    out_file_name = input_dict.get(
        "out_file_name",
        "df_dict.pickle")

    db_path = os.path.join(db_path, out_file_name)

    
    # bash_comm = "rclone copy out_data/df_dict.pickle raul_dropbox:" + db_path
    bash_comm = "rclone copyto " + filename_i + " raul_dropbox:" + db_path
    print(bash_comm)
    os.system(bash_comm)
# -
# Deleting temp files
os.system("rm init.cif run_vasp.py")

# + {"active": ""}
#
#
#

# + {"jupyter": {"source_hidden": true}}
# Coping over job folders to easily recover state

# shutil.rmtree("./__test__/job_folders")

# shutil.copytree(
#     "./__test__/job_folders.orig",
#     "./__test__/job_folders",
#     )

# df_new_jobs

# # %load_ext autoreload
# # %autoreload 2
# from methods import parse_job_err

# row_i = df.iloc[0]
# path_i = row_i["path"]
# path_i



# parse_job_err(path_i)
