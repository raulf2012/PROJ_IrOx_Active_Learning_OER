# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.2.1
#   kernelspec:
#     display_name: Python [conda env:PROJ_IrOx_Active_Learning_OER]
#     language: python
#     name: conda-env-PROJ_IrOx_Active_Learning_OER-py
# ---

# # Import Modules

# +
#| - Import Modules
import os
import pickle
import pandas as pd
import shutil

from methods import parse_job_err
from methods import parse_finished_file
from methods import parse_job_state
from methods import is_job_submitted
from methods import get_isif_from_incar
from methods import read_write_CONTCAR

from methods import set_up__submit__new_job
#__|
# -

# # Script Inputs

# +
cwd = os.getcwd()

# Set to True when you're actually ready
change_file_sys = False

root_dir = "__test__/job_folders"
# -

# # Coping over job folders to easily recover state

# +
shutil.rmtree("./__test__/job_folders")

shutil.copytree(
    "./__test__/job_folders.orig",
    "./__test__/job_folders",
    )

# + {"jupyter": {"source_hidden": true}}
data_list = []

job_dirs = []


for subdir, dirs, files in os.walk(root_dir):
    last_dir = subdir.split("/")[-1]
    cond_0 = last_dir[0] == "_"
    cond_1 = last_dir[1:].isdigit()
    if cond_0 and cond_1:
        revision_i = int(last_dir[1:])

        job_pre_path_i = "/".join(subdir.split("/")[:-1])

        out_dict = dict(
            path=subdir,
            pre_path=job_pre_path_i,
            revision=revision_i,
            )
        data_list.append(out_dict)

        job_dirs.append(subdir)

df = pd.DataFrame(data_list)


# -

# # Parsing dirs to get job state info

# + {"jupyter": {"source_hidden": true}}
#| - Parsing dirs to get job state info
# #############################################################################
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
#__|
# -

# # Save dataframe to file and rclone to Dropbox

# +
# directory = "out_data"
# if not os.path.exists(directory):
#     os.makedirs(directory)

# with open("out_data/df.pickle", "wb") as fle:
#     pickle.dump(df, fle)


# if os.environ["USER"] == "flores12":
#     print("On NERSC probably")

#     db_path = os.path.join(
#         "01_norskov/00_git_repos/PROJ_IrOx_Active_Learning_OER",
#         "run_nersc_vasp/ml_bulk_opt/run_all_bulks/out_data/")

#     os.system("rclone copy out_data/df.pickle raul_dropbox:" + db_path)
# -


df

# # TEST ---------------

# +
unique_pre_paths = df["pre_path"].unique()

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
#     print(latest_rev["path"])




    # #########################################################################
    # #########################################################################
    # #########################################################################





    path = latest_rev["path"]
    job_state = latest_rev["job_state"]
    timed_out = latest_rev["timed_out"]
    isif = latest_rev["isif"]
    completed = latest_rev["completed"]
    pre_path = latest_rev["pre_path"]


    data_dict_i["pre_path"] = pre_path


    failed = False
    skip_job = False

    new_job_file_dict = dict()

    cond_0 = job_state == "RUNNING"
    cond_1 = job_state == "PENDING"
    cond_2 = job_state == "CONFIGURING"
    if cond_0 or cond_1 or cond_2:
        #| - If job is either running, pending or being configured
        # (whatever that means), just move on for now
        mess_i = "Job is busy, will skip"
        data_dict_i["action"] = mess_i

        skip_job = True
        pass
        #__|

    elif job_state == "SUCCEEDED" or completed:
        #| - SUCCEEDED
        # Picking the model.py script to use
        mess_i = "Job succeeded"
        data_dict_i["action"] = mess_i

        read_write_CONTCAR(path, new_job_file_dict)

        if isif == 7:
            model_file_path = os.path.join(
                os.environ["PROJ_irox"],
                "run_nersc_vasp/ml_bulk_opt",
                "bulk_opt_last.py")
            new_job_file_dict[model_file_path] = "model.py"

            data_dict_i["action"] += " | ISIF 7 calc finished, moving to ISIF 3"

        elif isif == 3:
            data_dict_i["action"] += " | ISIF 3 calc finished, moving to "

            if num_completed_isif_3 < 3:
                model_file_path = os.path.join(
                    os.environ["PROJ_irox"],
                    "run_nersc_vasp/ml_bulk_opt",
                    "bulk_opt_last.py")
                new_job_file_dict[model_file_path] = "model.py"
                data_dict_i["action"] += " | Running isif 3 again"
            else:
                model_file_path = os.path.join(
                    os.environ["PROJ_irox"],
                    "run_nersc_vasp/ml_bulk_opt",
                    "bulk_opt_init_final_isif_2.py")
                new_job_file_dict[model_file_path] = "model.py"
                data_dict_i["action"] += " | Moving to isif 2"


        elif isif == 2:
            data_dict_i["action"] = "Final ISIF 2 finished!"

            pass
        #__|

    elif timed_out or failed:
        #| - Timed out of failed
        print("timed out or failed")

        data_dict_i["action"] = "Time our or failed"

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
        set_up__submit__new_job(latest_rev, new_job_file_dict, )

    data_list.append(data_dict_i)




df_new_jobs = pd.DataFrame(data_list)
df_new_jobs
# -

data_list

# + {"active": ""}
# # Job rules
#
# ## If the job is running or pending or configuring, ignore for now
# ## If the job timed_out or has a filed job_state, then resubmit
# ## If the job has succeeded job_state, then check if isif is 7 or 3,
#   ## if 7 resubmit with 3, if 3 then resubmit 1 more time with 3, then finally submit one more time with isif 2
#

# + {"active": ""}
#
#
#
#
#
# + {}
directory = "out_data"
if not os.path.exists(directory):
    os.makedirs(directory)

df_dict = {
    "df": df,
    "df_new_jobs": df_new_jobs,
    }

with open("out_data/df_dict.pickle", "wb") as fle:
    pickle.dump(df, fle)


if os.environ["USER"] == "flores12":
    print("On NERSC probably")

    db_path = os.path.join(
        "01_norskov/00_git_repos/PROJ_IrOx_Active_Learning_OER",
        "run_nersc_vasp/ml_bulk_opt/run_all_bulks/out_data/")

    os.system("rclone copy out_data/df_dict.pickle raul_dropbox:" + db_path)
# -


df

# + {"jupyter": {"source_hidden": true}}
# path


# contcar_file_path = os.path.join(path, "CONTCAR")
# my_file = Path(contcar)
# if my_file.is_file():
#     atoms_for_next_job = io.read(contcar_file_path)
