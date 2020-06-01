# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Import Modules

# + jupyter={}
import os
import sys

import json
import pickle
from shutil import copyfile

from pathlib import Path

# #############################################################################
import pandas as pd

from ase import io
from ase.visualize import view

# #############################################################################
from dft_job_automat.compute_env import ComputerCluster

from methods import calc_wall_time


sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import proj_dir_name

# +
model_file = os.path.join(
    os.environ["PROJ_irox"],
    "run_nersc_vasp/ml_bulk_opt",
    "bulk_opt_init.py",
    )

rootdir = os.getcwd()

# +
structure_files_dir = os.path.join(
    os.environ["PROJ_DATA"],
    "04_IrOx_surfaces_OER/ml_bulk_static_preopt_structures",
    "fixed_prototypes_iro2",
    )

data_list = []
for subdir, dirs, files in os.walk(structure_files_dir):
    for file in files:
        file_path_i = os.path.join(subdir, file)

        path_short = file_path_i.replace(os.environ["PROJ_irox"] + "/", "")

        atoms_i = io.read(file_path_i)

        num_atoms_i = atoms_i.get_number_of_atoms()

        data_list.append({
            "path": file_path_i,
            "path_short": path_short,
            "file_name": file,
            "init_atoms": atoms_i,
            "number_of_atoms": num_atoms_i,
            "id": int(file.split("_")[0]),
            })

df = pd.DataFrame(data_list)
df = df.sort_values("id")

# Filtering large structures (more than 100 atoms in cell)
df = df[df["number_of_atoms"] < 100]
# -

# # Filtering with id list

from ids_to_run import ids_to_run
df = df[df["id"].isin(ids_to_run)]

directory = rootdir
if not os.path.exists(directory):
    os.makedirs(directory)


def submit_job(
    wall_time_i=None,
    nodes_i=None,
    job_0_dir_i=None,
    ):
    CC = ComputerCluster()

    if os.environ["COMPENV"] == "sherlock":
        def_params = {
            "wall_time": wall_time_i,
            "nodes": nodes_i,
            "path_i": job_0_dir_i}

    elif os.environ["COMPENV"] == "slac":
        def_params = {
            "wall_time": wall_time_i,
            "cpus": 12,
            "queue": "suncat2",
            "path_i": job_0_dir_i}

    else:
        def_params = {
            "wall_time": wall_time_i,
            # "queue": "premium",
            "queue": "regular",
            "architecture": "knl",
            "nodes": nodes_i,
            "path_i": job_0_dir_i}

    CC.submit_job(**def_params)


for i_cnt, row_i in df.iterrows():
    folder_name = str(row_i["id"]).zfill(3)

    job_root_dir_i = os.path.join(
        rootdir,
        "job_folders",
        folder_name)

    job_0_dir_i = os.path.join(job_root_dir_i, "_1")


    # Checkif job dir is already present
    job_folder = Path(job_root_dir_i)
    if job_folder.is_dir():
        print("Job dir already exists, skipping")

    else:
        try:
            os.makedirs(job_0_dir_i)
        except:
            pass

        # Write atoms object
        atoms_i = row_i["init_atoms"]
        io.write(os.path.join(job_0_dir_i, "init.cif"), row_i["init_atoms"])


        dft_params_dict = {
            # "encut": 600,
            # "kpar": 5,
            # "ediffg": 5e-3,
            # "ediff": 1e-6
            }



        num_atoms = atoms_i.get_number_of_atoms()
        wall_time_i = calc_wall_time(num_atoms, factor=1.4)
        wall_time_i = int(wall_time_i)

        if num_atoms > 100:
            nodes_i = 10
            dft_params_dict["kpar"] = 10
        else:
            nodes_i = 5
            dft_params_dict["kpar"] = 5

        if os.environ["COMPENV"] == "sherlock":
            print("SDIJFIDSJIFJDISJFIJSDIFJIDSJF")
            dft_params_dict["npar"] = 4

        if os.environ["COMPENV"] == "slac":
            dft_params_dict["kpar"] = 3
            dft_params_dict["npar"] = 4

        if os.environ["COMPENV"] != "slac":
            if wall_time_i > 600:
                wall_time_i = 600
        else:
            wall_time_i = 8. * wall_time_i

            if wall_time_i > 2760:
                wall_time_i = 2760

        
        # Write dft paramters json file to job dir
        with open(os.path.join(job_0_dir_i, "dft-params.json"), "w+") as fle:
            json.dump(dft_params_dict, fle, indent=2, skipkeys=True)

        # Copy model job script
        copyfile(model_file, os.path.join(job_0_dir_i, "model.py"))


        # Submit job ##############################################################
        submit_job(
            wall_time_i=int(wall_time_i),
            nodes_i=nodes_i,
            job_0_dir_i=job_0_dir_i)

# + active=""
#
#
#
# -

# # OLD

# + jupyter={}
#         CC = ComputerCluster()

#         if os.environ["COMPENV"] == "sherlock":
#             def_params = {
#                 "wall_time": wall_time_i,
#                 "nodes": nodes_i,
#                 "path_i": job_0_dir_i}
#         elif os.environ["COMPENV"] == "slac":
#             def_params = {
#                 "wall_time": wall_time_i,
#                 "cpus": 30,
#                 "path_i": job_0_dir_i}

#         else:
#             def_params = {
#                 "wall_time": wall_time_i,
#                 # "queue": "premium",
#                 "queue": "regular",
#                 "architecture": "knl",
#                 "nodes": nodes_i,
#                 "path_i": job_0_dir_i}

#         CC.submit_job(**def_params)
