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
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import os
import sys

import json
from shutil import copyfile

# #############################################################################
import pandas as pd
from ase import io
from ase.visualize import view

# #############################################################################
from dft_job_automat.compute_env import ComputerCluster

# +
# model_file = "./model_0.py"

model_file = os.path.join(
    os.environ["PROJ_irox"],
    "run_nersc_vasp/ml_bulk_opt",
    "bulk_opt_init.py",
    )
# -

rootdir = os.getcwd()

# +
structure_files_dir = os.path.join(
    os.environ["PROJ_irox"],
    "chris_prototypes_structures/fixed_prototypes_iro2")

data_list = []
for subdir, dirs, files in os.walk(structure_files_dir):
    for file in files:
        file_path_i = os.path.join(subdir, file)

        atoms_i = io.read(file_path_i)

        num_atoms_i = atoms_i.get_number_of_atoms()

        data_list.append({
            "path": file_path_i,
            "file_name": file,
            "init_atoms": atoms_i,
            "number_of_atoms": num_atoms_i,
            "id": int(file.split("_")[0]),
            })

df = pd.DataFrame(data_list)
df = df.sort_values("id")

# Filtering large structures (more than 100 atoms in cell)
df = df[df["number_of_atoms"] < 100]

# +
# TEMP
# df = df[0:3]
# -

# ## Methods

from methods import calc_wall_time

# +
# row_i = df.iloc[0]

for i_cnt, row_i in df.iterrows():
    folder_name = str(row_i["id"]).zfill(3)

    job_root_dir_i = os.path.join(
        rootdir,
        "job_folders",
        folder_name)

    job_0_dir_i = os.path.join(job_root_dir_i, "_1")

    try:
        os.makedirs(job_0_dir_i)
    except:
        pass

    
    dft_params_dict = {
        "encut": 600,
#         "kpar": 5,
        "ediffg": 5e-3,
        "ediff": 1e-6
        }

    # Write atoms object
    atoms_i = row_i["init_atoms"]
    io.write(os.path.join(job_0_dir_i, "init.cif"), row_i["init_atoms"])

    num_atoms = atoms_i.get_number_of_atoms()
    wall_time_i = calc_wall_time(num_atoms, factor=1.4)
    wall_time_i = int(wall_time_i)

    if num_atoms > 100:
        nodes_i = 10
        dft_params_dict["kpar"] = 10
    else:
        nodes_i = 5
        dft_params_dict["kpar"] = 5


    # Write dft paramters json file to job dir
    with open(os.path.join(job_0_dir_i, "dft-params.json"), "w+") as fle:
        json.dump(
            dft_params_dict, fle,
            indent=2, skipkeys=True)

    # Copy model job script
    copyfile(model_file, os.path.join(job_0_dir_i, "model.py"))
    

    # Submit job ##############################################################
    CC = ComputerCluster()

    def_params = {
        "wall_time": wall_time_i,
        "queue": "premium",
        "architecture": "knl",
        "nodes": nodes_i,
        "path_i": job_0_dir_i}

#     CC.submit_job(**def_params)

# +
# assert False

# +
# os.system("rm -r job_folders")

# + {"jupyter": {"source_hidden": true}}
# def calc_wall_time(num_atoms, factor=1.0):
#     """
#     """
#     term_2 = 0.0091 * (num_atoms ** 2)
#     term_1 = 2.2415 * (num_atoms ** 1)
#     term_0 = 29.760 * (num_atoms ** 0)
#     wall_time = factor * (term_2 + term_1 + term_0)
    
#     return(wall_time)


# wall_time = calc_wall_time(num_atoms, factor=1.4)

# row_i = df.iloc[0]

# atoms_i = row_i["init_atoms"]
# num_atoms = atoms_i.get_number_of_atoms()

# atoms_i

# with open(os.path.join(job_0_dir_i, "dft-params.json"), "w+") as fle:
#     json.dump(
#         dft_params_dict, fle,
#         indent=2, skipkeys=True)

# +
# df_tmp = df.sort_values("number_of_atoms", ascending=False)

# view(df_tmp.iloc[2]["init_atoms"])

# +
import chart_studio.plotly as py
import plotly.graph_objs as go
import os

trace = go.Scatter(
#     x=x_array,
    y=df["number_of_atoms"],
    mode="markers",
    )

data = [trace]

fig = go.Figure(data=data)
fig.show()

# + {"active": ""}
#
#
#

# + {"jupyter": {"source_hidden": true}}
# num_atoms_list = []
# for i_cnt, row_i in df.iterrows():
#     atoms_i = row_i["init_atoms"]
    
#     num_atoms_i = atoms_i.get_number_of_atoms()
    
#     num_atoms_list.append(num_atoms_i)

# def method(row_i):
#     """
#     """
#     atoms_i = row_i["init_atoms"]
#     num_atoms_i = atoms_i.get_number_of_atoms()

#     return(num_atoms_i)

# df["number_of_atoms"] = df.apply(
#     method,
#     axis=1,
#     )

# + {"jupyter": {"source_hidden": true}}
# def calc_wall_time(num_atoms, factor=1.0):
#     """
#     """
#     term_2 = 0.0091 * (num_atoms ** 2)
#     term_1 = 2.2415 * (num_atoms ** 1)
#     term_0 = 29.760 * (num_atoms ** 0)
#     wall_time = factor * (term_2 + term_1 + term_0)
    
#     return(wall_time)
