#!/usr/bin/env python

"""TEMP.

Author(s): Raul A. Flores
"""

#| - IMPORT MODULES
import os
# import sys
# import pickle

# import pandas as pd
from ase import io

from pathlib import Path

# #############################################################################
from dft_job_automat.compute_env import ComputerCluster

from vasp.vasp_methods import parse_incar







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
#__|

#| - Job parsing methods
# #############################################################################
def parse_job_err(path):
    status_dict = {
        "timed_out": None,
        }

    job_err_file_path = os.path.join(path, "job.err")
    my_file = Path(job_err_file_path)
    if my_file.is_file():
        with open(job_err_file_path, 'r') as f:
            lines = f.readlines()

        for line in lines:
            if "DUE TO TIME LIMIT" in line:
                status_dict["timed_out"] = True

    return(status_dict)

# #############################################################################
def parse_finished_file(path):
    status_dict = {
        "completed": None,
        }

    job_err_file_path = os.path.join(path, ".FINISHED")
    my_file = Path(job_err_file_path)
    if my_file.is_file():
        status_dict["completed"] = True

    return(status_dict)

# #############################################################################
def parse_job_state(path):
    print("")
    print(path)
    CC = ComputerCluster()
    job_state = CC.cluster.job_state(path_i=path)

    return({"job_state": job_state})

# #############################################################################
def is_job_submitted(path):
    status_dict = {
        "submitted": None,
        }

    submitted_file_path = os.path.join(path, ".SUBMITTED")
    my_file = Path(submitted_file_path)
    if my_file.is_file():
        status_dict["submitted"] = True
    return(status_dict)

# #############################################################################
def get_isif_from_incar(path):
    isif = None

    incar_file_path = os.path.join(path, "INCAR")
    my_file = Path(incar_file_path)
    if my_file.is_file():
        with open(incar_file_path, 'r') as f:
            lines = f.readlines()

            isif = parse_incar(lines)["ISIF"]
    return({"isif": isif})

#__|

def calc_wall_time(num_atoms, factor=1.0):
    """
    """
    term_2 = 0.0091 * (num_atoms ** 2)
    term_1 = 2.2415 * (num_atoms ** 1)
    term_0 = 29.760 * (num_atoms ** 0)
    wall_time = factor * (term_2 + term_1 + term_0)

    return(wall_time)

def set_up__submit__new_job(
    latest_revision,
    new_job_file_dict,
    run_calc=False,
    ):
    """
    """
    #| - set_up__submit__new_job
    path = latest_revision["path"]
    pre_path = latest_revision["pre_path"]
    current_revision = latest_revision["revision"]

    job_state = latest_revision["job_state"]
    timed_out = latest_revision["timed_out"]
    isif = latest_revision["isif"]


    new_revision = current_revision + 1
    job_0_dir_i = os.path.join(pre_path, "_" + str(new_revision))


    try:
        os.makedirs(job_0_dir_i)
    except:
        pass



    #| - Copy files into new directory
    for file_path_i, dist_name in new_job_file_dict.items():
        copyfile(file_path_i, os.path.join(job_0_dir_i, dist_name))
    #__|


    dft_params_dict = {
        # "encut": 600,
        # "kpar": 5,
        # "ediffg": 5e-3,
        # "ediff": 1e-6
        }

    # Write atoms object
    # atoms_i = row_i["init_atoms"]
    # io.write(os.path.join(job_0_dir_i, "init.cif"), row_i["init_atoms"])

    atoms_i = io.read(os.path.join(job_0_dir_i, "init.cif"))

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

    # # Copy model job script
    # copyfile(model_file, os.path.join(job_0_dir_i, "model.py"))


    # Submit job ##############################################################
    CC = ComputerCluster()

    def_params = {
        "wall_time": wall_time_i,
        # "queue": "premium",
        "queue": "regular",
        "architecture": "knl",
        "nodes": nodes_i,
        "path_i": job_0_dir_i}

    if run_calc:
        CC.submit_job(**def_params)

    #__|


def read_write_CONTCAR(path, new_job_file_list):
    # Getting the latest atoms object from CONTCAR
    contcar_file_path = os.path.join(path, "CONTCAR")
    my_file = Path(contcar_file_path)
    if my_file.is_file():
        atoms_for_next_job = io.read(contcar_file_path)
        atoms_for_next_job.write("init.cif")
        new_job_file_list["./init.cif"] = "init.cif"
        # new_job_file_list.append("init.cif")
    else:
        print("Not good! Can't find CONTCAR to restart the job from")
