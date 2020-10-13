#!/usr/bin/env python

"""TEMP.

Author(s): Raul A. Flores
"""

# | - IMPORT MODULES
import os
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

# | - Job parsing methods
# #########################################################
def parse_job_err(path):
    status_dict = {
        "timed_out": None,
        "error": None,
        "error_type": None,
        }

    # print(2 * "")
    # print("methds | parse_job_err")

    compenv = os.environ["COMPENV"]
    # print("compenv:", compenv)

    # | - Parsing SLAC job
    if compenv == "slac":
        job_out_file_path = os.path.join(path, "job.out")
        my_file = Path(job_out_file_path)
        if my_file.is_file():
            with open(job_out_file_path, 'r') as f:
                lines = f.readlines()

            for line in lines:
                if "job killed after reaching LSF run time limit" in line:
                    print("Found following line in job.err")
                    print("job killed after reaching LSF run time limit")
                    status_dict["timed_out"] = True
                    break
    #__|

    # | - Parsing error file
    job_err_file_path = os.path.join(path, "job.err")
    my_file = Path(job_err_file_path)
    if my_file.is_file():
        with open(job_err_file_path, 'r') as f:
            lines = f.readlines()

        # else:
        for line in lines:
            if "DUE TO TIME LIMIT" in line:
                status_dict["timed_out"] = True

            if "Traceback (most recent call last):" in line:
                status_dict["error"] = True

        # if compenv == "slac":
        #     for line in lines:
        #         if "job killed after reaching LSF run time limit" in line:
        #             print("Found following line in job.err")
        #             print("job killed after reaching LSF run time limit")
        #             status_dict["timed_out"] = True

    #__|


    # | - Parsing out file
    if status_dict["error"] is True:
        job_out_file_path = os.path.join(path, "job.out")
        my_file = Path(job_out_file_path)
        if my_file.is_file():
            with open(job_out_file_path, 'r') as f:
                lines = f.readlines()

            for line in lines:
                err_i = "VERY BAD NEWS! internal error in subroutine SGRCON:"
                if err_i in line:
                    status_dict["error_type"] = "Error in SGRCON (symm error)"
                    break
    #__|

    # print(2 * "")

    return(status_dict)

# #########################################################
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
    # print("")
    # print(path)
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

# #############################################################################
def get_number_of_ionic_steps(path):
    status_dict = {"num_steps": None}

    outcar_path = os.path.join(path, "OUTCAR")
    try:
        traj = io.read(outcar_path, index=":")
        status_dict["num_steps"] = len(traj)
    except:
        status_dict["num_steps"] = 0
        pass

    return(status_dict)

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
    # | - set_up__submit__new_job
    path = latest_revision["path"]
    pre_path = latest_revision["pre_path"]
    current_revision = latest_revision["revision"]
    num_prev_steps = latest_revision["num_steps"]

    job_state = latest_revision["job_state"]
    timed_out = latest_revision["timed_out"]
    isif = latest_revision["isif"]


    new_revision = current_revision + 1
    job_0_dir_i = os.path.join(pre_path, "_" + str(new_revision))


    try:
        os.makedirs(job_0_dir_i)
    except:
        pass



    # | - Copy files into new directory
    for file_path_i, dist_name in new_job_file_dict.items():
        copyfile(file_path_i, os.path.join(job_0_dir_i, dist_name))
    #__|


    dft_params_dict = {
        # "encut": 600,
        # "kpar": 5,
        # "ediffg": 5e-3,
        # "ediff": 1e-6
        }

    if os.environ["COMPENV"] == "sherlock":
        print("iksfijsijfisddfi8998y0934389 | TEMP")
        print("In Sherlock")
        dft_params_dict["npar"] = 4


    # Write atoms object
    # atoms_i = row_i["init_atoms"]
    # io.write(os.path.join(job_0_dir_i, "init.cif"), row_i["init_atoms"])

    atoms_i = io.read(os.path.join(job_0_dir_i, "init.cif"))

    if num_prev_steps < 10:
        wall_time_factor = 2.5
    elif num_prev_steps < 4:
        wall_time_factor = 3
    else:
        wall_time_factor = 1.8

    num_atoms = atoms_i.get_number_of_atoms()
    wall_time_i = calc_wall_time(num_atoms, factor=wall_time_factor)
    wall_time_i = int(wall_time_i)

    if os.environ["COMPENV"] != "slac":
        if wall_time_i > 600:
            wall_time_i = 600
    else:
        wall_time_i = 8. * wall_time_i

        if wall_time_i > 2760:
            wall_time_i = 2760


    if num_atoms > 100:
        nodes_i = 10
        dft_params_dict["kpar"] = 10
    else:
        nodes_i = 5
        dft_params_dict["kpar"] = 5

    if os.environ["COMPENV"] == "slac":
        dft_params_dict["kpar"] = 3
        dft_params_dict["npar"] = 4



    # Write dft paramters json file to job dir
    with open(os.path.join(job_0_dir_i, "dft-params.json"), "w+") as fle:
        json.dump(
            dft_params_dict, fle,
            indent=2, skipkeys=True)

    # # Copy model job script
    # copyfile(model_file, os.path.join(job_0_dir_i, "model.py"))


    # Submit job ##############################################################
    CC = ComputerCluster()

    wall_time_i = int(wall_time_i)

    def_params = {
        "wall_time": wall_time_i,
        # "queue": "premium",
        "queue": "regular",
        "architecture": "knl",
        "nodes": nodes_i,
        "priority": "scavenger",
        "path_i": job_0_dir_i}

    if os.environ["COMPENV"] == "sherlock":
        print("Found SHERLOCK env var")

        def_params = {
            "wall_time": wall_time_i,
            "nodes": nodes_i,
            "queue": "iric",
            "path_i": job_0_dir_i}
    elif os.environ["COMPENV"] == "slac":
        print("Found SLAC env var")

        def_params = {
            "wall_time": wall_time_i,
            "cpus": 12,
            "queue": "suncat2",
            "path_i": job_0_dir_i}

    if run_calc:
        CC.submit_job(**def_params)

    #__|


def read_write_CONTCAR(path, new_job_file_list):
    """
    """
    # | - read_write_CONTCAR
    # Getting the latest atoms object from CONTCAR
    contcar_file_path = os.path.join(path, "CONTCAR")
    contcar_file = Path(contcar_file_path)

    init_file_path = os.path.join(path, "init.cif")
    init_file = Path(init_file_path)

    contcar_parsed = False
    if contcar_file.is_file():
        try:
            atoms_for_next_job = io.read(contcar_file_path)
            atoms_for_next_job.write("init.cif")
            new_job_file_list["./init.cif"] = "init.cif"
            contcar_parsed = True
        except:
            pass

    if init_file.is_file() and not contcar_parsed:
        try:
            print("Using previous init.cif not CONTCAR")
            atoms_for_next_job = io.read(init_file_path)
            atoms_for_next_job.write("init.cif")
            new_job_file_list["./init.cif"] = "init.cif"
        except:
            print("Attempted to read init.cif but failed")


    # else:
    #     print("Failed reading CONTCAR and init.cif not parsable for some reasons")
    #
    # else:
    #     print("Not good! Can't find CONTCAR or init.cif to restart from")
    #__|
