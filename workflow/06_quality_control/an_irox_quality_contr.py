# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_IrOx_Active_Learning_OER]
#     language: python
#     name: conda-env-PROJ_IrOx_Active_Learning_OER-py
# ---

# + [markdown] Collapsed="false"
# # Checking Job Sets for Consistency and Quality
#
# ***

# + Collapsed="false" active=""
#   The jobs within a Job set (defined as a set of DFT jobs containing bare, ooh, o, oh species)
# should be consistent with one another.
#
# 1. All atoms (other than adsorbate atoms) should be reasonably consistent across all jobs
# 2. Magnetic structures should be reasonably consistent
# 3. Free from unwanted dissocation or other unexpected atom movement
#
#
# ***************************************************************************************
#
# The atoms objects need the magnetic moments set to the init_magmom attribute
#
# ads_magmoms = ads.get_magnetic_moments()
# ads.set_initial_magnetic_moments(ads_magmoms)
#
# ***************************************************************************************

# + [markdown] Collapsed="false"
# # Import Modules

# + Collapsed="false"
# %%capture
# %load_ext autoreload
# %autoreload 2

# + Collapsed="false"
import sys
import os

sys.path.insert(
    0, os.path.join(
    os.environ["PROJ_irox"],
    "workflow"))

sys.path.insert(
    0, os.path.join(
    os.environ["PROJ_irox"],
    "data"))

# + Collapsed="false"
###############################################################################
# Local Imports ###############################################################
from proj_data_irox import data_dir
from an_data_processing import load_df


# #############################################################################
# Python Modules # ############################################################
import subprocess
import datetime

import numpy as np
import pandas as pd

from ase.visualize import view
from ase.io import Trajectory
from ase import io


# #############################################################################
# IPython Imports #############################################################
from IPython.display import display
# -

pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

# + [markdown] Collapsed="false"
# # Script Inputs

# + Collapsed="false"
prop_name_list = [
    'bulk_system',
#     'coverage',
    'coverage_type',
    'facet',
    ]

# + [markdown] Collapsed="false"
# # Read and Process Data Frame

# + Collapsed="false"
df_pourbaix, df_ads, df_surf = load_df(
    from_file=True,
    root_dir=data_dir,
    data_dir=data_dir,
    file_name="df_master.pickle",
    process_df=True,
    )

df_m = df_ads
df_m.loc[df_m["coverage_type"] == "O-4_OH-0", "coverage_type"] = "o_covered"
df_m.loc[df_m["coverage_type"] == "O-2_OH-0", "coverage_type"] = "o_covered_2"
df_m.loc[df_m["coverage_type"] == "O-2_OH-2", "coverage_type"] = "h_covered"


# + Collapsed="false"
def create_energetics_table(
    atoms_dir,
    sys_name_i,
    root_dir,
    ):
    """
    """
    # | - create_energetics_table
    os.chdir(
        os.path.join(
            atoms_dir,
            sys_name_i,
            ),
        )

    bash_comm = os.path.join(
        os.environ["PROJ_irox"],
        "scripts",
        "dft_scripts",
        "get_energies_top_sites.py",
        )

    result = subprocess.run(
        [bash_comm],
        stdout=subprocess.PIPE,
        )

    out_i = result.stdout.decode("utf-8")

    with open("energetics_summary.txt", "w") as text_file:
        text_file.write(out_i)
#     out_i = [i for i in out_i.split("\n") if i != '']

    os.chdir(root_dir)
    #__|


# + [markdown] Collapsed="false"
# # Creating Traj Objects for Each 'Job Set'

# + Collapsed="false"
d = datetime.datetime.today()
date_i = str(d.year) + str(d.month) + str(d.day)

root_dir = os.getcwd()

num_atoms_dirs = len([i for i in os.listdir(".") if "job_sets_atoms" in i])

# atoms_dir = "00_job_sets_atoms/" + str(num_atoms_dirs).zfill(2) + "_job_sets_atoms"
atoms_dir = "out_data/00_job_sets_atoms_" + date_i + "/" + str(num_atoms_dirs).zfill(2) + "_job_sets_atoms"

if not os.path.exists(atoms_dir):
    os.makedirs(atoms_dir)

# + Collapsed="false"
groups = df_m.groupby(by=prop_name_list)
for i_ind, (props_i, group_i) in enumerate(groups):
    sys_name_i = str(i_ind).zfill(2) + "_" + "___".join(props_i)

    sys_i_dir = os.path.join(
        atoms_dir,
        sys_name_i)

    # #########################################################################
    # Create folder ###########################################################
    if not os.path.exists(sys_i_dir): os.makedirs(sys_i_dir)

    for j_cnt, row_j in group_i.iterrows():
        name_j = row_j["adsorbate"]
        if not pd.isna(row_j["ooh_direction"]):
            name_j += "_" + row_j["ooh_direction"]

        name_j += ".json"

        if row_j["atoms_object"] is None:
            continue
        else:
            atoms_j = row_j["atoms_object"][-1]

            atoms_j_full_path = os.path.join(
                atoms_dir,
                sys_name_i,
                name_j)

            print(atoms_j_full_path)

            io.write(
                atoms_j_full_path,
                atoms_j)


    create_energetics_table(
        atoms_dir,
        sys_name_i,
        root_dir)

# +
# images_list = [i for i in group_i["atoms_object"].tolist() if i is not None]
# final_atoms_list = [i[-1] for i in images_list]
