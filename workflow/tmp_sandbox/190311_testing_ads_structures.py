# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Checking OER Adsorbate Structures
# ___

# # Import Modules

# %%capture
# %load_ext autoreload
# %autoreload 2

# +
# Setting Custom Paths ********************************************************
# *****************************************************************************
import os
import sys

sys.path.insert(0, os.path.join(
        os.environ["PROJ_irox"],
        "data"))

sys.path.insert(0, os.path.join(
        os.environ["PROJ_irox"],
        "workflow"))

# Python Modules **************************************************************
# *****************************************************************************
from ase.visualize import view

# My Modules ******************************************************************
# *****************************************************************************

# Local Imports ***************************************************************
# *****************************************************************************
from an_data_processing import load_df

# Project Data
from proj_data_irox import (
    data_dir,
    )
# -

df_pourbaix, df_ads, df_surf = load_df(
    from_file=False,
    root_dir=data_dir,
    data_dir=data_dir + "/190103_new_job_df",
    file_name="df_master.pickle",
    process_df=True,
    )

# +
df_oh = df_ads[df_ads["adsorbate"] == "oh"]


atoms_list = []
for traj_i in df_oh["atoms_object"].tolist():
    atoms_list.append(traj_i[-1])
# -

view(atoms_list)
