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
#     display_name: Python [conda env:research-new]
#     language: python
#     name: conda-env-research-new-py
# ---

# # Surface Energy Convergence of IrOx Systems
# ---
#
# Procedure:
# * TMP
# * TMP2
# * TMP3

# + {"toc-hr-collapsed": true, "cell_type": "markdown"}
# # Notebook Setup

# + {"toc-hr-collapsed": true, "cell_type": "markdown"}
# ## Import Modules
# -

# ### Notebook Magik Commands

# %load_ext autoreload
# %autoreload 2

# ### Python Modules

# +
# %%capture
# TEMP
import os
import sys

# #############################################################################
sys.path.insert(0,
    os.path.join(
        os.environ["PROJ_irox"],
        "data"))
from proj_data_irox import irox_bulk_color_map

# #############################################################################
sys.path.insert(0,
    os.path.join(
        os.environ["PROJ_irox"],
        "workflow"))
from an_data_processing import load_df
from an_data_processing import oxy_ref, hyd_ref

# #############################################################################
import pickle
import pandas as pd
import numpy as np
import plotly.graph_objs as go
import plotly.express as px

# #############################################################################
from misc_modules.pandas_methods import drop_columns
from surface_energy.surface_energy import SurfaceEnergy

from plotting.my_plotly import my_plotly_plot
from surface_energy.surface_energy import SurfaceEnergyConvergence as SE_Conv

# #############################################################################
pd.set_option("display.max_columns", None)

# #############################################################################
from IPython.display import display
# -

# ## Script Inputs

verbose = False

# + {"toc-hr-collapsed": true, "cell_type": "markdown"}
# ## Read Data
# -

# ### Read surface energy dataframe

# +
dataframe_dir = os.path.join(
    os.environ["PROJ_DATA"],
    "04_IrOx_surfaces_OER/190321_new_job_df")

df_pourbaix, df_ads, df_surf = load_df(
    from_file=True,
    root_dir=dataframe_dir,
    data_dir=dataframe_dir,
    file_name="df_master.pickle",
    process_df=True)
df_m = df_surf


# Filter the jobs that were unsuccessful
df_m = df_m[[not i for i in pd.isna(df_m["elec_energy"].tolist())]]
df_m = df_m[df_m["job_type"] == "surface_energy"]


cols_to_keep = [
    'facet',
    'job_type',
    'layers',
    'surface_type',
    'elec_energy',
    'atoms_object',
    'bulk_system',
    'coverage_type',
    ]

df_m = drop_columns(df=df_m, columns=cols_to_keep, keep_or_drop="keep")
# -

# ### Read bulk systems data

bulk_data_path = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/an_bulk_systems",
    "bulk_systems.pickle")
with open(bulk_data_path, "rb") as fle:
    bulk_data = pickle.load(fle)


# # -------------------------------------

# # Instantiate SurfaceEnergy to DataFrame

# +
def method(row_i):
    """
    """
    SE = SurfaceEnergy(
        atoms=row_i["atoms_object"][-1],
        bulk_atoms=bulk_data[row_i["bulk_system"]],
        H_ref_electronic_energy=hyd_ref,
        O_ref_electronic_energy=oxy_ref,
        verbose=verbose,
        )

    return(SE)

df_m["SurfaceEnergy"] = df_m.apply(
    method,
    axis=1,
    )

# +
SE_i = df_m.iloc[0]["SurfaceEnergy"]

dir(SE_i)

SE_i.bulk_electronic_energy

SE_i.bulk_electronic_energy_per_atom

SE_i.bulk_atoms.get_number_of_atoms()
# -

assert False
