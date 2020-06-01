# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_IrOx_Active_Learning_OER]
#     language: python
#     name: conda-env-PROJ_IrOx_Active_Learning_OER-py
# ---

# %load_ext autoreload
# %autoreload 2

# +
import os
print(os.getcwd())
import sys

sys.path.insert(
    0, os.path.join(
    os.environ["PROJ_irox"],
    "workflow"))

sys.path.insert(
    0, os.path.join(
    os.environ["PROJ_irox"],
    "data"))

import numpy as np

# #############################################################################
# My IMports ##################################################################
from misc_modules.pandas_methods import drop_columns

# #############################################################################
# Project Data ################################################################
from proj_data_irox import (
    smart_format_dict,
    data_dir,
    groupby_props,
    )

# #############################################################################
# Local Imports ###############################################################
from an_data_processing import load_df
# __|

# +
prop_name_list = [
    'bulk_system',
    'coverage_type',
    'facet',
    ]

SC_PLT_share_props = dict(
    num_round=2)

fit_lines_shared = dict(width=1)

# +
# %%capture

df_pourbaix, df_ads, df_surf = load_df(
    from_file=False,
    root_dir=data_dir,
    data_dir=data_dir,
    file_name="df_master.pickle",
    process_df=True)

df_m = df_ads

# +
from misc_modules.pandas_methods import drop_columns

columns_to_keep = [
    'bulk_system',
    'facet',
    'adsorbate',
    'coverage_type',
    'ooh_direction',
    'ads_e',
    # 'elec_energy',
    # 'total_magmom',
    # 'abs_magmom',
    # 'path_short',
    # 'name_i',
    # 'max_force',
    # 'sum_force',
    # 'elem_num_dict',
    # 'incar_parsed',
    'init_atoms',
    # 'atoms_object',
    # 'N_atoms',
    # 'dipole_correction',
    # 'path',
    # 'name_i_2',
    # 'name_i_3',
    # 'priority',
    # 'surface_type',
    ]

df_m = drop_columns(
    df=df_m,
    columns=columns_to_keep,            
    keep_or_drop="keep")


# +
# def method(row_i, argument_0, optional_arg=None):
def method(row_i):
    """
    """
    atoms = row_i.init_atoms

    positions = atoms.get_positions()
    z_array = positions[:, 2]

    cell = atoms.get_cell()
    cell_z = cell[2][2]

    slab_thickness = np.abs(z_array.max() - z_array.min())

    vacuum = cell_z - slab_thickness

    return(vacuum)

df_i = df_m
df_i["vacuum_z"] = df_i.apply(
    method,
    axis=1)
df_m = df_i

# + jupyter={"outputs_hidden": true}
df_m = df_m[df_m.adsorbate != "ooh"]

# +
print("min:", df_m.vacuum_z.min())
print("max:", df_m.vacuum_z.max())

print("mean:", df_m.vacuum_z.mean())

df_m.sort_values("vacuum_z")

# + active=""
#
#
#
#

# + jupyter={}
# row_i = df_ads.iloc[0]

# atoms = row_i.init_atoms

# positions = atoms.get_positions()

# z_array = positions[:, 2]

# cell = atoms.get_cell()
# cell_z = cell[2][2]

# slab_thickness = np.abs(z_array.max() - z_array.min())

# vacuum = cell_z - slab_thickness

# vacuum

# + jupyter={}

# dir(atoms)

# + jupyter={}
# import plotly.graph_objects as go
# import numpy as np

# trace = go.Histogram(
#     x=z_array,
#     )

# data = [trace]

# fig = go.Figure(data=data)

# fig.show()

# + jupyter={}
# #############################################################################
# Python Modules ##############################################################
# import numpy as np
# import plotly.graph_objs as go

# #############################################################################
# My Modules ##################################################################
#  from oxr_reaction.oxr_rxn import ORR_Free_E_Plot
#  from oxr_reaction.oxr_plotting_classes.oxr_plot_scaling import (
#      Scaling_Relations_Plot)

# from plotting.my_plotly import my_plotly_plot, add_duplicate_axes
