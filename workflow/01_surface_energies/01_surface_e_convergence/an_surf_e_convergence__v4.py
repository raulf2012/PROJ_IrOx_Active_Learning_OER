# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
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

# # Surface Energy Convergence of IrOx Systems
# ---
#
# Procedure:
# * TMP
# * TMP2
# * TMP3

# + [markdown] toc-hr-collapsed=true
# # Notebook Setup

# + [markdown] toc-hr-collapsed=true
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

# + [markdown] toc-hr-collapsed=true
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
# -

# # TEMP | Select subset of data

# +
# # TEMP TEMP
# df_m = df_m[
#     (df_m["bulk_system"] == "IrO2") &
# #     (df_m["bulk_system"] == "IrO3_rutile-like") &
# #     (df_m["facet"] == "001") &
#     [True for i in range(len(df_m))]
#     ]

# + [markdown] toc-hr-collapsed=true
# # Surface Energy vs Slab Width Plot <------------------------
#
# Explain what's happenging here
# -

SEC_data = []
data = []
grouped = df_m.groupby(["bulk_system", "facet"])
for i_cnt, (name, group) in enumerate(grouped):

    SEC = SE_Conv(
        SurfaceEnergy_instances=group["SurfaceEnergy"].tolist(),
        verbose=verbose,
        )
    self = SEC

    # Fit bulk and then use to recalculate the surface energies for all slabs
    self.fit_bulk_energy()
    self.calculate_surface_energies(bulk_energy=self.fitted_bulk_energy)

    color_i = irox_bulk_color_map[name[0]]
    name_i = "_".join(list(name))

    SEC_data.append({
        "name": name_i,
        "bulk_system": name[0],
        "facet": name[1],
        "SEC": SEC})

    data_i = self.plot_surface_energy(
        name_i=name_i,
        color_i=color_i)

    data += data_i

# ## Plot

my_plotly_plot(plot_name="TEMP_PLOT",
    save_dir=None, data=data, upload_plot=False)

# + active=""
#
#
#
#
# -

# # -------------------------------------

# + [markdown] toc-hr-collapsed=true
# # Averaging the fitted bulk energies across different facets
# -

# ## Setting up new DataFrame

# +
df_new = pd.DataFrame(SEC_data)

# #############################################################################
def method(row_i):
    new_se = row_i["SEC"].fitted_bulk_energy
    return(new_se)
df_new["fitted_bulk_energy"] = df_new.apply(
    method, axis=1)

# #############################################################################
def method(row_i):
    new_se = row_i["SEC"].ave_surface_energy_per_area
    return(new_se)
df_new["dft_bulk_ave_surface_e"] = df_new.apply(
    method, axis=1)

def method(row_i):
    new_se = row_i["SEC"].new_ave_surface_energy_per_area
    return(new_se)
df_new["fitted_bulk_ave_surface_e"] = df_new.apply(
    method, axis=1)
# -

# ## Group bulk fitted energies, average and std dev.
#
# The fitted bulk energies should hopefully be the same across fits done for different surfaces

final_averaged_fitted_bulk_energies = dict()
grouped = df_new.groupby(["bulk_system"])
for i_cnt, (name, group) in enumerate(grouped):
    # display(group)

    print(80 * "_")
    name_tmp = "|  " + name + "  |"
    print(name_tmp)
    print(len(name_tmp) * "-")

    ave_fitted_bulk_energy = group["fitted_bulk_energy"].mean()
    print("ave_fitted_bulk_energy:", ave_fitted_bulk_energy)

    std_fitted_bulk_energy = group["fitted_bulk_energy"].std()
    print("std_fitted_bulk_energy:", std_fitted_bulk_energy)

    final_averaged_fitted_bulk_energies[name] = ave_fitted_bulk_energy

# + active=""
# ________________________________________________________________________________
# |  IrO2  |
# ----------
# ave_fitted_bulk_energy: -7.044324480309523
# std_fitted_bulk_energy: 5.383171936804735e-05
# ________________________________________________________________________________
# |  IrO3  |
# ----------
# ave_fitted_bulk_energy: -6.463874089016009
# std_fitted_bulk_energy: 0.0017695334929991133
# ________________________________________________________________________________
# |  IrO3_rutile-like  |
# ----------------------
# ave_fitted_bulk_energy: -6.458260171334821
# std_fitted_bulk_energy: 0.0001088472764921181

# + active=""
# These results look pretty good
# The fitted bulk energies are very close to one another for each bulk system

# + active=""
#
#
#
# -

# # Recalculate all surface energies with newly average fitted bulk energetics

SEC_data = []
data = []
grouped = df_m.groupby(["bulk_system", "facet"])
for i_cnt, (name, group) in enumerate(grouped):
    bulk_system_i = name[0]
    facet_i = name[1]

    SEC = SE_Conv(
        SurfaceEnergy_instances=group["SurfaceEnergy"].tolist(),
        verbose=verbose,
        )
    self = SEC

    ave_fitt_bulk_energy = final_averaged_fitted_bulk_energies[bulk_system_i]

    # Use the averaged fitted bulk  energies
    self.calculate_surface_energies(
        bulk_energy=ave_fitt_bulk_energy
        )
    print(name_i)
    print(SEC.new_ave_surface_energy_per_area)

    color_i = irox_bulk_color_map[name[0]]
    name_i = "_".join(list(name))

    SEC_data.append({
        "name": name_i,
        "bulk_system": name[0],
        "facet": name[1],
        "SEC": SEC})

    data_i = self.plot_surface_energy(
        name_i=name_i,
        color_i=color_i)

    data += data_i

my_plotly_plot(plot_name="TEMP_PLOT",
    save_dir=None, data=data, upload_plot=False)

# + active=""
#
#
#
#
#

# + jupyter={}
#     # #########################################################################
#     # | - Surface Energy (DFT Bulk)
#     y_surface_e = []; x_slab_thickness = []
#     for SE_inst_i in self.SurfaceEnergy_instances:
#         y_surface_e.append(SE_inst_i.surface_e_per_area)
#         x_slab_thickness.append(SE_inst_i.slab_thickness)

#     trace_i = go.Scatter(
#         x=x_slab_thickness,
#         y=y_surface_e,
#         mode='markers+lines',
#         name=name_i,
#         legendgroup=name_i,
#         showlegend=True,
#         line=dict(
#             width=1.5,
#             color=color_i,
#             dash='dash',
#             ),
#         marker=dict(
#             symbol="square",
#             size=8,
#             color=color_i,
#             line=dict(
#                 width=1.0,
#                 color="black",
#                 ),
#             ),

#         )
#     data.append(trace_i)
#     #__|

#     # #########################################################################
#     # | - Surface Energy (Fitted Bulk)
#     y_surface_e = []; x_slab_thickness = []
#     for SE_inst_i in self.new_SurfaceEnergy_instances:
#         y_surface_e.append(SE_inst_i.surface_e_per_area)
#         x_slab_thickness.append(SE_inst_i.slab_thickness)

#     trace_i = go.Scatter(
#         x=x_slab_thickness,
#         y=y_surface_e,
#         mode='markers+lines',
#         name=name_i,
#         # legendgroup=name_i,
#         showlegend=False,
#         line=dict(
#             width=1.5,
#             color=color_i,
#             dash='solid',
#             ),
#         marker=dict(
#             symbol="circle",
#             size=10,
#             color=color_i,
#             line=dict(
#                 width=1.5,
#                 color="black",
# #                 dash='solid',
#                 ),
#             ),
#         )
#     data.append(trace_i)
#     #__|

#     # #########################################################################
#     # | - Average Surface Energy (DFT Bulk)
#     ave_surface_energy = self.ave_surface_energy_per_area
#     trace_i = go.Scatter(
#         x=[0, 30],
#         y=[ave_surface_energy, ave_surface_energy],
#         name=name_i,
#         # legendgroup=name_i,
#         mode='lines',
#         showlegend=False,
#         line=dict(
#             width=1.5,
#             color=color_i,
#             dash='dash',
#             ),
#         marker=dict(
#             symbol="square",
#             size=10,
#             color=color_i,
#             line=dict(
#                 width=1,
#                 # color='rgb(0, 0, 0)',
#                 color=color_i,
#                 ),
#             ),
#         )
#     data.append(trace_i)
#     #__|

#     # #########################################################################
#     # | - Average Surface Energy (Fitted Bulk)
#     ave_surface_energy = self.new_ave_surface_energy_per_area
#     trace_i = go.Scatter(
#         x=[0, 30],
#         y=[ave_surface_energy, ave_surface_energy],
#         name=name_i,
#         # legendgroup=name_i,
#         mode='lines',
#         showlegend=False,
#         line=dict(
#             width=1.5,
#             color=color_i,
#             dash='solid',
#             ),
#         marker=dict(
#             symbol="square",
#             size=10,
#             color=color_i,
#             line=dict(
#                 width=1,
#                 color="black",
#                 ),
#             ),
#         )
#     data.append(trace_i)
#     #__|

# + jupyter={}
# TEMP
# # self.fit_bulk_energy()
# # self.calculate_surface_energies(bulk_energy=self.fitted_bulk_energy)
# # self.calculate_surface_energies(bulk_energy=-10)

# self.new_SurfaceEnergy_instances
# self.ave_surface_energy_per_area
