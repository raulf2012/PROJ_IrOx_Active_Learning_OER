# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernel_info:
#     name: python3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Import Modules

# %%capture
# %load_ext autoreload
# %autoreload 2

# +
import sys
import os

sys.path.insert(
    0,
    os.path.join(
        os.environ["PROJ_irox"],
        "scripts",
        ),
    )

sys.path.insert(
    0,
    os.path.join(
        os.environ["PROJ_irox"],
        "workflow/an_analysis_dir",
        ),
    )

sys.path.insert(
    0,
    os.path.join(
        os.environ["PROJ_irox"],
        "workflow/data",
        ),
    )

from Pourbaix_simple_generic import plot_Pourbaix
from an_data_processing import load_df

from ase import io
from ase.visualize import view

from ase_modules.ase_methods import create_species_element_dict

# +
from data import h2_ref as h2
from data import h2o_ref as h2o

from data import zpe_h2o
from data import cv_h2o
from data import ts_h2o

from data import zpe_h2
from data import cv_h2
from data import ts_h2

from data import zpe_ooh
from data import cv_ooh
from data import ts_ooh

from data import zpe_o
from data import cv_o
from data import ts_o

from data import zpe_oh
from data import cv_oh
from data import ts_oh
# -

# # Script Inputs

# +
close_plt = False

save_dir = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/pl_master_plots/pl_pourbaix"
    )

Umin=0.0
Umax=2.2
# -

Pourbaix_arg_dict = {
#     "surfs": surfs,
    "h2": h2,
    "zpe_h2": zpe_h2,
    "ts_h2": ts_h2,
    "cv_h2": cv_h2,
    "h2o": h2o,
    "zpe_h2o": zpe_h2o,
    "ts_h2o": ts_h2o,
    "cv_h2o": cv_h2o,
    "zpe_o": zpe_o,
    "ts_o": ts_o,
    "cv_o": cv_o,
    "zpe_oh": zpe_oh,
    "ts_oh": ts_oh,
    "cv_oh": cv_oh,
    "zpe_ooh": zpe_ooh,
    "ts_ooh": ts_ooh,
    "cv_ooh": cv_ooh,
    "Umin": Umin,
    "Umax": Umax,
    "print_out": False,
    "save_dir": save_dir,
#     "file_name": "_".join(list(key_i)) + ".pdf",
    "close_plt": close_plt,
    }

# # Read/Process DataFrame

# +
df_pourbaix, df_ads, df_surf = load_df(
    from_file=False,
    root_dir="../data",
    data_dir="../data",
    file_name="df_master.pickle",
    process_df=True,
    )

# df_m = df_pourbaix

# # Elimate structures that aren't converged w.r.t. forces
# df_m = df_m[df_m["max_force"] < 0.05]

# df_m["name_i"] = df_m["name_i"].str.replace("_", " ")
# df_m["name_i"] = df_m["name_i"].str.replace("|", ",")

# grouped = df_m.groupby(["facet", "bulk_system"])
# group_dict = {}
# for i_ind, (name, group) in enumerate(grouped):
#     df_i = group
#     group_dict[name] = group

# +
atoms_i = df_ads[
    (df_ads["coverage_type"] == "o_covered") & \
    (df_ads["adsorbate"] == "ooh") & \
    (df_ads["bulk_system"] == "IrO2") & \
    (df_ads["facet"] == "110")
#     (df_ads[""] == ) & \
    ].iloc[0].atoms_object[-1]

io.write("iro2_110.traj", atoms_i)


atoms_i = df_ads[
    (df_ads["coverage_type"] == "o_covered") & \
    (df_ads["adsorbate"] == "ooh") & \
    (df_ads["bulk_system"] == "IrO2") & \
    (df_ads["facet"] == "100")
#     (df_ads[""] == ) & \
    ].iloc[0].atoms_object[-1]

io.write("iro2_110.traj", atoms_i)


# view(atoms_i)

# +
import pandas as pd


frames = [df_ads, df_surf]
result = pd.concat(frames, sort=False)

result
# -

df_surf[df_surf["max_force"] > 0.02].path.tolist()


df_pourbaix[df_pourbaix["max_force"] > 0.02].path.tolist()

df_ads[df_ads["max_force"] > 0.02]

# # Import Modules

# %%capture
# %load_ext autoreload
# %autoreload 2

# +
import sys
import os

sys.path.insert(
    0,
    os.path.join(
        os.environ["PROJ_irox"],
#         "scripts/01_Michal_OER_Plot_Script",
        "scripts",
        ),
    )

sys.path.insert(
    0,
    os.path.join(
        os.environ["PROJ_irox"],
        "workflow/an_analysis_dir",
        ),
    )

sys.path.insert(
    0,
    os.path.join(
        os.environ["PROJ_irox"],
        "workflow/data",
        ),
    )

from Pourbaix_simple_generic import plot_Pourbaix
from an_data_processing import load_df

from ase_modules.ase_methods import create_species_element_dict

# +
import pickle

from ase.visualize import view
from ase import io

import plotly as py

# +
process_df = True




df_pourbaix, df_ads, df_surf = load_df(
    from_file=False,
    root_dir="../data",
    data_dir="../data",
    file_name="df_master.pickle",
    process_df=True,
    )
df_m = df_pourbaix


# # df_pourbaix, df_ads = load_df(
# df_m = load_df(
#     from_file=False,
#     root_dir=".",
#     data_dir=".",
#     file_name="df_master.pickle",
#     process_df=process_df,
#     )

# if process_df:
#     df_pourbaix, df_ads = df_m

# df_m = df_ads

# df_m = df_m[[
#     "adsorbate",
#     "bulk_system",
#     "coverage",
#     "coverage_type",
#     "facet",
#     "job_type",
#     "path_short",
#     "elec_energy",
#     "name_i",
#     "ads_e",
#     "max_force",
#     "sum_force",
#     ]]

# +
df_m = df_pourbaix


df_i = df_m[
    df_m["bulk_system"] == "IrO3_rutile-like"
    ]

df_i["path_short"].tolist()
# -

df_m[
    (df_m["bulk_system"] == "IrO3") & \
#     (df_m["bulk_system"] == "IrO3")
    (df_m["facet"] == "111")
    
    ]

# +
# grouped = df_m.groupby(["facet", "coverage_type", "bulk_system"])
# group_df_list = []
# for i_ind, (name, group) in enumerate(grouped):
#     df_i = group
# #     print(df_i)
#     print(len(df_i))
#     print(20 * "*")
# -

# **********
# **********
# **********
# **********

# # Figuring out the discrepancy between Michals IrO2 110 and mine

# +
df_tmp = df_m[
    (df_m["facet"] == "100") & \
    (df_m["bulk_system"] == "IrO2") & \
    (df_m["coverage_type"] == "h_covered")
    ]

df_tmp
# -

# # Reference State Energetics

# +
from energetics.dft_energy import Element_Refs

# VASP Gas-phase References
h2_ref = -6.77014123
h2o_ref = -14.21744725

# h2_ref =  -6.759300
# h2o_ref = -14.019771

# Free Energy Corrections
corrections_dict = {
    "ooh": 0.34475,
    "o": -0.0145,
    "oh": 0.30225,
    "bare": 0.,
    }

Elem_Refs = Element_Refs(
    H2O_dict={
        "gibbs_e": h2o_ref,
        "electronic_e": h2o_ref,
        },

    H2_dict={
        "gibbs_e": h2_ref,
        "electronic_e": h2_ref,
        },
    )

oxy_ref, hyd_ref = Elem_Refs.calc_ref_energies()

oxy_ref = oxy_ref.gibbs_e
hyd_ref = hyd_ref.gibbs_e
# -

# # My Calculations

# +
# atoms_OH = df_tmp.loc[25]["atoms_object"][-1]
# atoms_O = df_tmp.loc[23]["atoms_object"][-1]
# atoms_bare = df_tmp.loc[19]["atoms_object"][-1]

# E_O_mine = atoms_O.get_potential_energy()
# E_OH_mine = atoms_OH.get_potential_energy()
# E_bare_mine = atoms_bare.get_potential_energy()

# # E_ads_OH_mine = E_OH_mine - E_bare_mine - oxy_ref - hyd_ref # + corrections_dict["oh"]
# E_ads_OH_mine = -346.7727701 - E_bare_mine - oxy_ref - hyd_ref #+ corrections_dict["oh"]
# E_ads_O_mine = E_O_mine - E_bare_mine - oxy_ref #+ corrections_dict["o"]

# print("E_ads_OH: ")
# print(E_ads_OH_mine)

# print(" ")
# print("E_ads_O: ")
# print(E_ads_O_mine)

# print(" ")
# print("E_O - E_OH: ")
# E_ads_O_mine - E_ads_OH_mine
# -

# # Michal's IrO2 (110) Calculation

# +
# proj_root = os.environ["PROJ_irox"]
# Mich_dir = os.path.join(proj_root, "workflow/Michal_IrO2_110_calc")

# atoms_O_Mich = io.read(os.path.join(Mich_dir, "03_o/IrO2_PBE_opt_Oads_OUTCAR"))
# atoms_OH_Mich = io.read(os.path.join(Mich_dir, "04_oh/IrO2_PBE_opt_OHads_OUTCAR"))
# atoms_bare_Mich = io.read(os.path.join(Mich_dir, "01_bare/IrO2_PBE_opt_clean_OUTCAR"))

# E_O_Mich = atoms_O_Mich.get_potential_energy()
# E_OH_Mich = atoms_OH_Mich.get_potential_energy()
# E_bare_Mich = atoms_bare_Mich.get_potential_energy()

# E_ads_OH_Mich = E_OH_Mich - E_bare_Mich - oxy_ref - hyd_ref #+ corrections_dict["oh"]
# E_ads_O_Mich = E_O_Mich - E_bare_Mich - oxy_ref #+ corrections_dict["o"]

# print("E_ads_OH: ")
# print(E_ads_OH_Mich)

# print(" ")
# print("E_ads_O: ")
# print(E_ads_O_Mich)

# print(" ")
# print("E_O - E_OH: ")
# E_ads_O_Mich - E_ads_OH_Mich

# +
# for i_ind, row_i in df_m.iterrows():

#     atoms_i = row_i.atoms_object[-1]
#     name_i = row_i.name_i_2
    
#     atoms_i.write("atoms_objects_ads/" + name_i + ".cif")
# -

name_i = row_i.name_i_2

