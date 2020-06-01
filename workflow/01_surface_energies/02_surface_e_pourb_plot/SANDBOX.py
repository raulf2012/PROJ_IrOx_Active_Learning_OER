# ---
# jupyter:
#   jupytext:
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

ids = [
  194,
  196,
  197,
  200,
  202,
  294,
  ]

# +
# %%capture
import os
import sys

# #############################################################################
sys.path.insert(0,
    os.path.join(
        os.environ["PROJ_irox"],
        "workflow"))
from an_data_processing import load_df
from an_data_processing import oxy_ref, hyd_ref

# # #############################################################################
import pickle
import pandas as pd

# # #############################################################################
from misc_modules.pandas_methods import drop_columns
from surface_energy.surface_energy import SurfaceEnergy

from proj_data_irox import bulk_e_per_atom_dict
# -

from ase.visualize import view

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
df_m = df_m[df_m["job_type"] == "surface_coverage_energy"]

cols_to_keep = [
    'facet',
    'job_type',
    'layers',
    'surface_type',
    'elec_energy',
    'atoms_object',
    'bulk_system',
    'coverage_type',
    'nonstoich_Os',
    ]

df_m = drop_columns(df=df_m, columns=cols_to_keep, keep_or_drop="keep")

# #############################################################################

bulk_data_path = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/an_bulk_systems",
    "bulk_systems.pickle")
with open(bulk_data_path, "rb") as fle:
    bulk_data = pickle.load(fle)


# +
# df_m = df_m[
#     (df_m["bulk_system"] == "IrO3") &
#     (df_m["facet"] == "111") &
# #     (df_m[""] == "") &
#     [True for i in range(len(df_m))]
#     ]

# +
# # %%capture

def method(row_i):
    """
    """
    print(row_i["bulk_system"] + "_" + row_i["facet"] + "_" + row_i["coverage_type"])


    bulk_e_per_atom = bulk_e_per_atom_dict[row_i["bulk_system"]]

    SE = SurfaceEnergy(
        atoms=row_i["atoms_object"][-1],
        bulk_atoms=bulk_data[row_i["bulk_system"]],
        bulk_electronic_energy_per_atom=bulk_e_per_atom,
        H_ref_electronic_energy=hyd_ref,
        O_ref_electronic_energy=oxy_ref,
        verbose=True,
        )

    print("index: ", row_i.name)
    print("")

    return(SE)

df_m["SurfaceEnergy"] = df_m.apply(
    method,
    axis=1,
    )

# +
import copy

non_stoich_comp = df_m.iloc[1]["SurfaceEnergy"].non_stoich_comp

non_stoich_comp_new = copy.copy(non_stoich_comp)

print(non_stoich_comp)

special_species_dict = dict()
if "O" in non_stoich_comp.keys():

    num_Os = non_stoich_comp.get("O")
    
    if "H" in non_stoich_comp.keys():
        num_Hs = non_stoich_comp.get("H")
        
        min_num = min([num_Os, num_Hs])
        
        num_OHs = min_num

        left_over_Hs = num_Hs - min_num
        left_over_Os = num_Os - min_num

        special_species_dict["*OH"] = num_OHs
        special_species_dict["*O"] = left_over_Os

        non_stoich_comp_new["O"] = 0  # All nonstoich Os will be *O species
        non_stoich_comp_new["H"] = left_over_Hs
    else:
        num_OHs = 0
        special_species_dict["*OH"] = num_OHs

        left_over_Hs = 0
        left_over_Os = num_Os

        special_species_dict["*O"] = left_over_Os
        special_species_dict["*OH"] = 0

        non_stoich_comp_new["O"] = 0  # All nonstoich Os will be *O species
        non_stoich_comp_new["H"] = left_over_Hs
        
else:
    num_OHs = 0
    left_over_Os = num_Os
    left_over_Hs = 0

    if "H" in non_stoich_comp.keys():    
        if non_stoich_comp.get("H") > 0:
            raise ValueError("NOT GOOD HERE, THERE IS AN *H WITHOUT and *OH")
            
print("----")
print(non_stoich_comp_new)
print(special_species_dict)

# +
list0 = [
{
    "index": "orig",
    "O": non_stoich_comp.get("O", 0),
    "H": non_stoich_comp.get("H", 0),
    },
    ]

df = pd.DataFrame(list0)

# df["H_O"] = 

df
# -

min([3, 4])

# + active=""
# O3 H4
#
# O5 H4

# +
# df_i = df_m.loc[ids]

# atoms_list = [i[-1] for i in df_i["atoms_object"].tolist()]
# df_i["atoms"] = atoms_list

# + jupyter={}
# import tempfile
# import shutil

# dirpath = tempfile.mkdtemp(
#     suffix=None,
#     prefix="RAUL_TEMP_DIR_",
#     )

# # dirpath = "/tmp/RAUL_TEMP_DIR_i6m7jdtp"
# # print(dirpath)


# def method(row_i):
#     row_i["atoms"].write(
#         dirpath + "/" + str(row_i.name).zfill(4) + ".cif")

# df_i.apply(
#     method,
#     axis=1)



# # shutil.rmtree(dirpath)

# + active=""
#
#
#
#
#

# +
"IrO3_100_h_covered"

df_tmp = df_m[
    (df_m["bulk_system"] == "IrO3") &
    (df_m["facet"] == "111") &
#     (df_m[""] == "") &
    [True for i in range(len(df_m))]
    ]

df_tmp.loc[202]["SurfaceEnergy"].non_stoich_comp

print(df_tmp["atoms_object"].loc[200][-1])
print(df_tmp["atoms_object"].loc[201][-1])
print(df_tmp["atoms_object"].loc[202][-1])

df_tmp
# -

for i_cnt, row_i in df_tmp.iterrows():
    row_i["SurfaceEnergy"]
#     df_tmp["SurfaceEnergy"]

    non_stoich_comp = row_i["SurfaceEnergy"].non_stoich_comp
    print(non_stoich_comp)

from ase_modules.ase_methods import view_in_vesta

# +
# atoms_0 = df_tmp["atoms_object"].loc[200][-1]
# atoms_1 = df_tmp["atoms_object"].loc[201][-1]

# atoms_list = [
#     df_tmp["atoms_object"].loc[200][-1],
#     df_tmp["atoms_object"].loc[201][-1],
#     df_tmp["atoms_object"].loc[202][-1],
#     ]


# view_in_vesta(
#     atoms_list,
#     ase_gui=True,
#     )

# +
df_tmp = df_m[
    (df_m["bulk_system"] == "IrO3") &
    (df_m["facet"] == "100") &
#     (df_m[""] == "") &
    [True for i in range(len(df_m))]
    ]; df_tmp

atoms_list = []
for i_cnt, row_i in df_tmp.iterrows():
    print(i_cnt)
    print(row_i["SurfaceEnergy"].non_stoich_comp)
    
    atoms_i = row_i["atoms_object"][-1]
    
    atoms_list.append(atoms_i)
    
    
# view_in_vesta(atoms_list, ase_gui=True)

# +
# atoms_i.write("iro3_100_o_covered.cif")
# -

view_in_vesta



# +
row_i = df_m.loc[194]

self = row_i["SurfaceEnergy"]
self.non_stoich
# -

view_in_vesta(row_i["atoms_object"][-1])

# + active=""
#
#
#
#
#
#
#
# -

assert False

# + jupyter={}

# | - IMPORT MODULES
import numpy as np
import pandas as pd

import plotly.graph_objs as go

from ase_modules.ase_methods import create_species_element_dict

from pymatgen.core.composition import Composition
#__|


# + jupyter={}
main_atom = "Ir"  # Make this a class attribute
find_bulk_form_units_method = "main_atom"  # 'gcm' (greatest common multiple)

bulk_atoms = self.bulk_atoms
atoms = self.atoms


comp0 = Composition(bulk_atoms.get_chemical_formula())

df = pd.DataFrame([
    create_species_element_dict(atoms,
        elems_to_always_include=["O", "H"]),
    dict(comp0.to_data_dict["reduced_cell_composition"])],
    index=["slab", "bulk"])

# Replace NaNs with 0.
df = df.replace(np.nan, 0.0, regex=True)

# Removingg columns with 0
df = df.loc[:, (df != 0).any(axis=0)]

# slab_comp_array = np.array(list(df.loc["slab"]))
# bulk_comp_array = np.array(list(df.loc["bulk"]))
# df.loc["slab"].to_numpy()
# df.loc["bulk"].to_numpy()

# Number of unit of the bulk's reduced formula that fit into the slab
if find_bulk_form_units_method == "main_atom": 
    bulk_formula_units_in_slab = int(df.loc["slab"]["Ir"] / df.loc["bulk"]["Ir"])

elif find_bulk_form_units_method == "gcm":
    bulk_formula_units_in_slab = int(min(
        df.loc["slab"].to_numpy() / df.loc["bulk"].to_numpy()
        ))
bfuis = bulk_formula_units_in_slab

# #####################################################################
# Getting the non-stoicheometric atoms composition
df.loc["nonstoich"] = df.loc["slab"] - bfuis * df.loc["bulk"]
non_stoich_comp = df.loc["nonstoich"].to_dict()
self.non_stoich_comp = non_stoich_comp

print(bulk_formula_units_in_slab)
print(non_stoich_comp)
# return(bulk_formula_units_in_slab)


# + jupyter={}
bulk_formula_units_in_slab = int(df.loc["slab"]["Ir"] / df.loc["bulk"]["Ir"])
bfuis = bulk_formula_units_in_slab

df.loc["tmp"] = df.loc["slab"] - bfuis * df.loc["bulk"]
non_stoich_comp = df.loc["tmp"].to_dict()

print(non_stoich_comp)
print(bulk_formula_units_in_slab)

# + jupyter={}
df.loc["slab"].to_numpy()

df.loc["bulk"].to_numpy()
