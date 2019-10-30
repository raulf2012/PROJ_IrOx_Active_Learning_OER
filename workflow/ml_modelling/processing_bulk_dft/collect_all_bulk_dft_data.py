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
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Import Modules

# +
import os
import sys

import pickle

import numpy as np
import pandas as pd

# #############################################################################
sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))

from proj_data_irox import (
    bulk_dft_data_path,
    unique_ids_path,
    prototypes_data_path,
    static_irox_structures_path)

print(os.getcwd())
# -

# # Read Data

path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling",
    "processing_bulk_dft/parse_chris_bulk_dft/out_data",
    "df_dft_calcs.pickle")
with open(path_i, "rb") as fle:
    df_chris = pickle.load(fle)
df_chris["source"] = "chris"

# +
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/processing_bulk_dft",
    "parse_my_bulk_dft/out_data",
    "df_bulk_raul_irox.pickle")
with open(path_i, "rb") as fle:
    df_raul_irox = pickle.load(fle)
    
df_raul_irox["source"] = "raul"

# +
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/processing_bulk_dft",
    "parse_my_oer_bulk_dft/out_data",
    "df_oer_bulk.pickle")

with open(path_i, "rb") as fle:
    df_oer_bulk = pickle.load(fle)

# df_raul_irox["source"] = "raul"

# +
from proj_data_irox import oqmd_irox_data_path
with open(oqmd_irox_data_path, "rb") as fle:
    df_oqmd_data = pickle.load(fle)

df_oqmd_data = df_oqmd_data.drop(
    labels=[
        "source",
#         "id_unique",
        ],
    axis=1,
    )

df_oqmd_data["source"] = "oqmd"
# -

# # Combining Chris and Raul data

# +
frames = [
    df_raul_irox,
    df_oer_bulk,
    df_chris,
    df_oqmd_data,
    ]

df_m = pd.concat(frames)
# -

# # Mapping unique ID scheme

# +
df_id = pd.read_csv(unique_ids_path)

id_mapp_iro2 = dict(zip(
    df_id[df_id["stoich"] == "AB2"]["id"],
    df_id[df_id["stoich"] == "AB2"]["unique_ids"]))

id_mapp_iro3 = dict(zip(
    df_id[df_id["stoich"] == "AB3"]["id"],
    df_id[df_id["stoich"] == "AB3"]["unique_ids"]))

# #############################################################################

def method(row_i):
    """
    """

    if row_i["source"] == "raul_oer":
        id_unique_i = row_i.name
    else:
        id_i = row_i["id_old"]

        if row_i["stoich"] == "IrO2" or row_i["stoich"] == "AB2":
            mapping_dict = id_mapp_iro2
        elif row_i["stoich"] == "IrO3" or row_i["stoich"] == "AB3":
            mapping_dict = id_mapp_iro3
        else:
            print("BAD BAD | Couldn't process id: ", row_i)

        id_unique_i = mapping_dict[id_i]

    return(id_unique_i)

df_m["id_unique"] = df_m.apply(method, axis=1)
df_m.set_index("id_unique", inplace=True)


# -

# # Adding energy per atom column

# +
def method(row_i):
    atoms_i = row_i["atoms"]
    energy = None
    if atoms_i is None:
        energy = None
    else:
        try:
            energy = atoms_i.get_potential_energy()
        except:
            energy = None
    return(energy)
df_m["energy"] = df_m.apply(method, axis=1)

def method(row_i):
    energy_norm_i = None

    atoms_i = row_i["atoms"]
    # energy_pa = row_i["energy_pa"]
    energy_pa = row_i.get("energy_pa", np.nan)
    

    if not np.isnan(energy_pa) and row_i["source"] == "oqmd":
        energy_norm_i = energy_pa

    else:
        if atoms_i is None:
            energy_norm_i = None
        else:
            num_atoms_i = len(atoms_i.get_atomic_numbers())
            energy_norm_i = row_i["energy"] / num_atoms_i

    return(energy_norm_i)
df_m["energy_pa"] = df_m.apply(method, axis=1)
# -

# # Adding Formation Enthalpy and Gibbs Free Energy

# +
from proj_data_irox import calc_dH


def method(row_i, calc_dH):
    energy_pa = row_i["energy_pa"]
    stoich = row_i["stoich"]
    
    dH = calc_dH(energy_pa, stoich=stoich)
    
    return(dH)

df_m["dH"] = df_m.apply(method, args=(calc_dH, ), axis=1)

# +
# def method(row_i, argument_0, optional_arg=None):
#     """
#     """
#     return(argument_0)

# arg1 = "TEMP_0"
# df_i = model_i
# df_i["column_name"] = df_i.apply(
#     method,
#     axis=1,
#     args=(arg1, ),
#     optional_arg="TEMP_1"
#     )
# -

# # Removing rows with missing atoms objects

df_m = df_m[df_m["atoms"].notnull()]


# # Count number of atoms

# +
def method(row_i):
    atoms_i = row_i["atoms"]
    num_atoms_i = atoms_i.get_number_of_atoms()
    return(num_atoms_i)


df_i = df_m
df_i["num_atoms"] = df_i.apply(
    method,
    axis=1)
df_m = df_i
# -

# # Save data

# +
directory = "out_data"
if not os.path.exists(directory):
    os.makedirs(directory)

with open(os.path.join(directory, "df_bulk_dft.pickle"), "wb") as fle:
    pickle.dump(df_m, fle)
# -
3 * -7.049 - (2 * -4.657947279999998 + -9.304929736367313)

-7.04

# +
df_m[
    (df_m["stoich"] == "AB2") & \
#     (df_m["stoich"] == "AB2")
    (df_m["source"] == "raul")
    ].sort_values("energy_pa").loc[
    
    [
        
"cubqbpzd7k",
"6r716sxr9t",        
    ]
]

# +
print(-7.049062 - -7.040906)

print(-2.515127 - -2.490658)

# + {"active": ""}
#
#
#
#


# + {"jupyter": {"source_hidden": true}}
# df_m[
#     (df_m["source"] == "raul") & \
#     (df_m["stoich"] == "AB2")
#     ]["source"]

# # df_m.head()

# + {"jupyter": {"source_hidden": true}}
# df_m[df_m["source"] == "raul_oer"].index.tolist()

# + {"jupyter": {"source_hidden": true}}
# def method(row_i):
#     atoms_i = row_i["atoms"]
#     num_atoms_i = atoms_i.get_number_of_atoms()
#     return(num_atoms_i)


# df_i = df_m
# df_i["num_atoms"] = df_i.apply(
#     method,
#     axis=1)
# df_m = df_i

# + {"jupyter": {"source_hidden": true}}
# df_m.sort_values("num_atoms").iloc[0]["path"]

# df_m.loc["9qzl9t7l84"]

# + {"jupyter": {"source_hidden": true}}
# df_m = df_m[df_m["stoich"] == "AB3"]
# df_m = df_m[df_m["source"] == "chris"]

# # if row_i["energy_pa"]

# # row_i = df_m.loc["v4zonyzw7d"]
# row_i = df_oqmd_data.iloc[0]

# energy_pa = row_i["energy_pa"]
# if not np.isnan(energy_pa) and row_i["source"] == "oqmd":
#     energy_norm_i = energy_pa

# # row_i

# def method(row_i):
#     """
#     """
#     atoms_i = row_i["atoms"]    
#     num_atoms_i = atoms_i.get_number_of_atoms()
#     return(num_atoms_i)

# df_i = df_m
# df_i["num_atoms"] = df_i.apply(
#     method,
#     axis=1,
#     )
# df_m = df_i