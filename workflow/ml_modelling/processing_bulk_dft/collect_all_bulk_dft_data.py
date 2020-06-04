# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_irox]
#     language: python
#     name: conda-env-PROJ_irox-py
# ---

# # Import Modules

# +
import os
print(os.getcwd())
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

# -

# # Read Data

# +
# path_i = os.path.join(
#     os.environ["PROJ_irox"],
#     "workflow/ml_modelling",
#     "processing_bulk_dft/parse_chris_bulk_dft/out_data",
#     "df_dft_calcs.pickle")
# with open(path_i, "rb") as fle:
#     df_chris = pickle.load(fle)
# df_chris["source"] = "chris"

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
    # df_chris,
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

# + active=""
#
#
#
#
#
#
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
# -

# # Removing rows with missing atoms objects

df_m = df_m[df_m["atoms"].notnull()]


# # Count number of atoms

# +
def method(row_i):
    atoms_i = row_i["atoms"]
    num_atoms_i = atoms_i.get_number_of_atoms()
    return(num_atoms_i)

df_m["num_atoms"] = df_m.apply(
    method,
    axis=1)


# +
def method(row_i):
    atoms_i = row_i["atoms"]
    volume = atoms_i.get_volume()

    return(volume)

df_m["volume"] = df_m.apply(
    method,
    axis=1)
# -

df_m["volume_pa"] = df_m.volume / df_m.num_atoms

# # Save data

# +
directory = "out_data"
if not os.path.exists(directory):
    os.makedirs(directory)

with open(os.path.join(directory, "df_bulk_dft.pickle"), "wb") as fle:
    pickle.dump(df_m, fle)
# -
print(20 * "# # ")
print("All done!")
assert False
