# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.2.1
#   kernelspec:
#     display_name: Python [conda env:PROJ_IrOx_Active_Learning_OER]
#     language: python
#     name: conda-env-PROJ_IrOx_Active_Learning_OER-py
# ---

# # Parse Atoms Objects for IrO2 and IrO3 Unique Prototypes
# ---

# %%capture
# %load_ext autoreload
# %autoreload 2

# +
import os
import sys

import pickle


import numpy as np

from ase import io
from ase.visualize import view
import pandas as pd


import bulk_enumerator as be
import time

from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.ase import AseAtomsAdaptor

# pd.set_option('display.max_rows', None)
# -

# # Reading Structures

# +
root_path = os.path.join(
    os.environ["PROJ_DATA"],
    "04_IrOx_surfaces_OER/ml_bulk_static_preopt_structures")


master_list = []
for root, dirs, files in os.walk(root_path):
    if ".ipynb_checkpoints" in root:
        continue

    if "iro2" in root:
        stoich_i = "AB2"
    elif "iro3" in root:
        stoich_i = "AB3"
    else:
        stoich_i = None

    if "oqmd" in root:
        source_i = "oqmd"
    else:
        source_i = "chris"

    for file_i in files:
        if ".POSCAR" in file_i or ".cif" in file_i:
            id_i = file_i.split("_")[0]

            path_i = root.replace("/mnt/c/Users/raulf/Dropbox/01_norskov/00_projects/", "")

            atoms_i = io.read(
                os.path.join(root, file_i))

            sys_i = {
                "id_old": int(id_i),
                "atoms": atoms_i,
                "stoich": stoich_i,
                "path": path_i,
                "source": source_i,
                }
            master_list.append(sys_i)

df_struct = pd.DataFrame(master_list)
# -

# # Setting Unique ID Tag

# +
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "data/ml_irox_data",
    "unique_ids.csv")
df_id = pd.read_csv(path_i)


id_mapp_iro2 = dict(zip(
    df_id[df_id["stoich"] == "AB2"]["id"],
    df_id[df_id["stoich"] == "AB2"]["unique_ids"]))

id_mapp_iro3 = dict(zip(
    df_id[df_id["stoich"] == "AB3"]["id"],
    df_id[df_id["stoich"] == "AB3"]["unique_ids"]))


# +
def method(row_i):
    id_i = row_i["id_old"]

    if row_i["stoich"] == "AB2":
        unique_id_i = id_mapp_iro2[id_i]
    elif row_i["stoich"] == "AB3":
        unique_id_i = id_mapp_iro3[id_i]
    else:
        print("BADDDDD!!!!! fsdfjisajids")
        unique_id_i = None

    return(unique_id_i)

df_struct["id_unique"] = df_struct.apply(
    method,
    axis=1,
    )

df_struct.set_index("id_unique", inplace=True)
# -

# # Adding secondary index row that is unique and separate from the regular id_unique

# +
df_static_unique_ids = pd.read_csv("static_unique_ids.csv")
id_mapp_static_unique = dict(zip(
    df_static_unique_ids["unique_ids"],
    df_static_unique_ids["static_unique_ids"]))

def method(row_i):
    id_i = row_i.name
    static_id_i = id_mapp_static_unique[id_i]
    return(static_id_i)

df_struct["static_id"] = df_struct.apply(
    method,
    axis=1,
    )

# df_struct.set_index("id_unique", inplace=True)
# -

df_struct.head()

# # Analyzing Structures with Bulk Enumerator

# +
# t0 = time.time()

# data_list = []
# for id_i, row_i in df_struct.iterrows():   
#     atoms_i = row_i["atoms"]

#     structure_i = AseAtomsAdaptor.get_structure(atoms_i)
#     poscar_str_i = Poscar(structure_i).get_string()

#     b = be.bulk.BULK()
#     b.set_structure_from_file(poscar_str_i)

#     spacegroup_i = b.get_spacegroup()
#     species_i = b.get_species()
#     wyckoff_i = b.get_wyckoff()
#     name_i = b.get_name()
#     parameter_values_i = b.get_parameter_values()

#     row_dict_i = {
#         "id": id_i,
#         "spacegroup_i": spacegroup_i,
#         "species_i": species_i,
#         "wyckoff_i": wyckoff_i,
#         "name_i": name_i,
#         "parameter_values_i": parameter_values_i,
#         }
#     data_list.append(row_dict_i)


# t1 = time.time()
# print("time to complete for loop: ")
# print(t1 - t0)

# df_proto = pd.DataFrame(data_list)
# df_proto.set_index("id", inplace=True)

# print(
#     "Number of entries processed: ",
#     len(df_proto["name_i"].to_list())
#     )

# print(len(
#     "Unique entries (some systems with the same prototype): ", 
#     set(df_proto["name_i"].to_list())
#     ))
# -

# # Save data to pickle

directory = "out_data"
if not os.path.exists(directory):
    os.makedirs(directory)

# +
with open("out_data/data_structures.pickle", "wb") as fle:
    pickle.dump(df_struct, fle)

# with open("out_data/data_prototypes.pickle", "wb") as fle:
#     pickle.dump(df_proto, fle)

# + {"active": ""}
#
#
#
#
# -

# # TEMP | Number of atoms in structures

# +
def method(row_i, argument_0, optional_arg=None):
    new_column_values_dict = {"num_atoms": None}

    new_column_values_dict["num_atoms"] = row_i["atoms"].get_number_of_atoms()

    # #########################################################################
    for key, value in new_column_values_dict.items():
        row_i[key] = value
    return(row_i)

df_i = df_struct

arg1 = "TEMP_0"
df_i = df_i.apply(
    method,
    axis=1,
    args=(arg1, ),
    optional_arg="TEMP_1"
    )
df_struct = df_i

df_struct_ab3 = df_struct[df_struct["stoich"] == "AB3"]

df_struct_ab3

print(df_struct_ab3.shape)
print(df_struct_ab3[df_struct_ab3["num_atoms"] > 100].shape)

df_struct_ab3[df_struct_ab3["num_atoms"] > 100]
