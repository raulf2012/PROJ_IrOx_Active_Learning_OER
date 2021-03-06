# ---
# jupyter:
#   jupytext:
#     cell_metadata_json: true
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

# + [markdown] {"Collapsed": "false"}
# # Parse Atoms Objects for IrO2 and IrO3 Unique Prototypes
# ---

# + {"Collapsed": "false"}
# %%capture
# %load_ext autoreload
# %autoreload 2

# + {"Collapsed": "false"}
import os
print(os.getcwd())
import sys

import time

import pickle

import numpy as np
import pandas as pd

from ase import io
from ase.visualize import view

from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.ase import AseAtomsAdaptor

# #############################################################################
from IPython.display import display

# + [markdown] {"Collapsed": "false"}
# # Reading Structures

# + {"Collapsed": "false"}
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

            path_i = root

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

# + [markdown] {"Collapsed": "false"}
# # Setting Unique ID Tag

# + {"Collapsed": "false"}
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


# + {"Collapsed": "false"}
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

# + [markdown] {"Collapsed": "false"}
# # Adding secondary index row that is unique and separate from the regular id_unique

# + {"Collapsed": "false"}
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

# + [markdown] {"Collapsed": "false"}
# # Analyzing Structures with Bulk Enumerator
# -

is_bulk_enumerator_installed = False
try:
    import bulk_enumerator as be
    is_bulk_enumerator_installed = True
except:
    print("bulk_enumerator is not installed/importable")
    print("Contact ankitjain.me.iitk@gmail.com to be added as a guest so that you can install the Enumerator package")
    print("https://gitlab.com/ankitjainmeiitk/Enumerator")

# + {"Collapsed": "false"}
if is_bulk_enumerator_installed:
    t0 = time.time()

    data_list = []
    for id_i, row_i in df_struct.iterrows():
        atoms_i = row_i["atoms"]

        structure_i = AseAtomsAdaptor.get_structure(atoms_i)
        poscar_str_i = Poscar(structure_i).get_string()

        b = be.bulk.BULK()
        b.set_structure_from_file(poscar_str_i)

        spacegroup_i = b.get_spacegroup()
        species_i = b.get_species()
        wyckoff_i = b.get_wyckoff()
        name_i = b.get_name()
        parameter_values_i = b.get_parameter_values()

        row_dict_i = {
            "id": id_i,
            "spacegroup_i": spacegroup_i,
            "species_i": species_i,
            "wyckoff_i": wyckoff_i,
            "name_i": name_i,
            "parameter_values_i": parameter_values_i,
            }
        data_list.append(row_dict_i)


    t1 = time.time()
    print("time to complete for loop: ")
    print(t1 - t0)

    df_proto = pd.DataFrame(data_list)
    df_proto.set_index("id", inplace=True)

    print(
        "Number of entries processed: ",
        len(df_proto["name_i"].to_list())
        )

    print(len(
        "Unique entries (some systems with the same prototype): ", 
        set(df_proto["name_i"].to_list())
        ))

# + [markdown] {"Collapsed": "false"}
# # Save data to pickle

# + {"Collapsed": "false"}
directory = "out_data"
if not os.path.exists(directory):
    os.makedirs(directory)

# + {"Collapsed": "false"}
if True:
# if False:
    with open("out_data/data_structures.pickle", "wb") as fle:
        pickle.dump(df_struct, fle)

    if is_bulk_enumerator_installed:
        with open("out_data/data_prototypes.pickle", "wb") as fle:
            pickle.dump(df_proto, fle)

# + {"Collapsed": "false"}
# #############################################################################
path_i = os.path.join(
    "out_data",
    "data_structures.pickle")
with open(path_i, "rb") as fle:
    df_struct = pickle.load(fle)
# #############################################################################

# #############################################################################
if is_bulk_enumerator_installed:
    path_i = os.path.join(
        "out_data",
        "data_prototypes.pickle")
    with open(path_i, "rb") as fle:
        df_proto = pickle.load(fle)
else:
    # COMBAK Read from PROJ_DATA instead
    df_proto = None

    path_i = os.path.join(
        os.environ["PROJ_DATA"],
        "04_IrOx_surfaces_OER/PROJECT_COMPUTED_OUT_DATA/PROJ_IrOx_Active_Learning_OER",
        "workflow/ml_modelling/processing_bulk_dft/static_prototypes_structures",
        "out_data/data_prototypes.pickle")
    with open(path_i, "rb") as fle:
        df_proto = pickle.load(fle)
# #############################################################################

# + {"active": ""}
#
#
#
#
# -

print(20 * "# # ")
print("All done!")
assert False
