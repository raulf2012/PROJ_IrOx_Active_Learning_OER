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

# # Import Modules

# +
import os
print(os.getcwd())

import sys

import pickle

import pandas as pd

import bulk_enumerator as be
import time

from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.ase import AseAtomsAdaptor

# #############################################################################
sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))     

from proj_data_irox import (
    bulk_dft_data_path,
    unique_ids_path,
    prototypes_data_path,
    static_irox_structures_path,
    oqmd_irox_data_path,
    )
# -

# # Read Data

# +
# #############################################################################
with open(bulk_dft_data_path, "rb") as fle:
    df_bulk_dft = pickle.load(fle)
# #############################################################################

df_bulk_dft = df_bulk_dft[df_bulk_dft.source == "raul"]
# -

# # Classify prototype info

# tolerance = 1e-12
# tolerance = 1e-9
# tolerance = 1e-8
# tolerance = 1e-7
# tolerance = 1e-6
# tolerance = 1e-5
tolerance = 1e-4
# tolerance = 1e-3
# tolerance = 1e-2
# tolerance = 1e-1


# + jupyter={}
t0 = time.time()

data_list = []
# for id_i, row_i in df_bulk_dft.iloc[0:20].iterrows():
for id_i, row_i in df_bulk_dft.iterrows():

    atoms_i = row_i["atoms"]

    structure_i = AseAtomsAdaptor.get_structure(atoms_i)
    poscar_str_i = Poscar(structure_i).get_string()

    b = be.bulk.BULK(
        tolerance=tolerance, 
        )
    b.set_structure_from_file(poscar_str_i)

    spacegroup_i = b.get_spacegroup()
    species_i = b.get_species()
    wyckoff_i = b.get_wyckoff()
    name_i = b.get_name()
    parameter_values_i = b.get_parameter_values()

    primitive_natom = b.get_primitive_natom()
    std_natom = b.get_std_natom()

    row_dict_i = {
        "id": id_i,
        "spacegroup_i": spacegroup_i,
        "species_i": species_i,
        "wyckoff_i": wyckoff_i,
        "name_i": name_i,
        "parameter_values_i": parameter_values_i,
        "primitive_natoms": primitive_natom,
        "std_natom": std_natom,
        }
    data_list.append(row_dict_i)


t1 = time.time()
print("time to complete for loop: ", t1 - t0, "sec")
print("time to complete for loop (per iter): ", (t1 - t0) / len(data_list), "sec")
print("")

df_proto = pd.DataFrame(data_list)
df_proto.set_index("id", inplace=True)

print(
    "Number of entries processed: ",
    len(df_proto["name_i"].to_list())
    )

print(
    "Unique entries (some systems with the same prototype): ", 
    len(set(df_proto["name_i"].tolist())),
    )


# +
num_atoms_removed = (df_proto.std_natom - df_proto.primitive_natoms).sum()

print("num_atoms_removed:", num_atoms_removed)

# + jupyter={"outputs_hidden": true}
df_proto

# + active=""
# 1e-1 | 0
# 1e-2 | 105
# 1e-3 | 126
# 1e-4 | 114
# 1e-5 | 117
# 1e-6 | 114
# 1e-7 | 114
# 1e-8 | 114
# 1e-9 | 114

# + active=""
# 1e-3 | 5784
# 1e-4 | 5323

# + active=""
#
#
#
#

# + jupyter={}
# b.get_name()
# b.get_parameter_gradients()
# b.get_parameter_values()
# b.get_parameters()
# b.get_primitive_natom()
# b.get_primitive_poscar()
# b.get_spacegroup()
# b.get_species()
# b.get_species_permutations()
# b.get_std_natom()
# b.get_std_poscar()
# b.get_wyckoff()
# b.get_wyckoff_list()
# b.get_wyckoff_structure_map()

# + jupyter={}
# b.get_primitive_natom()
# b.get_std_natom()
# b.get_std_poscar()
# b.get_wyckoff()
