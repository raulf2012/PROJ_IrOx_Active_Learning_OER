# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_IrOx_Active_Learning_OER]
#     language: python
#     name: conda-env-PROJ_IrOx_Active_Learning_OER-py
# ---

# # Analysing prototypes of structures computed for OER section of paper
# ---

# +
import os
import sys

import pickle

# ase
from ase.visualize import view
from ase import io

# pymatgen
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.ase import AseAtomsAdaptor


import bulk_enumerator as be

# + [markdown] toc-hr-collapsed=true
# # Analysing Bulk Structures
#
# 01_IrO2  02_IrO3  03_iro3_rutile-like  04_iro3_battery
# -

prototype_names = {}

# ## rutile-IrO2

# +
file_path_i = os.path.join(
    os.environ["PROJ_irox"],
    "01_bulk_structures",
    "01_IrO2",
    "bulk_symm.cif",
    )

atoms_i = io.read(file_path_i)
structure_i = AseAtomsAdaptor.get_structure(atoms_i)

b = be.bulk.BULK()
b.set_structure_from_file(Poscar(structure_i).get_string())

prototype_names["rutile_iro2"] = b.get_name()
print(prototype_names["rutile_iro2"])
# -

# ## a-AlF3 IrO3

# +
file_path_i = os.path.join(
    os.environ["PROJ_irox"],
    "01_bulk_structures",
    "02_IrO3",
    "bulk.cif",
    )

atoms_i = io.read(file_path_i)
structure_i = AseAtomsAdaptor.get_structure(atoms_i)

b = be.bulk.BULK()
b.set_structure_from_file(Poscar(structure_i).get_string())

prototype_names["a-AlF3_iro3"] = b.get_name()
print(prototype_names["a-AlF3_iro3"])
# -

# ## rutile-like IrO3

# +
file_path_i = os.path.join(
    os.environ["PROJ_irox"],
    "01_bulk_structures",
    "03_iro3_rutile-like",
    "bulk.cif",
    )

atoms_i = io.read(file_path_i)
structure_i = AseAtomsAdaptor.get_structure(atoms_i)

b = be.bulk.BULK()
b.set_structure_from_file(Poscar(structure_i).get_string())

prototype_names["rutile_iro3"] = b.get_name()
print(prototype_names["rutile_iro3"])
# -

# ## battery IrO3

# +
file_path_i = os.path.join(
    os.environ["PROJ_irox"],
    "01_bulk_structures",
    "04_iro3_battery",
    "bulk.cif",
    )

atoms_i = io.read(file_path_i)
structure_i = AseAtomsAdaptor.get_structure(atoms_i)

b = be.bulk.BULK()
b.set_structure_from_file(Poscar(structure_i).get_string())

prototype_names["battery_iro3"] = b.get_name()
print(prototype_names["battery_iro3"])
# -

# # Results

# + active=""
# rutile IrO2: AB2_2_a_f_136
# a-AlF3 IrO3: AB3_2_b_e_167
# rutile IrO3: AB3_4_e_fj_136
# batter IrO3: AB3_8_ik_glm2_66
# -

# # Comparing to Static Prototypes List

# +
# file_i = os.path.join(
#     os.environ["PROJ_irox"],
#     "workflow/ml_modelling/00_ml_workflow/outdata",
#     "01_irox_data_featurized.pickle")
# with open(file_i, "rb") as fle:
#     df_m = pickle.load(fle)

# +
import pickle
import os


sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import prototypes_data_path, unique_ids_path

with open(prototypes_data_path, "rb") as fle:
    df_proto = pickle.load(fle)

import pandas as pd
df_ids = pd.read_csv(unique_ids_path)


df_ids_ab2 = df_ids[df_ids["stoich"] == "AB2"]
df_ids_ab3 = df_ids[df_ids["stoich"] == "AB3"]
# -

df_iro2 = df_proto.loc[df_ids_ab2["unique_ids"]]
df_iro3 = df_proto.loc[df_ids_ab3["unique_ids"]]

# +
# df_iro2 = df_m[df_m["default_columns"]["stoich"] == "AB2"]

# df_iro3 = df_m[df_m["default_columns"]["stoich"] == "AB3"]

# +
prototype_names

#IrO2: 200
#a-AlF3 IrO3: 231
#rutile IrO3: 164
# -

# ## rutile IrO2

# +
name_i = prototype_names["rutile_iro2"]
name_i in df_iro2["name_i"].tolist()

df_iro2.loc[df_iro2["name_i"] == name_i]
# -

# ## a-AlF3 IrO3

# +
name_i = prototype_names["a-AlF3_iro3"]
name_i in df_iro3["name_i"].tolist()

df_iro3.loc[df_iro3["name_i"] == name_i]
# -

# ## rutile IrO3

# +
name_i = prototype_names["rutile_iro3"]
name_i in df_iro3["name_i"].tolist()

df_iro3.loc[df_iro3["name_i"] == name_i]
# -

# ## battery IrO3

# +
name_i = prototype_names["battery_iro3"]
name_i in df_iro3["name_i"].tolist()

df_iro3.loc[df_iro3["name_i"] == name_i]
