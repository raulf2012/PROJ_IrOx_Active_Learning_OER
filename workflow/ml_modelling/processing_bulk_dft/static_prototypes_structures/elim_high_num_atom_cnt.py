# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_irox] *
#     language: python
#     name: conda-env-PROJ_irox-py
# ---

# # Import Modules

# +
import os
print(os.getcwd())
import sys
import pickle

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import (
    bulk_dft_data_path,
    static_irox_structures_path)
# -

# # Script Inputs

atom_cutoff = 75

# # Read Data

# +
with open(bulk_dft_data_path, "rb") as fle:
    df_bulk_dft = pickle.load(fle)

with open(static_irox_structures_path, "rb") as fle:
    df_static_irox = pickle.load(fle)


# -

# # Add 'number of atoms' column to dataframe

# +
def method(row_i):
    atoms_i = row_i["atoms"]
    num_atoms_i = atoms_i.get_number_of_atoms()
    return(num_atoms_i)

df_i = df_static_irox
df_i["num_atoms"] = df_i.apply(method, axis=1)
df_static_irox = df_i
# -

df_discard = df_static_irox[df_static_irox["num_atoms"] > atom_cutoff]
ids_to_discard = df_discard.index.unique().tolist()

print(
    "Number of IrO2 structures that are above the atom_number cutoff",
    "\n",
    df_discard[df_discard["stoich"] == "AB2"].shape[0],
    )
print("")
print(
    "Number of IrO3 structures that are above the atom_number cutoff",
    "\n",
    df_discard[df_discard["stoich"] == "AB3"].shape[0],
    )

698 - 131

# # Saving data

directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)
with open(os.path.join(directory, "ids_to_discard__too_many_atoms.pickle"), "wb") as fle:
    pickle.dump(ids_to_discard, fle)

# + active=""
#
#
#
#

# + jupyter={}
# df_static_iro2 = df_static_irox[df_static_irox["stoich"] == "AB2"]

# df_static_iro3 = df_static_irox[df_static_irox["stoich"] == "AB3"]

# + jupyter={}
# atom_cutoff = 100

# print(df_static_iro3[df_static_iro3["num_atoms"] > atom_cutoff].shape)
# print(df_static_iro2[df_static_iro2["num_atoms"] > atom_cutoff].shape)

# df_iro3_discard = df_static_iro3[df_static_iro3["num_atoms"] > atom_cutoff]
# df_iro2_discard = df_static_iro2[df_static_iro2["num_atoms"] > atom_cutoff]

# df_iro2_discard

# + jupyter={}
# row_i = df_static_irox.iloc[0]

# atoms_i = row_i["atoms"]

# atoms_i.get_number_of_atoms()
