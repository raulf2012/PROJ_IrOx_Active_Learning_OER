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

# +
import pickle
import os

path_i = "out_data/data_structures.pickle"
with open(path_i, "rb") as fle:
    df_static = pickle.load(fle)

# +
row_i = df_static.iloc[0]

atoms = row_i["atoms"]

num_atoms = atoms.get_number_of_atoms()


# +
def method(row_i, argument_0, optional_arg=None):
    """
    """    
#     row_i = df_static.iloc[0]

    atoms = row_i["atoms"]

    num_atoms = atoms.get_number_of_atoms()

    return(num_atoms)

arg1 = "TEMP_0"
df_i = df_static
df_i["num_atoms"] = df_i.apply(
    method,
    axis=1,
    args=(arg1, ),
    optional_arg="TEMP_1"
    )
df_static = df_i

# +
df_ab2 = df_static[df_static["stoich"] == "AB2"]

df_ab3 = df_static[df_static["stoich"] == "AB3"]
# -

num_atoms_cutoff = 70
print(df_ab2[df_ab2["num_atoms"] > num_atoms_cutoff].shape)
print(df_ab3[df_ab3["num_atoms"] > num_atoms_cutoff].shape)
