# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_irox] *
#     language: python
#     name: conda-env-PROJ_irox-py
# ---

# # Import Modules

# +
import os

import pandas as pd

from ase import io

os.getcwd()
# -

# # Read Atoms Objects

# +
rootdir = os.path.join(
    os.environ["PROJ_DATA"],
    "04_IrOx_surfaces_OER/bulk_systems")

bulk_atoms_dict = {
    "IrO3_rutile-like": os.path.join(rootdir, "r_iro3/_1/OUTCAR"),
    "IrO3": os.path.join(rootdir, "a_iro3/_2/OUTCAR"),
    "IrO3_battery": os.path.join(rootdir, "b_iro3/_3/OUTCAR"),
    "IrO2": os.path.join(rootdir, "r_iro2/_1/OUTCAR"),
    }
# -

data_list = []
for key, value in bulk_atoms_dict.items():

    if "IrO3" in key:
        stoich_i = "AB3"
    elif "IrO2" in key:
        stoich_i = "AB2"
    else:
        stoich_i = None

    atoms_i = io.read(value)

    data_dict_i = {
        "atoms": atoms_i,
        "id": key,
        "path": value,
        "stoich": stoich_i,
        }

    data_list.append(data_dict_i)

# +
df_oer_bulk = pd.DataFrame(data_list)

df_oer_bulk["source"] = "raul_oer"
df_oer_bulk["id_unique"] = df_oer_bulk["id"]

df_oer_bulk = df_oer_bulk.set_index("id_unique")
# -

# # Save Data

# +
directory = "out_data"
if not os.path.exists(directory):
    os.makedirs(directory)

import pickle
with open(os.path.join(directory, "df_oer_bulk.pickle"), "wb") as fle:
    pickle.dump(df_oer_bulk, fle)
# -

print(20 * "# # ")
print("All done!")
assert False
