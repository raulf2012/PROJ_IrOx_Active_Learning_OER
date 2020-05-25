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

# +
import os
print(os.getcwd())
import sys

import ase.db
import pandas as pd

# +
db_path = os.path.join(
    os.environ["PROJ_DATA"],
    "04_IrOx_surfaces_OER/ml_bulk_static_preopt_structures/reoptimized_structures_from_kirsten",
    "OptimizedStructures.db")

db = ase.db.connect(db_path)
# -

data_list = []
for row in db.select():
    atoms_i = row.toatoms()
    id_unique_i = row.structure_id
    
    data_dict_i = {
        "atoms": atoms_i,
        "id_unique": id_unique_i}
    data_list.append(data_dict_i)

df_struct = pd.DataFrame(data_list)
df_struct = df_struct.set_index("id_unique")

# Pickling data ######################################################
import os; import pickle
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)
with open(os.path.join(directory, "data_structures_kirsten.pickle"), "wb") as fle:
    pickle.dump(df_struct, fle)
# #####################################################################

from ase_modules.ase_methods import view_in_vesta

# +
df_struct = df_struct.sample(n=10)

atoms = df_struct["atoms"].tolist()
names_list = df_struct.index.tolist()

# view_in_vesta(atoms, ase_gui=False, name_list=names_list)
