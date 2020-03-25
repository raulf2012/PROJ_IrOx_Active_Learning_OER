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

# # Import Modules

# +
import os
import sys

import pandas as pd

from ase.db import connect

# import ase
# ase.__version__
# -

# # Read df_bulk_dft dataframe

# +
# %%capture

sys.path.insert(0,
    os.path.join(
        os.environ["PROJ_irox"],    
        "workflow/ml_modelling"))

from ml_methods import get_data_for_al

# #####################################
data_dict = get_data_for_al(
    stoich="AB2",
    verbose=True,
    drop_too_many_atoms=True)
df_bulk_dft = data_dict["df_bulk_dft"]
df_bulk_dft_ab2 = df_bulk_dft

# #####################################
data_dict = get_data_for_al(
    stoich="AB3",
    verbose=True,
    drop_too_many_atoms=True)
df_bulk_dft = data_dict["df_bulk_dft"]
df_bulk_dft_ab3 = df_bulk_dft


# #############################################################################
df_bulk_dft = pd.concat(
    [
        df_bulk_dft_ab2,
        df_bulk_dft_ab3])

df_bulk_dft = df_bulk_dft[df_bulk_dft.source == "raul"]
df_i = df_bulk_dft
#  path_i = os.path.join(
#      os.environ["PROJ_irox"],
#      "workflow/ml_modelling/processing_bulk_dft/static_prototypes_structures/out_data",
#      "ids_to_discard__proto_dupl.pickle")
# TEMP
# -

df_i.shape

assert False

db = connect('out_data/FinalStructures_1.db')

# # WARNING: DON"T RUN THIS AGAIN BECAUSE IT WILL APPEND TO DATABASE
#
# I will comment out this code block to avoid issues, uncomment and run only if you're starting from scratch

# +
# for i_cnt, row_i in df_i.iterrows():
#     atoms = row_i.atoms
#     structure_id = row_i.name

#     key_value_pairs = dict(
#         stoich=row_i.stoich,
#         id_old=row_i.id_old,
#         structure_id=row_i.name,
#         )

#     db.write(
#         atoms,
#         key_value_pairs=key_value_pairs,
#         # data={},
#         # id=,
#         )
# -

assert False

# +
db = connect('out_data/FinalStructures_1.db')

tmp_list = []
for row in db.select():
    tmp = 42   
    id_i = row.key_value_pairs["structure_id"]
    tmp_list.append(id_i)

print("Unique number of ids in database:", len(set(tmp_list)))

# + active=""
#
#
#

# + jupyter={"source_hidden": true}
# # db.write?

# for row in db.select():
#     tmp = 42


# row.key_value_pairs

# # db.write?

# for row in db.select(0):
#     tmp = 42
    
# row

# import pickle
# from pandas.core.internals import managers
# data = pickle.load(open('df_bulk_dft.pickle', 'rb'))

# indices = df_i.index.values
# for i, atoms in enumerate(df_i.atoms.values):
#     db.write(atoms, structure_id=indices[i])
