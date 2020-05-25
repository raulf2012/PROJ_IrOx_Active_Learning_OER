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
print(os.getcwd())
import sys

import pandas as pd

from ase.db import connect
# -

# # Script Inputs

filename = "FinalStructures_1.db"

# # Read df_bulk_dft dataframe

# +
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"], "workflow/ml_modelling"))
from ml_methods import get_ml_dataframes

DF_dict = get_ml_dataframes(names=[
    "df_dft_final_final_path",
    # "",
    ])

df_bulk_dft = DF_dict["df_dft_final_final"]

df_i = df_bulk_dft

# +
# df_i = df_i.loc[[
#     'vovgximhm2',
#     '8dce6kz2vf',
#     'vhv39q6e9j',
#     '8ymh8qnl6o',
#     '6fcdbh9fz2',
#     '7qm56wxj8s',
#     'mu6omk6k9l',
#     '6dzhcimdxs',
#     ]]

# +
# assert False
# -

# # WARNING: DON"T RUN THIS AGAIN BECAUSE IT WILL APPEND TO DATABASE
#
# I will comment out this code block to avoid issues, uncomment and run only if you're starting from scratch

# +
# # db = connect('out_data/FinalStructures_1.db')
# db = connect(os.path.join("out_data", filename))

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

# # Sandbox

# +
# db = connect('out_data/FinalStructures_1.db')
db = connect(os.path.join("out_data", filename))

tmp_list = []
for row in db.select():
    tmp = 42   
    id_i = row.key_value_pairs["structure_id"]
    tmp_list.append(id_i)

print("Unique number of ids in database:", len(set(tmp_list)))

# +
id_list = [
    'vovgximhm2',
    '8dce6kz2vf',
    'vhv39q6e9j',
    '8ymh8qnl6o',
    '6fcdbh9fz2',
    '7qm56wxj8s',
    'mu6omk6k9l',
    '6dzhcimdxs',
    ]

for id_i in id_list:
    print(id_i in tmp_list)

# + active=""
#
#
#

# + jupyter={}
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

# + jupyter={}
# df_i.shape

# df_bulk_dft_ab2.shape

# df_bulk_dft_ab2.source.unique()

# + jupyter={}
# # #############################################################################
# import pickle; import os
# path_i = os.path.join(
#      os.environ["PROJ_irox"],
#      "workflow/ml_modelling/processing_bulk_dft/static_prototypes_structures/out_data",
#      "ids_to_discard__proto_dupl.pickle")
# with open(path_i, "rb") as fle:
#     duplicate_ids = pickle.load(fle)
# # #############################################################################

# len(duplicate_ids)

# +
# # %%capture

# sys.path.insert(0, os.path.join(
#     os.environ["PROJ_irox"],    
#     "workflow/ml_modelling"))

# # #####################################
# data_dict = get_data_for_al(
#     stoich="AB2",
#     verbose=True,
#     drop_too_many_atoms=True)
# df_bulk_dft = data_dict["df_bulk_dft"]
# df_bulk_dft_ab2 = df_bulk_dft
# df_bulk_dft_ab2 = df_bulk_dft

# # #####################################
# data_dict = get_data_for_al(
#     stoich="AB3",
#     verbose=True,
#     drop_too_many_atoms=True)
# df_bulk_dft = data_dict["df_bulk_dft"]
# df_bulk_dft_ab3 = df_bulk_dft


# # #############################################################################
# df_bulk_dft = pd.concat(
#     [
#         df_bulk_dft_ab2,
#         df_bulk_dft_ab3])

# df_bulk_dft = df_bulk_dft[df_bulk_dft.source == "raul"]
# df_i = df_bulk_dft
