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

import yaml

import pandas as pd

from pymatgen.ext.matproj import MPRester
from pymatgen.io.ase import AseAtomsAdaptor
# -

# # Read Data

# +
path_i = os.path.join(os.environ["PROJ_irox"], "config", "config.yml")
with open(path_i) as file:
    config_dict = yaml.load(file, Loader=yaml.FullLoader)

api_key = config_dict['materials_project']['api_key']

MPR = MPRester(
    api_key=api_key,
    endpoint=None,
    include_user_agent=True,
    )

# +
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling"))
from ml_methods import get_ml_dataframes


DF_dict = get_ml_dataframes()

df_dij = DF_dict['df_dij']
df_dft_final_final = DF_dict['df_dft_final_final']
df_dft = df_dft_final_final

# +
# df_dft = df_dft.iloc[0:20]

df_dft = df_dft.loc[[
    'cg8p7fxq65',
    '64cg6j9any',
    '85z4msnl6o',
    'xozr8f7p7g',
    '949rnem5z2',
    'mkmsvkcyc5',
    'vwxfn3blxi',
    'nrml6dms9l',
    ]]

# +
tmp_list = []

data_dict_list = []
for i_cnt, row_i in df_dft.iterrows():
    data_dict_i = dict()

    # #####################################################
    name_i = row_i.name
    stoich_i = row_i.stoich

    # #####################################################
    data_dict_i["id"] = name_i
    data_dict_i["stoich"] = stoich_i

    atoms_i = row_i.atoms
    struct_i = AseAtomsAdaptor.get_structure(atoms_i)

    duplicates_tmp = MPR.find_structure(struct_i)
    tmp_list.append(duplicates_tmp)

    data_dict_i["mp_duplicates"] = duplicates_tmp
    
    data_dict_list.append(data_dict_i)

# +
df = pd.DataFrame(data_dict_list)

df_mp_dupl = df[[True if len(i) != 0 else False for i in df.mp_duplicates.tolist()]]
df_mp_dupl = df_mp_dupl.set_index("id")

df_mp_dupl
# -

# Pickling data ###########################################
import os; import pickle
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)
with open(os.path.join(directory, "df_mp_dupl.pickle"), "wb") as fle:
    pickle.dump(df_mp_dupl, fle)
# #########################################################

# #########################################################
import pickle; import os
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "CatHub_MPContribs_upload/MPContribs_upload/duplicate_MP_entries",
    "out_data/df_mp_dupl.pickle")
with open(path_i, "rb") as fle:
    df_mp_dupl = pickle.load(fle)
# #########################################################

# + active=""
#
#
#

# +
# bulk_dft_data = DF_dict['bulk_dft_data']
# unique_ids = DF_dict['unique_ids']
# prototypes_data = DF_dict['prototypes_data']
# static_irox_structures = DF_dict['static_irox_structures']
# static_irox_structures_kirsten = DF_dict['static_irox_structures_kirsten']
# oqmd_irox_data = DF_dict['oqmd_irox_data']
# df_features_pre_opt = DF_dict['df_features_pre_opt']
# df_features_pre_opt_kirsten = DF_dict['df_features_pre_opt_kirsten']
# df_features_post_opt = DF_dict['df_features_post_opt']
# oer_bulk_structures = DF_dict['oer_bulk_structures']
# df_ccf = DF_dict['df_ccf']
# ids_to_discard__too_many_atoms = DF_dict['ids_to_discard__too_many_atoms']
# ids_duplicates = DF_dict['ids_duplicates']

# + jupyter={}
# df_dij.loc[
#     [
#         "mkmsvkcyc5",
#         "xozr8f7p7g",
#         ],
#     [
#         "mkmsvkcyc5",
#         "xozr8f7p7g",
#         ]
#     ]

# "xozr8f7p7g" in df_dft.index
# "mkmsvkcyc5" in df_dft.index

# df_dft.loc["mkmsvkcyc5"]

# df_dft.loc["xozr8f7p7g"]

# + jupyter={}
# assert False

# + jupyter={}
# # df_w_duplicates.to_csv()

# # df_w_duplicates.id.tolist()
# # df_w_duplicates.mp_duplicates.tolist()

# df_w_duplicates
