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
import os
print(os.getcwd())
import sys

import pandas as pd

from ase.db import connect

from protosearch.build_bulk.cell_parameters import CellParameters
from protosearch.build_bulk.classification import get_classification

# +
# # %%capture

# sys.path.insert(0,
#     os.path.join(
#         os.environ["PROJ_irox"],
#         "workflow/ml_modelling"))
# from ml_methods import get_ml_dataframes, get_data_for_al

# DF_dict = get_data_for_al()
# # list(DF_dict.keys())

# df_features_post = DF_dict["df_features_post"]
# df_features_pre = DF_dict["df_features_pre"]
# df_bulk_dft = DF_dict["df_bulk_dft"]
# df_ids = DF_dict["df_ids"]
# df_dij = DF_dict["df_dij"]
# df_static_irox = DF_dict["df_static_irox"]

# +
# %%capture

sys.path.insert(0, 
    os.path.join(
        os.environ["PROJ_irox"],
        "workflow/ml_modelling"))

from ml_methods import get_ml_dataframes
DF_dict = get_ml_dataframes(
    names=[
        'bulk_dft_data_path',
        'unique_ids_path',
        'prototypes_data_path',
        'static_irox_structures_path',
        'static_irox_structures_kirsten_path',
        'oqmd_irox_data_path',
        'df_features_pre_opt_path',
        'df_features_pre_opt_kirsten_path',
        'df_features_post_opt_path',
        'oer_bulk_structures_path',
        'df_ccf_path',
        'df_dij_path',
        'ids_to_discard__too_many_atoms_path',
        'df_prototype_dft_path',
        'df_prototype_static_path',
        ]
    )

df_bulk_dft = DF_dict["bulk_dft_data"]
df_static_irox = DF_dict["static_irox_structures"]


# bulk_dft_data = DF_dict["bulk_dft_data"]
# unique_ids = DF_dict["unique_ids"]
# prototypes_data = DF_dict["prototypes_data"]
# static_irox_structures = DF_dict["static_irox_structures"]
# static_irox_structures_kirsten = DF_dict["static_irox_structures_kirsten"]
# oqmd_irox_data = DF_dict["oqmd_irox_data"]
# df_features_pre_opt = DF_dict["df_features_pre_opt"]
# df_features_pre_opt_kirsten = DF_dict["df_features_pre_opt_kirsten"]
# df_features_post_opt = DF_dict["df_features_post_opt"]
# oer_bulk_structures = DF_dict["oer_bulk_structures"]
# df_ccf = DF_dict["df_ccf"]
# df_dij = DF_dict["df_dij"]
# ids_to_discard__too_many_atoms = DF_dict["ids_to_discard__too_many_atoms"]
# df_prototype_static = DF_dict["df_prototype_static"]
# df_prototype_dft = DF_dict["df_prototype_dft"]
# -

df_bulk_dft = df_bulk_dft[df_bulk_dft.source == "raul"]


# + active=""
#
#
#
#
# -

def process_row(atoms, id_i):
    """
    """
    # #####################################################
    out_dict = dict()

    out_dict["id_unique"] = id_i

    # #####################################################
    prototype = get_classification(atoms)

    p_name = prototype[0]["p_name"]
    structure_name = prototype[0]["structure_name"]
    spacegroup = prototype[0]["spacegroup"]
    wyckoffs = prototype[0]["wyckoffs"]
    species = prototype[0]["species"]

    # #####################################################    
    out_dict["p_name"] = p_name
    out_dict["structure_name"] = structure_name
    out_dict["spacegroup"] = spacegroup
    out_dict["wyckoffs"] = wyckoffs
    out_dict["species"] = species

    # #####################################################    
    wyckoff_params = prototype[1]


    # #####################################################
    return(out_dict)

# # DFT Opt. Structures Processing

# +
# # df_bulk_dft = df_bulk_dft.sample(n=20)

# data_dict_list = []
# for id_i, row_i in df_bulk_dft.iterrows():
#     atoms = row_i.atoms

#     out_dict = process_row(atoms, id_i)
#     data_dict_list.append(out_dict)

# df_prototype_dft = pd.DataFrame(data_dict_list)
# df_prototype_dft = df_prototype_dft.set_index("id_unique")

# df_prototype_dft.head()
# -

# # Static Structures Processing

assert False

df_static_irox.loc["8p8evt9pcg"]

# +
# df_static_irox = df_static_irox.sample(n=20)

data_dict_list = []
for id_i, row_i in df_static_irox.iterrows():
    atoms = row_i.atoms

    out_dict = process_row(atoms, id_i)
    data_dict_list.append(out_dict)

df_prototype_static = pd.DataFrame(data_dict_list)
df_prototype_static = df_prototype_static.set_index("id_unique")

df_prototype_static.head()
# -

# df_prototype_static.loc["8p8evt9pcg"]
df_prototype_static.shape

# + active=""
#
#
#
#

# +
# import os; import pickle
# directory = "out_data"
# if not os.path.exists(directory): os.makedirs(directory)

# # Pickling data ###########################################
# with open(os.path.join(directory, "df_prototype_static.pickle"), "wb") as fle:
#     pickle.dump(df_prototype_static, fle)

# with open(os.path.join(directory, "df_prototype_dft.pickle"), "wb") as fle:
#     pickle.dump(df_prototype_dft, fle)
# #########################################################

# + active=""
#
#
#
#

# + jupyter={}
# df_static_irox.

# data_dict_list = []

# for id_i, row_i in df_bulk_dft.iterrows():
#     atoms = row_i.atoms

#     # #####################################################
#     out_dict = dict()

#     out_dict["id_unique"] = id_i

#     # #####################################################
#     prototype = get_classification(atoms)

#     p_name = prototype[0]["p_name"]
#     structure_name = prototype[0]["structure_name"]
#     spacegroup = prototype[0]["spacegroup"]
#     wyckoffs = prototype[0]["wyckoffs"]
#     species = prototype[0]["species"]

#     # #####################################################    
#     out_dict["p_name"] = p_name
#     out_dict["structure_name"] = structure_name
#     out_dict["spacegroup"] = spacegroup
#     out_dict["wyckoffs"] = wyckoffs
#     out_dict["species"] = species

#     # #####################################################    
#     wyckoff_params = prototype[1]


#     # #####################################################        
#     data_dict_list.append(out_dict)

# df_prototype = pd.DataFrame(data_dict_list)
# df_prototype = df_prototype.set_index("id_unique")
