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

# + [markdown] Collapsed="false"
# # Import Modules

# + Collapsed="false"
import os
print("os.getcwd():", os.getcwd())

import sys

# Add paths to $PATH
# sys.path.insert(0, "/home/raulf2012/Dropbox/01_norskov/00_git_repos/PROJ_IrOx_Active_Learning_OER/workflow/ml_modelling")
sys.path.insert(0, 
    os.path.join(
        os.environ["PROJ_irox"],
        "workflow/ml_modelling"
        )                
    )

# #############################################################################
from ml_methods import get_ml_dataframes

# #############################################################################
from IPython.display import display

# + [markdown] Collapsed="false"
# # Read Data

# + [markdown] Collapsed="false"
# ## Read dataframes

# + Collapsed="false"
DF_dict = get_ml_dataframes(
    names=[
        "bulk_dft_data_path",
        "unique_ids_path",
        "prototypes_data_path",
        "static_irox_structures_path",
        "static_irox_structures_kirsten_path",
        "oqmd_irox_data_path",
        "df_features_pre_opt_path",
        "df_features_pre_opt_kirsten_path",
        "df_features_post_opt_path",
        # "df_features_path",
        # "df_features_cleaned_path",
        # "df_features_cleaned_pca_path",
        "oer_bulk_structures_path",
        "df_ccf_path",
        "df_dij_path",
        "ids_to_discard__too_many_atoms_path",
        ],

    )

# + Collapsed="false"
DF_dict.keys()

df_bulk_dft = DF_dict["bulk_dft_data"]
df_bulk_dft = df_bulk_dft[df_bulk_dft.source == "raul"]

df_dij = DF_dict["df_dij"]

# + jupyter={"outputs_hidden": true}
# #########################################################
import pickle; import os
path_i = os.path.join(
    # os.environ[""],
    "out_data",
    "ids_to_discard__proto_dupl.pickle")
with open(path_i, "rb") as fle:
    ids_to_drop = pickle.load(fle)
# #########################################################

# ids_to_drop

# + [markdown] Collapsed="false"
# ## Read prototype duplicate data

# + Collapsed="false"
# #############################################################################
import pickle; import os
path_i = os.path.join(
    "out_data",
    "duplicates_proto.pickle")
with open(path_i, "rb") as fle:
    duplicates_proto = pickle.load(fle)
# #############################################################################

# +
# df_dij.loc[["6qca6qnfme", "826imfvjm5"]]
# df_dij.loc["6qca6qnfme"]["919hvh8qx3"]
# df_dij.loc["6qca6qnfme"]["826imfvjm5"]

# df_dij.loc["826imfvjm5"]

df_dij.shape

# ["826imfvjm5"]
# ["6qca6qnfme", "826imfvjm5"]

# +
# df_bulk_dft.loc["826imfvjm5"]
# -

duplicates_proto

assert False

# + [markdown] Collapsed="false"
# # TEMP

# + Collapsed="false"
len(duplicates_proto)

# + Collapsed="false"
duplicates_proto

ids_to_drop = []
for duplicate_list in duplicates_proto:
    PRINT = False

    if "8snlxpnhmq" in duplicate_list:
        PRINT = True

    num_dupl = len(duplicate_list)

    dupl_ids_in_df = df_bulk_dft.index.intersection(duplicate_list)

    if PRINT:
        print(dupl_ids_in_df)
    
    if len(dupl_ids_in_df) == 0:
        ids_to_drop.extend(dupl_ids_not_in_df[1:])

    else:
        dupl_ids_not_in_df = [i for i in duplicate_list if i not in df_bulk_dft.index.tolist()]
        ids_to_drop.extend(dupl_ids_not_in_df)

    # print(dupl_ids_not_in_df)

    df_i = df_bulk_dft.loc[dupl_ids_in_df]

    df_not_most_stable = df_i[~(df_i.dH == df_i.dH.min())]
    dupl_ids_not_most_stable = df_not_most_stable.index.tolist()
    ids_to_drop.extend(dupl_ids_not_most_stable)

    # if df_i.shape[0] > 1:
    #      display(df_i)

# + Collapsed="false"
# Pickling data ######################################################
import os; import pickle
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)
with open(os.path.join(directory, "ids_to_discard__proto_dupl.pickle"), "wb") as fle:
    pickle.dump(ids_to_drop, fle)
# #####################################################################

# + Collapsed="false" active=""
#
#
#

# + Collapsed="false" jupyter={}
# for duplicate_list in duplicates_proto:
# #     print(duplicate_list)
    
#     tmp = [i for i in duplicate_list if i not in ids_to_drop]
# #     print(tmp)
    
#     if len(tmp) == 0:
#         print(duplicate_list)
