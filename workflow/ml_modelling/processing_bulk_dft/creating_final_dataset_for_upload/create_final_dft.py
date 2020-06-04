# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_irox]
#     language: python
#     name: conda-env-PROJ_irox-py
# ---

# +
import os
print(os.getcwd())
import sys

import pandas as pd

# +
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling"))

from ml_methods import get_ml_dataframes

# +
DF_dict = get_ml_dataframes()

bulk_dft_data = DF_dict['bulk_dft_data']
unique_ids = DF_dict['unique_ids']
prototypes_data = DF_dict['prototypes_data']
static_irox_structures = DF_dict['static_irox_structures']
static_irox_structures_kirsten = DF_dict['static_irox_structures_kirsten']
oqmd_irox_data = DF_dict['oqmd_irox_data']
df_features_pre_opt = DF_dict['df_features_pre_opt']
df_features_pre_opt_kirsten = DF_dict['df_features_pre_opt_kirsten']
df_features_post_opt = DF_dict['df_features_post_opt']
oer_bulk_structures = DF_dict['oer_bulk_structures']
df_ccf = DF_dict['df_ccf']
df_dij = DF_dict['df_dij']
ids_to_discard__too_many_atoms = DF_dict['ids_to_discard__too_many_atoms']


duplicates = DF_dict.get("ids_duplicates")
# duplicates = ids_duplicates_dict

# +
bulk_dft_data = bulk_dft_data[bulk_dft_data.source == "raul"]
# bulk_dft_data.shape

static_irox_structures = static_irox_structures[static_irox_structures.source == "chris"]
# static_irox_structures.shape

# +
# %%capture

# static_irox_structures["num_of_atoms"] = [i.get_number_of_atoms() for i in static_irox_structures.atoms]
static_irox_structures.loc[:, "num_of_atoms"] = [i.get_number_of_atoms() for i in static_irox_structures.atoms]

static_irox_structures = static_irox_structures[static_irox_structures.num_of_atoms <= 75]
# static_irox_structures = static_irox_structures[static_irox_structures.num_of_atoms < 75] 
# -

print("TEMP:", static_irox_structures.shape[0])

static_irox_structures[static_irox_structures.stoich == "AB2"].shape

# +
static_ind = static_irox_structures.index

dft_ind = bulk_dft_data.index

intersection_indices = static_ind.intersection(dft_ind)

bulk_dft_data = bulk_dft_data.loc[intersection_indices]

bulk_dft_data.shape

# +
bulk_dft_data_ab2 = bulk_dft_data[bulk_dft_data.stoich == "AB2"]
bulk_dft_data_ab3 = bulk_dft_data[bulk_dft_data.stoich == "AB3"]

print("bulk_dft_data_ab2.shape:", bulk_dft_data_ab2.shape)
print("bulk_dft_data_ab3.shape:", bulk_dft_data_ab3.shape)
# -

# # Removing Post-DFT Duplicates

# +
bulk_dft_data_ab2 = bulk_dft_data_ab2.drop(
    bulk_dft_data_ab2.index.intersection(duplicates["AB2"])
    )
print(bulk_dft_data_ab2.shape)

bulk_dft_data_ab3 = bulk_dft_data_ab3.drop(
    bulk_dft_data_ab3.index.intersection(duplicates["AB3"])
    )
print(bulk_dft_data_ab3.shape)

# +
print(len(duplicates["AB2"]))

print(len(list(set(duplicates["AB2"]))))

# +
df_dft = pd.concat([
    bulk_dft_data_ab2,
    bulk_dft_data_ab3,
    ])

# df_dft
# -

df_dft[df_dft.stoich == "AB3"].sort_values("dH")

# Pickling data ###########################################
import os; import pickle
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)
with open(os.path.join(directory, "df_dft_final_no_dupl.pickle"), "wb") as fle:
    pickle.dump(df_dft, fle)
# #########################################################

print(20 * "# # ")
print("All done!")
assert False
