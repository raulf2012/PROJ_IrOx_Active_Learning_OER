# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:hydrogen
#     text_representation:
#       extension: .py
#       format_name: hydrogen
#       format_version: '1.3'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_IrOx_Active_Learning_OER]
#     language: python
#     name: conda-env-PROJ_IrOx_Active_Learning_OER-py
# ---

# %% [markdown]
# # Import Modules

# %% jupyter={}
import os
print(os.getcwd())
import sys

import pandas as pd

import pickle

# %%
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow"))
    # "workflow/ml_modelling/00_ml_workflow/191102_new_workflow"))
from al_data import al_data_files_dict, main_AB2_run, main_AB3_run

# %%
main_AB2_run
main_AB3_run

# %% [markdown]
# # Script Input

# %%
ab2_file_list = [
    main_AB2_run,
    ]

ab3_file_list = [
    main_AB3_run,
    ]

file_list_dict = dict(
    AB2=ab2_file_list,
    AB3=ab3_file_list,
    )


# %%
# #############################################################################
def get_duplicates_list(stoich_i, file_list_dict=None):
    """
    """
    file_list = file_list_dict[stoich_i]

    duplicates_lists = []
    for file in file_list:
        path_i = file

        print(path_i)

        with open(path_i, "rb") as fle:
            AL = pickle.load(fle)


        last_gen = list(AL.al_gen_dict.keys())[-1]
        AL_i = AL.al_gen_dict[last_gen]

        model = AL_i.model
        model.sort_values("y_real")

        duplicates_i = model[model.duplicate == True].index.tolist()
        # print(len(duplicates_i))

        duplicates_lists.append(duplicates_i)

    # #########################################################################
    # Checking that all duplicates lists are the same #########################
    duplicates_are_the_same_list = []
    for duplicates_i in duplicates_lists:
        for duplicates_j in duplicates_lists:
            duplicates_are_the_same = duplicates_j == duplicates_i
            duplicates_are_the_same_list.append(duplicates_are_the_same)
    duplicates_are_the_same_final = all(duplicates_are_the_same_list)
    assert duplicates_are_the_same_final, "IJDSFIISD"


    return(duplicates_i)

# %%
duplicates_ab2 = get_duplicates_list("AB2", file_list_dict=file_list_dict)
duplicates_ab3 = get_duplicates_list("AB3", file_list_dict=file_list_dict)

duplicates_dict = dict(
    AB2=duplicates_ab2,
    AB3=duplicates_ab3,
    )

# %%
"6fcdbh9fz2" in duplicates_dict["AB3"]

# %%
# # Pickling data ######################################################
# directory = "out_data"
# if not os.path.exists(directory): os.makedirs(directory)
# with open(os.path.join(directory, "duplicates.pickle"), "wb") as fle:
#     pickle.dump(duplicates_dict, fle)
# # #####################################################################

# %% [markdown]
# # Constructing Duplicates Manually (Without AL Run)

# %%
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling"))
from ml_methods import get_ml_dataframes
from ccf_similarity.ccf import CCF

# %%
DF_dict = get_ml_dataframes()

df_dft = DF_dict.get("bulk_dft_data")
df_dij = DF_dict.get("df_dij")
ids_to_discard__too_many_atoms = DF_dict.get("ids_to_discard__too_many_atoms")


df_dft = df_dft[df_dft.source == "raul"]
df_dft = df_dft.drop(columns=["path", "id_old", "id", "form_e_chris", "atoms"])

# %%
id_i = "cubqbpzd7k"
id_i in df_dft.index

# %%
# assert False

# %%
ids_to_discard__too_many_atoms

# %%
df_dij

# %%
print(df_dft.shape)
df_dft = df_dft.drop(
    index=df_dft.index.intersection(ids_to_discard__too_many_atoms)
    )
print(df_dft.shape)

df_dij = df_dij.drop(
    index=df_dij.index.intersection(ids_to_discard__too_many_atoms)
    )


# %%
def get_duplicates_list_manually(
    stoich_i=None, 
    df_dft=None,
    df_dij=None,
    ):
    # #########################################################
    # df_dft = df_dft[df_dft.source == "raul"]
    df_dft = df_dft[df_dft.stoich == stoich_i]

    df_dij = df_dij.loc[df_dft.index, df_dft.index]


    CCF_i = CCF(df_dij=df_dij, d_thresh=0.02)


    ids_to_drop = []
    for id_i in df_dft.index.tolist():
        simil_dict_i = CCF_i.i_all_similar(id_i)
        if simil_dict_i is not None:
            similar_ids = [id_i] + list(simil_dict_i.keys())
            df_i = df_dft.loc[similar_ids]
            ids_to_drop_i = df_i.sort_values("energy_pa").iloc[1:].index.tolist()
            ids_to_drop.extend(ids_to_drop_i)

    ids_to_drop__duplicates = ids_to_drop
    ids_to_drop__duplicates = list(set(ids_to_drop__duplicates))
    
    return(ids_to_drop__duplicates)

# %%
# assert False

# %%
duplicates_ab2_manual = get_duplicates_list_manually(stoich_i="AB2", df_dft=df_dft, df_dij=df_dij)

duplicates_ab3_manual = get_duplicates_list_manually(stoich_i="AB3", df_dft=df_dft, df_dij=df_dij)

duplicates_dict_manual = dict(
    AB2=duplicates_ab2_manual,
    AB3=duplicates_ab3_manual,
    )

# %%
# Pickling data ######################################################
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)
with open(os.path.join(directory, "duplicates.pickle"), "wb") as fle:
    pickle.dump(duplicates_dict_manual, fle)
# #####################################################################

# %% [raw]
#
#
#

# %% [markdown]
# # Comparing duplicate lists constructed from AL to those constructed manually

# %%
print(len(duplicates_ab3_manual))
print(len(set(duplicates_ab3)))

print("")
for i in duplicates_ab3_manual:
    if i not in duplicates_ab3:
        print(i)

print("")

for i in duplicates_ab3:
    if i not in duplicates_ab3_manual:
        print(i)

# %%
#           duplicates_ab2_manual
print(len(duplicates_ab2_manual))
print(len(set(duplicates_ab2)))

print("")
for i in duplicates_ab2_manual:
    if i not in duplicates_ab2:
        print(i)

print("")

for i in duplicates_ab2:
    if i not in duplicates_ab2_manual:
        print(i)

# %% [markdown]
# I think that it is best to go with the duplicates processed manually
#
# It looks like the AL runs are missing something

# %% [markdown]
# # TEST

# %%
stoich_i = "AB2"

# stoich_i=None,
df_dft=df_dft
df_dij=df_dij

# %%
# def get_duplicates_list_manually(
# stoich_i=None, 
# df_dft=None,
# df_dij=None,
# ):

# #########################################################
# df_dft = df_dft[df_dft.source == "raul"]
df_dft = df_dft[df_dft.stoich == stoich_i]

df_dij = df_dij.loc[df_dft.index, df_dft.index]


CCF_i = CCF(df_dij=df_dij, d_thresh=0.02)


ids_to_drop = []
for id_i in df_dft.index.tolist():
    simil_dict_i = CCF_i.i_all_similar(id_i)
    if simil_dict_i is not None:
        similar_ids = [id_i] + list(simil_dict_i.keys())
        df_i = df_dft.loc[similar_ids]
        ids_to_drop_i = df_i.sort_values("energy_pa").iloc[1:].index.tolist()
        ids_to_drop.extend(ids_to_drop_i)

ids_to_drop__duplicates = ids_to_drop
ids_to_drop__duplicates = list(set(ids_to_drop__duplicates))

# return(ids_to_drop__duplicates)

# %%
len(ids_to_drop__duplicates)

df_dft
df_dij.index

# ids_to_drop__duplicates

# df_dft.stoich.unique()

# %%
df_dft.sort_values("dH")

# %%
# df_i.sort_values("dH")

"64cg6j9any" in ids_to_drop__duplicates

# # CCF_i.i_all_similar("64cg6j9any")

# %%
assert False

# %% [markdown]
# # Comparing old and new duplicates

# %% jupyter={}
# #############################################################################
path_i = os.path.join(
    os.environ["PROJ_irox_2"],
    "FIGS_IrOx_Active_Learning_OER/01_figures/00_main_publ_figs/03_E_vs_V_coord/scripts",
    "old.duplicates.pickle")
with open(path_i, "rb") as fle:
    duplicates_dict_old = pickle.load(fle)
# #############################################################################

# #############################################################################
path_i = os.path.join(
    os.environ["PROJ_irox_2"],
    "FIGS_IrOx_Active_Learning_OER/01_figures/00_main_publ_figs/03_E_vs_V_coord/scripts",
    "new.duplicates.pickle")
with open(path_i, "rb") as fle:
    duplicates_dict_new = pickle.load(fle)
# #############################################################################

# %% jupyter={}
duplicates_dict_new["AB2"] == duplicates_dict_old["AB2"]

print(len(duplicates_dict_new["AB2"]))
print(len(duplicates_dict_old["AB2"]))

# %% [raw]
#
#
#
#

# %%
# for id in duplicates_dict_old["AB2"]:
#     if id not in duplicates_dict_new["AB2"]:
#         print(id)

# DF_dict = get_ml_dataframes()

# df_dft = DF_dict.get("bulk_dft_data")
# df_dij = DF_dict.get("df_dij")

# # #########################################################
# df_dft = df_dft[df_dft.source == "raul"]
# df_dft = df_dft[df_dft.stoich == "AB3"]

# df_dij = df_dij.loc[df_dft.index, df_dft.index]

# df_bulk_dft.index

# df_i

# df_bulk_dft.loc["pudomile_09"]

# sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
# from proj_data_irox import unique_ids_path

# stoich_i = "AB3"
# stoich_i = "AB2"

# drop_duplicates = False

# df_ids = pd.read_csv(unique_ids_path)

# al_data_files_dict["AB3"].keys()

# main_AB3_run

# "6fcdbh9fz2" in duplicates_ab3_manual

# "6fcdbh9fz2" in duplicates_dict_manual["AB3"]
