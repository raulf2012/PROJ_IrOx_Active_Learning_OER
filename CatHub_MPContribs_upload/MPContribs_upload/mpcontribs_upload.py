# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python [conda env:mpcontribs]
#     language: python
#     name: conda-env-mpcontribs-py
# ---

# # Import Modules

# + jupyter={}
import os
print(os.getcwd())
import sys

import json

import numpy as np
import pandas as pd

from pymatgen.io.ase import AseAtomsAdaptor

from mpcontribs.client import load_client
# -

# # Read IrOx DFT Data

# +
# #########################################################
import pickle; import os
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/processing_bulk_dft/creating_final_dataset_for_upload",
    "out_data/df_dft_final_no_dupl.pickle")
with open(path_i, "rb") as fle:
    df_dft = pickle.load(fle)
# #########################################################

df_dft = df_dft.drop(columns=["id", "form_e_chris", "path", "source"])

# +
df_dft = df_dft.sort_values(["stoich", "dH"])


# #########################################################
# df_dft = df_dft.iloc[0:16]
# df_dft = df_dft.sample(n=10)

# +
# %%capture

sys.path.insert(0, 
    os.path.join(
        os.environ["PROJ_irox"],
        "workflow/ml_modelling"))

from ml_methods import get_ml_dataframes
DF_dict = get_ml_dataframes(
    names=[
        'static_irox_structures_path',
        'df_prototype_dft_path',
        'df_prototype_static_path',
        ]
    )

df_prototype_static = DF_dict["df_prototype_static"]
df_prototype_dft = DF_dict["df_prototype_dft"]

static_irox_structures = DF_dict['static_irox_structures']


# -

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


# # Add 'Formula' Column to df

# +
def method(row_i):
    stoich_to_form_dict = {
        "AB2": "IrO2",
        "AB3": "IrO3"}

    stoich = row_i.stoich
    formula = stoich_to_form_dict.get(stoich)
    return(formula)

df_i = df_dft
df_i["formula"] = df_i.apply(
    method,
    axis=1)
# -

# # MPContribs

# +
import yaml
path_i = os.path.join(os.environ["PROJ_irox"], "config", "config.yml")
with open(path_i) as file:
    config_dict = yaml.load(file, Loader=yaml.FullLoader)

api_key = config_dict["mpcontrib"]["api_key"]

# +
project = 'active_learned_irox_polymorphs'

client = load_client(api_key)

# print(dir(client))
# print(dir(client.projects))
# -

# # Deleting all data to start over

# +
total_count = 1234
while total_count != 0:
    deleted = client.contributions.delete_entries(project=project).result()

    total_count = deleted["total_count"]
    num_deleted = deleted["count"]

    print(
        num_deleted, 'contribution(s) deleted',
        "|",
        total_count, "contribution(s) remaining")

# Delete entire project to start from scratch
# DON'T DO THIS (PATRICK HAS TO MANUALLY APPROVE PROJECT EVERYTIME)
if False:
    results = client.projects.delete_entry(pk=project).result()

# # client.projects.delete_entry?
# -

# # Create Project (Once)

# +
is_public = True
info = {"project": project}

# CREATE PROJECT FROM SCRATCH (DO THIS SELDOMLY)
# DON'T DO THIS (PATRICK HAS TO MANUALLY APPROVE PROJECT EVERYTIME)
if False:
# if True:
    client.projects.create_entry(project=info).result()

# +
if False:
    all_data = client.projects.get_entry(pk=project, _fields=['_all']).result()

# all_data
# -

# # Adding new data rows

contributions = dict()
for ind_i, row_i in df_dft.iterrows():
    stoich = row_i.stoich
    dH = row_i.dH
    formula = row_i.formula
    energy_pa = row_i.energy_pa
    num_atoms = row_i.num_atoms
    volume = row_i.volume
    volume_pa = row_i.volume_pa

    # #####################################################
    if ind_i in df_prototype_dft.index:
        row_proto_dft = df_prototype_dft.loc[ind_i]

        prototype_name_dft = row_proto_dft.p_name
        spacegroup_dft = int(row_proto_dft.spacegroup)
    else:
        print("Woops", ind_i)
        
        prototype_name_dft = ""
        spacegroup_dft = None

    # #########################################################
    if ind_i in df_prototype_static.index:
        row_proto_static_i = df_prototype_static.loc[ind_i]

        prototype_name_static = row_proto_static_i.p_name
        spacegroup_static = int(row_proto_static_i.spacegroup)
    else:
        print("Woops", ind_i)

        prototype_name_static = ""
        spacegroup_static = None

    # prototype_name_static = prototype_name_static

    # #####################################################
    # #####################################################
    dH = str(dH) + " eV/atom"
    dft_energy_per_atom = str(energy_pa) + " eV/atom"
    number_of_atoms = num_atoms
    volume = str(volume) + " angstrom**3"
    volume_pa = str(volume_pa) + " angstrom**3/atom"

    # #####################################################
    contributions[ind_i] = dict(
        contrib=dict(
            identifier=ind_i, project=project, is_public=is_public,
            # identifier="NA", project=project, is_public=is_public,
            data={
                "InternalID": row_i.name,

                "ΔH|formation": dH,
                "Formula": formula,
                "EnergyDFT": dft_energy_per_atom,
                "NumberOfAtoms": number_of_atoms,
                "Volume|UnitCell": volume,
                "Volume": volume_pa,

                "StructurePrototype|PreDFT": prototype_name_static,
                "StructurePrototype|PostDFT": prototype_name_dft,
                "SpaceGroupNumber|PreDFT": spacegroup_static,
                "SpaceGroupNumber|PostDFT": spacegroup_dft,
                },
            )
        )

# +
if True:
    contribs = []
    for key, val in contributions.items():
        contribs.append(val["contrib"])

    chunk_size = 20
    df_mp_list = []
    for contribs_chunk_i in chunks(contribs, chunk_size):

        created = client.contributions.create_entries(
            contributions=contribs_chunk_i).result()

        df_mp_i = pd.DataFrame(created["data"]).set_index("identifier")
        df_mp_list.append(df_mp_i)


df_mp = pd.concat(df_mp_list)
# -

# # Add Structures

for id_chunk_i in chunks(df_mp.index.tolist(), chunk_size):
    df_mp_i = df_mp.loc[id_chunk_i]

    structure_contribs = []
    for id_i, row_i in df_mp_i.iterrows():
        print(id_i)

        cid = row_i.id

        # #####################################################
        #  DFT Data ###########################################
        row_dft_i= df_dft.loc[id_i]

        formula = row_dft_i.formula
        atoms_final = row_dft_i.atoms

        # #####################################################
        # Static IrOx #########################################
        row_static_i = static_irox_structures.loc[id_i]

        atoms_init = row_static_i.atoms


        # #####################################################
        # #####################################################
        structure_final = AseAtomsAdaptor.get_structure(atoms_final)
        structure_init = AseAtomsAdaptor.get_structure(atoms_init)

        # #####################################################
        sdct = dict(contribution=cid,
            name=id_i + "_final",
            label="Final_DFT_Optimized",
            )
        sdct.update(structure_final.as_dict())
        structure_contribs.append(sdct)
        # print(id_i + "_final")

        sdct = dict(contribution=cid,
            name=id_i + "_init",
            # label="Initialstructuralprototype",
            label="Initial_Prototype",
            )
        sdct.update(structure_init.as_dict())
        structure_contribs.append(sdct)
        # print(id_i + "_init")

    sid = client.structures.create_entries(structures=structure_contribs).result()

assert False

# # MISC
# ---

# # Extracting data with `get_entries` methods

# +
# identifier = "mp-1234"

# client.contributions.get_entries(
#     project=project,
#     identifier=identifier,
#     # _fields=["formula"],
#     ).result()
#     # ).result()['data']

# # # client.projects.get_entries?

# + active=""
#
#
#
#
#

# + jupyter={}
# if False:
# # if True:
#     results = client.projects.update_entry(pk=project,
#         project={
#             # "is_public": False,
#             # "project": project,

#             "is_public": is_public,
#             "title": "Active Learned IrOx Polymorphs",
#             "owner": "raulf2012@gmail.com",
#             # "authors": "R. Flores, W. Kirsten",
#             "authors": "Raul A. Flores, Christopher Paolucci, Kirsten T. Winther, Ankit Jain, Jose Antonio Garrido Torres, Muratahan Aykol, Joseph Montoya, Jens K. Nørskov",

#             "description": " ".join("""
#                 Materials science is primarily concerned with the underlying relationship between a material's structure and functionality,
#                 where the knowledge of viable polymorphic forms of crystals plays an indispensable role.
#                 Machine-learning based surrogate models have the potential to accelerate this process of creating the knowledge-base for materials polymorphs for target applications in under-explored chemistries.
#                 Herein, we report on a readily generalizable active-learning (AL) accelerated algorithm for the targeted identification of novel and stable IrOx (x=2 or 3) polymorphs and subsequent thermochemical analyses of the activity of these discovered structures towards the oxygen evolution reaction (OER).
#                 We demonstrate that compared to a random search,
#                 the AL framework more than doubles the efficiency of using DFT to find stable polymorphs out of a large array of prototypical structures.
#                 We find nearly 195 IrO2 polymorphs within the thermodynamic synthesizability limit and reaffirm the rutile ground state.
#                 For IrO3, we find 74 unique synthesizable polymorphs and report a previously unknown FeF3-like ground state.
#                 The algorithm is exceptionally adept at quickly picking out the most stable polymorphs, with the most stable α-IrO3 phase discovered on average in only 4.3 generations.
#                 An analysis of the structural properties of these metastable polymorphs reveals that octahedral local coordination environments are preferred for all low energy structures.
#                 """.replace("\n", "").split()),
#             # "urls": None,
#             # "urls": dict(),
#             "urls": {
#                 "PaperGit":    "https://github.com/raulf2012/PAPER_IrOx_Active_Learning_OER",
#                 "ProjGit":     "https://github.com/raulf2012/PROJ_IrOx_Active_Learning_OER",
#                 "PaperURL":    "https://github.com/raulf2012/PROJ_IrOx_Active_Learning_OER",
#                 # "": "",
#                 },

#             "other": {
#                 "InternalID": "Unique ID used internally, including for posterity in case anybody wants to dig through the project's Git repo",
#                 "ΔH|formation": "Heat of formation",
#                 # "Formula": "",
#                 "EnergyDFT": "Raw DFT VASP energy",
#                 # "NumberOfAtoms": "",
#                 # "Volume|UnitCell": "",
#                 "Volume": "Total computational cell volume",
#                 "StructurePrototype|PreDFT": "Structural prototype of initial pre-DFT optimized structure candidate",
#                 "StructurePrototype|PostDFT": "Structural prototype of post-DFT relaxed structure",
#                 "SpaceGroupNumber|PreDFT": "Space group of pre-DFT structure candidate",
#                 "SpaceGroupNumber|PostDFT": "Space group of post-DFT structure candidate",

#                 },

#             }
#         )
#     results.result()

# + jupyter={}
# # #########################################################
# import pickle; import os
# path_i = os.path.join(
#     os.environ["PROJ_irox"],
#     "workflow/ml_modelling/processing_bulk_dft/creating_final_dataset_for_upload",
#     "out_data/df_dft_final_no_dupl.pickle")
# with open(path_i, "rb") as fle:
#     df_dft_new = pickle.load(fle)
# # #########################################################

# + jupyter={}
# df_dft_new[df_dft_new.stoich == "AB2"].sort_values("dH").iloc[4:10]

# # df_prototypes_dft.loc["xw9y6rbkxr"]
# # df_prototype_dft.loc["b5cgvsb16w"]
# # df_prototype_dft.loc["8p8evt9pcg"]
# # df_prototype_dft.loc["949rnem5z2"]
# # df_prototype_dft.loc["6pvt6r95ve"]
# # df_prototype_dft.loc["8l919k6s7p"]

# # df_prototype_dft.loc["8l919k6s7p"]
# # df_prototype_dft.loc["n36axdbw65"]
# # df_prototype_dft.loc["clc2b1mavs"]
# # df_prototype_dft.loc["ck638t75z3"]
# # df_prototype_dft.loc["mkbj6e6e9p"]
# # df_prototype_dft.loc["b49kx4c19q"]
# df_prototype_dft.loc["85z4msnl6o"]

# + jupyter={}
# assert False

# + jupyter={}
# print("Why are these numbers different")
# print(df_prototype_dft.shape)
# print(df_dft.shape)
