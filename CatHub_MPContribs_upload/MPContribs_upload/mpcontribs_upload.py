# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
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

# +
import ase

ase.__version__
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
df_dft = df_dft.sort_values(["stoich", "dH"])


# #########################################################
# df_dft = df_dft.iloc[0:16]
# df_dft = df_dft.sample(n=10)

# +
# df_dft = df_dft.loc[[
#     '94x5nkmjmf',
#     'cazivwbq94',
#     '8kvd8qnim4',
#     '656qniby7j',
#     '9yn4m16ux1',
#     'cfzam1mdbf',
#     'zy9dzknhnj',
#     'zanqv2xtvk',
#     'njntmu9w93',
#     'mrbine8k72',

#     'cg8p7fxq65',
#     '64cg6j9any',
#     '85z4msnl6o',
#     'xozr8f7p7g',
#     '949rnem5z2',
#     'mkmsvkcyc5',
#     'vwxfn3blxi',
#     'nrml6dms9l',
#     ]]

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

# #########################################################
import pickle; import os
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "CatHub_MPContribs_upload/MPContribs_upload/duplicate_MP_entries",
    "out_data/df_mp_dupl.pickle")
with open(path_i, "rb") as fle:
    df_mp_dupl = pickle.load(fle)
# #########################################################

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
# -

if False:
    all_data = client.projects.get_entry(pk=project, _fields=['_all']).result()
    # all_data

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
    # Setting the MP id as identifier if available ########
    identifier_i = ind_i
    if ind_i in df_mp_dupl.index:
        row_i = df_mp_dupl.loc[ind_i]
        # row_i = df_mp_dupl.loc["cg8p7fxq65"]

        mp_duplicates = row_i.mp_duplicates
        if len(mp_duplicates) > 1:
            print("There is more than 1 MP duplicate found")
        if len(mp_duplicates) == 0:
            print("There is no MP duplicate for this entry, but there should be")

        mp_duplicate = mp_duplicates[0]
        identifier_i = mp_duplicate    
    
    # #####################################################
    contributions[ind_i] = dict(
        contrib=dict(
            identifier=identifier_i,
            # identifier=ind_i,
            # identifier="NA", project=project, is_public=is_public,
            project=project, is_public=is_public,
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
# contributions

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

# + jupyter={}
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
