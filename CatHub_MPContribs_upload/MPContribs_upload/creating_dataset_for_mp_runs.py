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

from pymatgen.util.provenance import StructureNL
from pymatgen.io.ase import AseAtomsAdaptor

from datetime import datetime

import json

# + active=""
#
#
# -

date = datetime.today().strftime('%Y-%m-%d')
date

# # Read DFT Data

# #########################################################
import pickle; import os
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/processing_bulk_dft/creating_final_dataset_for_upload",
    "out_data/df_dft_final_no_dupl.pickle")
with open(path_i, "rb") as fle:
    df_dft_final_no_dupl = pickle.load(fle)
# #########################################################

# +
df_ab2 = df_dft_final_no_dupl[df_dft_final_no_dupl.stoich == "AB2"]
df_ab3 = df_dft_final_no_dupl[df_dft_final_no_dupl.stoich == "AB3"]

df_ab2 = df_ab2.sort_values("dH")
df_ab3 = df_ab3.sort_values("dH")

df_ab2["energy_order"] = [i for i in range(df_ab2.shape[0])]
df_ab3["energy_order"] = [i for i in range(df_ab3.shape[0])]

# Reconstructing dataframe
df_dft_final_no_dupl = pd.concat([
    df_ab2, df_ab3
    ], axis=0)

# +
authors = [
    {"name": "Raul A. Flores", "email": "raulf2012@gmail.com"},
    # {"name": "Chris Paolucci", "email": "cp9wx@virginia.edu"},
    # {"name": "Kirsten Winther", "email": "winther@stanford.edu"},
    # {"name": "Ankit Jain", "email": "ankitjain.me.iitk@gmail.com"},
    # {"name": "Jose Antonio Garrido Torres", "email": "jagt@stanford.edu"},
    # {"name": "Muratahan Aykol", "email": "muratahan.aykol@tri.global"},
    # {"name": "Joseph Montoya", "email": "joseph.montoya@tri.global"},
    # {"name": "Jens Kehlet Nørskov", "email": "jkno@dtu.dk"},
    {"name": "Michal Bajdich", "email": "bajdich@slac.stanford.edu"},
    # {"name": "Thomas Bligaard", "email": "tbli@dtu.dk"},
    ]

remarks = [
    "Structure is part of the `active_learned_irox_polymorphs` dataset on MPContribs",
    "https://portal.mpcontribs.org/active_learned_irox_polymorphs/",
    ]
# -

for id_i, row_i in df_dft_final_no_dupl.iterrows():
    # print(id_i)

    unique_id = row_i.name
    stoich = row_i.stoich
    energy_order = row_i.energy_order

    # #########################################################
    energy_order_prefix = str(energy_order).zfill(4)

    filename = stoich + "_" + energy_order_prefix + "_" + unique_id
    filename += ".json"

    
    atoms = row_i.atoms
    struct = AseAtomsAdaptor().get_structure(atoms)

    extra_data = {
        "_MPContribs_Internal_ID": unique_id,
        }

    struct_NL = StructureNL(
        struct,
        authors,
        projects=None,
        references="",
        remarks=remarks,
        data=extra_data,
        # history=extra_data,
        created_at=date,
        )

    path_i = os.path.join("out_data", filename)
    with open(path_i,"w") as file:
        json.dump(
            struct_NL.as_dict(),
            file,
            indent=2,
            )
