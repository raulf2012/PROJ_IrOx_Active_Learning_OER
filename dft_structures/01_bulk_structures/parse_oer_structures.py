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

# +
import os
import sys

from ase import io

# +
root_dir = os.path.join(
    os.environ["PROJ_irox"],
    "01_bulk_structures")

folders_list = [
    "01_IrO2/bulk_symm.cif",
    "02_IrO3/bulk.cif",
    "03_iro3_rutile-like/bulk.cif",
    "04_iro3_battery/bulk.cif",
    ]

main_systems = ["IrO2", "IrO3", "IrO3_battery", "IrO3_rutile-like"]


folders_dict = {
    "IrO2": "01_IrO2/bulk_symm.cif",
    "IrO3": "02_IrO3/bulk.cif",
    "IrO3_rutile-like": "03_iro3_rutile-like/bulk.cif",
    "IrO3_battery": "04_iro3_battery/bulk.cif",
    }


for folder_i in folders_list:

    bulk_path_i = os.path.join(root_dir, folder_i)
    atoms_i = io.read(bulk_path_i)
    
    
# -

for key, value in folders_dict.items():
    print(key)

atoms_i
