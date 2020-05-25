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
# | - Import Modules
import os
print(os.getcwd())
import sys

import pickle

import numpy as np
import pandas as pd
#__|

# +
# /mnt/f/Dropbox/01_norskov/00_git_repos/PROJ_IrOx_Active_Learning_OER/workflow/ml_modelling

sys.path.insert(0, "/mnt/f/Dropbox/01_norskov/00_git_repos/PROJ_IrOx_Active_Learning_OER/workflow/ml_modelling")

from ml_methods import get_ml_dataframes, get_data_for_al

# +
DF_dict = get_ml_dataframes()

DF_dict.keys()

static_irox_structures = DF_dict["static_irox_structures"]

# +
static_irox_structures = static_irox_structures[static_irox_structures.source == "chris"]


static_irox_structures[static_irox_structures.stoich == "AB3"]
static_irox_structures[static_irox_structures.stoich == "AB2"]
