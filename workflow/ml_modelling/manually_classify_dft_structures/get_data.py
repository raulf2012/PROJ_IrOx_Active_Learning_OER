#!/usr/bin/env python

"""TEMP.

Author(s): Raul A. Flores
"""

# | - IMPORT MODULES
import os
import sys
import pickle

import time
t0 = time.time()
sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import bulk_dft_data_path

# #############################################################################
# from ase_modules.ase_methods import view_in_vesta

# import pandas as pd

# from pymatgen.io.ase import AseAtomsAdaptor
# from pymatgen.analysis.local_env import NearNeighbors, VoronoiNN, site_is_of_motif_type
#__|


# | - Read Data
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling",
    "ccf_similarity_analysis/out_data",
    "all_ids_to_elim_1.pickle")
with open(path_i, "rb") as fle:
    all_ids_to_elim = pickle.load(fle)

with open(bulk_dft_data_path, "rb") as fle:
    df_bulk_dft = pickle.load(fle)
#__|


print("df_bulk_dft.shape:", df_bulk_dft.shape)

df_bulk_dft = df_bulk_dft[
    (df_bulk_dft["source"] != "oqmd") & \
    # (df_bulk_dft["source"] != "raul_oer") & \
    (df_bulk_dft["source"] != "chris") & \
    [True for i in range(len(df_bulk_dft))]
    ]

df_bulk_dft = df_bulk_dft.sort_values("energy_pa")

# print("df_bulk_dft.shape:", df_bulk_dft.shape)
# df_bulk_dft = df_bulk_dft.drop(all_ids_to_elim)
# print("df_bulk_dft.shape:", df_bulk_dft.shape)
