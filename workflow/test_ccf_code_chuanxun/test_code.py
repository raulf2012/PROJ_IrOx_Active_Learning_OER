# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.1.7
#   kernelspec:
#     display_name: Python [conda env:research-new]
#     language: python
#     name: conda-env-research-new-py
# ---

# + {"jupyter": {"source_hidden": true}}
import os
import sys

import pickle


# #############################################################################
from ccf import(
    struc2ccf,
    cal_ccf_d,
    cal_inter_atomic_d,
    d2ccf,
    weight_f,
    pearson_cc,
    gaussian_f,
    element_tag,
    cell_range,
    count_atoms_dict,
    )


sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import (
    bulk_dft_data_path,
    unique_ids_path,
    prototypes_data_path,
    static_irox_structures_path,
    oqmd_irox_data_path,
    voronoi_features_data_path,
    )
# -

# # Read Data

with open(static_irox_structures_path, "rb") as fle:
    df_structures = pickle.load(fle)

# +
atoms_i = df_structures.iloc[0]["atoms"]

atoms_i

# +

r_cut_off = 5.
# r_vector = [1., 5.,]
r_vector = 3.



struc2ccf(atoms_i, 1, [1, 2])


i_a_d =cal_inter_atomic_d(atoms_i, 5.)
d2ccf(i_a_d, r_cut_off, r_vector)
