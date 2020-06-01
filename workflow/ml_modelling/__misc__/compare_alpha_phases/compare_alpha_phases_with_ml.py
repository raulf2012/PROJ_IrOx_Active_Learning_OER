# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python [conda env:research-new]
#     language: python
#     name: conda-env-research-new-py
# ---

# # Import Modules

# +
import os
import sys

import numpy as np

from ase import io

from StructurePrototypeAnalysisPackage.ccf import struc2ccf
# -

from StructurePrototypeAnalysisPackage.ccf import (
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

# # Script Inputs

r_cut_off = 10.
r_vector = np.arange(1, 10, 0.02)

# # Read structure files

# +
struct_files = [
    "000__id-unique_8p8evt9pcg__id-short_207.cif",
    "22_179197-3835-5794-AlF3_F6Al2fixed.cif",
    "out.cif",
    ]

atoms_0 = io.read(struct_files[0])
atoms_1 = io.read(struct_files[1])
atoms_2 = io.read(struct_files[2])

# +
# atoms_i = row_i[atoms_key]

ccf_0 = struc2ccf(atoms_0, r_cut_off, r_vector)
ccf_1 = struc2ccf(atoms_1, r_cut_off, r_vector)
ccf_2 = struc2ccf(atoms_2, r_cut_off, r_vector)
# -

# d_ij = 
cal_ccf_d(ccf_0, ccf_1)
