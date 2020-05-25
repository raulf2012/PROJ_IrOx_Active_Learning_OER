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

# # Import Modules

# +
import os
print(os.getcwd())

import sys

import pickle

import pandas as pd

# import bulk_enumerator as be
import time

from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.ase import AseAtomsAdaptor


# #############################################################################
from ase import io

from catkit.build import bulk
from catkit.gen.symmetry import get_standardized_cell
from catkit.gen.symmetry import Symmetry

# #############################################################################
sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))     

from proj_data_irox import (
    bulk_dft_data_path,
    unique_ids_path,
    prototypes_data_path,
    static_irox_structures_path,
    oqmd_irox_data_path,
    )

# +
import ase

ase.__version__

ase
# -

# # Read Data

# +
with open(bulk_dft_data_path, "rb") as fle:
    df_bulk_dft = pickle.load(fle)

df_bulk_dft = df_bulk_dft[df_bulk_dft.source == "raul"]
df_bulk_dft = df_bulk_dft.drop(
    columns=["form_e_chris", "id", "id_old", "path", "energy"],
    )

# df_bulk_dft = df_bulk_dft.iloc[0:5]

# +
ids = [
    "8p8evt9pcg",
    "macixavwv3",
    "zimixdvdxd",
    "8ivkxwnhva",
    "9lmkmh8s8r",
    "vp7gxqv191",
    "9txdvicqcf",
    "8k7expx2bp",
    "xwvhnh9sx4",
    "vlbdnoxlnh",
    "xjxdzi73bg",
    "xg6exl6rmp",
    "6fcdbh9fz2",
    "9rz5nl9g6o",
    "xw9y6rbkxr",
    "v5cym3nycg",
    ]

df_bulk_dft = df_bulk_dft.loc[ids]


# + active=""
#
#
#
# -

def perform_symm_op(atoms, tol=1e-3, ang_tol=-1):
    """
    """
    number_of_atoms_0 = atoms.get_number_of_atoms()
    # #########################################################################
    Symm_0 = Symmetry(atoms, tol=tol, ang_tol=ang_tol)

    # #####################################
    data = Symm_0.data
    number = data["number"]
    hall_number = data["hall_number"]
    international = data["international"]
    hall = data["hall"]
    pointgroup = data["pointgroup"]

    # #####################################
    lattice_name = Symm_0.get_lattice_name()

    


    # #########################################################################
    # #########################################################################
    # Analyzing standardized atoms symmetry ###################################
    atoms_prim = get_standardized_cell(atoms, primitive=True, tol=tol)
    number_of_atoms_1 = atoms_prim.get_number_of_atoms()

    # Symm_1 = Symmetry(atoms_prim, tol=tol, ang_tol=ang_tol)
    # data = Symm_1.data
    # number_tmp = data["number"]
    # hall_number_tmp = data["hall_number"]
    # international_tmp = data["international"]
    # hall_tmp = data["hall"]
    # pointgroup_tmp = data["pointgroup"]
    # lattice_name_tmp = Symm_1.get_lattice_name()


    out_dict = dict(
        number=number,
        hall_number=hall_number,
        international=international,
        hall=hall,
        pointgroup=pointgroup,
        lattice_name=lattice_name,
        atoms_prim=atoms_prim,
        number_of_atoms_0=number_of_atoms_0,
        number_of_atoms_1=number_of_atoms_1,

        # number_2=number_tmp,
        # hall_number_2=hall_number_tmp,
        # international_2=international_tmp,
        # hall_2=hall_tmp,
        # pointgroup_2=pointgroup_tmp,
        # lattice_name_2=lattice_name_tmp,
        # # atoms_prim_2=atoms_prim,

        )

    return(out_dict)


# +
def method(row_i):
    atoms = row_i.atoms

    # tol = 9e-3
    tol = 9e-2
    # tol = 5e-4

    # ang_tol = -1
    ang_tol = -1

    symm_info = perform_symm_op(
        atoms,
        tol=tol,
        ang_tol=ang_tol)

    row_i = pd.Series()
    for key, value in symm_info.items():
        row_i[key] = value
    return(row_i)


df_symm = df_bulk_dft.apply(method, axis=1)

# + active=""
#
#
#
#
#

# +
ids_sub = [
    "8p8evt9pcg",
    "macixavwv3",
    "9lmkmh8s8r",
    "zimixdvdxd",
    ]

df_symm.loc[ids_sub]

# +
id_i = "9lmkmh8s8r"
atoms = df_symm.loc[id_i].atoms_prim

atoms.get_number_of_atoms()

atoms.write("out_data/" + id_i + "_reduced_2.cif")

# +
from ase_modules.ase_methods import view_in_vesta

df_i = df_symm.loc[ids_sub]
atoms_list = df_i.atoms_prim.tolist()
names = df_i.index.tolist()

# view_in_vesta(atoms_list, ase_gui=False, name_list=names)

# + active=""
#
#
#
#
#

# +
# row_i = df_bulk_dft.iloc[0]


# atoms = row_i.atoms
# tolerance = 1e-5

# Symm_0 = Symmetry(atoms, tol=tolerance, ang_tol=-1)
# atoms_stan = get_standardized_cell(atoms, primitive=True, tol=tolerance)

# print(atoms_stan.symbols)

# +
# df_bulk_dft.loc[ids].to_csv("tmp_df_bulk_dft.csv")

# df_symm.loc[ids].to_csv("tmp.csv")

# +
# import numpy as np
# tmp = np.array(df_symm.number - df_symm.number_2).tolist()

# [i for i in tmp if i != 0]
