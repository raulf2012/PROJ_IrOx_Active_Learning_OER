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
print(os.getcwd())

from ase.visualize import view
from ase import io
# -

atoms = io.read("iro3_100_o_covered.cif")

# +
new_cell = [
    2 * atoms.cell[0],
    2 * atoms.cell[1],
    atoms.cell[2],
    ]

atoms.set_cell(new_cell)
# -

atoms.center()

# +
# atoms.write("00_out.cif")

# +
# view(atoms, viewer="x3d")

# view(atoms, viewer="ngl")

# +
# view(atoms)

# +
# from ipywidgets import interact, interactive, fixed, interact_manual
# import ipywidgets as widgets

# def f(x):
#     return 2 * x

# interact(f, x=10);

# +
# import nglview as np

# from nglview.datafiles import ASE_Traj

# ASE_Traj

# +
# import MDAnalysis as mda
import nglview as nv
from nglview.datafiles import PDB, XTC

# u = mda.Universe(PDB, XTC)

# protein = u.select_atoms('protein')

# +
# w = nv.ASEStructure(atoms)
# nv.ASEStructure?

# w
# -

nv.__version__

nv.show_ase(atoms)

nv.demo()
