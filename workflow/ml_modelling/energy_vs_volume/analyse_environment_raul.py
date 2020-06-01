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

# + Collapsed="false" jupyter={}
import os
print(os.getcwd())

import sys

from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder

from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import SimplestChemenvStrategy, MultiWeightsChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments

import numpy as np

from pymatgen.core.structure import Structure

from pymatgen.io.ase import AseAtomsAdaptor

import ase
from ase.db import connect
from time import time
# -

# # Script Inputs

# +
t0 = time()

ids_to_run = [
    # 'vovgximhm2',
    # '8dce6kz2vf',
    # 'vhv39q6e9j',
    # '8ymh8qnl6o',
    # '6fcdbh9fz2',
    # '7qm56wxj8s',
    # 'mu6omk6k9l',
    # '6dzhcimdxs',
    
    # Had nan values
    # "vovl8d7wvi",
    ]

# + active=""
#
#
#

# + Collapsed="false"
db = connect('out_data/FinalStructures_1.db')
for row in db.select():
    # print('-------------------')
    # print(row.id)

    key_value_pairs = row.key_value_pairs
    structure_id = key_value_pairs["structure_id"]

    if structure_id in ids_to_run:
        print("RUNNING!")

        atoms = row.toatoms()
        Ir_indices = [i for i, s in enumerate(atoms.get_chemical_symbols())
                      if s == 'Ir']

        struct = AseAtomsAdaptor.get_structure(atoms)


        lgf = LocalGeometryFinder()
        lgf.setup_structure(structure=struct)

        se = lgf.compute_structure_environments(
            maximum_distance_factor=1.41,
            only_cations=False)


        strategy = MultiWeightsChemenvStrategy.stats_article_weights_parameters()

        lse = LightStructureEnvironments.from_structure_environments(
            strategy=strategy, structure_environments=se)

        isite = 0

        cor_env = []
        for isite in Ir_indices:
            c_env = lse.coordination_environments[isite]
            if len(c_env) == 0:
                continue
            cor_env += [c_env[0]['ce_symbol']]

        if len(cor_env) == 0:
            continue
        uniques, counts = np.unique(cor_env, return_counts=True)

        if len(uniques) > 1:
            coor = 'mixed'

            coor_l = [int(c.split(':')[-1]) for c in cor_env]

            mean_coor = np.mean(coor_l)
            #idx = np.argmax(counts)
            #coor = uniques[idx]
        else:
            coor = uniques[0]
            mean_coor = int(coor.split(':')[-1])

        print(coor, mean_coor)
        # db.update(id=row.id, coor_env=coor, mean_coor=mean_coor)
# -

print("Running time (s)", time() - t0)

# + active=""
#
#
#
#

# + Collapsed="false" jupyter={}
# struct

# se = lgf.compute_structure_environments(
#     maximum_distance_factor=1.41,
#     only_cations=False)
