#!/usr/bin/env python

"""Cut slabs from bulk structures using ASE, pymatgen, and CatKit.

Author: Raul A. Flores
"""

#| - Import Modules
from ase import io
# from ase.visualize import view

from atoms_objects.slab_generation import (
    cut_slab_ase,
    cut_slab_pymatgen,
    cut_slab_catkit,
    )

from catkit.gen import utils
#__|

#| - Script Inputs
facet = (1, 0, 0)
#__|

#| - Read Bulk Structure
# bulk = io.read("init.cif")
bulk = io.read("bulk.cif")

# bulk = utils.get_spglib_cell(bulk)
# bulk.write("bulk_standard.cif")
#__|

#| - ASE
slab_ase = cut_slab_ase(
    bulk,
    facet,
    layers=6,
    vacuum=8,
    )
slab_ase.write("out_slab_ase.cif")
#__|

#| - pymatgen
slab_pymatgen = cut_slab_pymatgen(
    bulk,
    facet,
    min_slab_size=18.,
    min_vacuum_size=10.,
    )
slab_pymatgen.write("out_slab_pymatgen.cif")
#__|

#| - CatKit
slab_catkit = cut_slab_catkit(
    bulk,
    facet,
    slab_thickness=25,
    vacuum=8.,
    )
slab_catkit = slab_catkit[0]
slab_catkit.write("out_slab_catkit.cif")
#__|
