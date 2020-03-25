from spglib import get_symmetry_dataset
from ase import Atoms
from ase.io import read
from sys import argv
import numpy as np

from ase.visualize import view

name = argv[1]
atoms = read(name)

images = [atoms.copy()]

spglibdata = get_symmetry_dataset((atoms.get_cell(),
                                   atoms.get_scaled_positions(),
                                   atoms.get_atomic_numbers()),
                                  symprec=1e-3)

spacegroup = spglibdata['number']

wyckoffs = spglibdata['wyckoffs']
print(spglibdata['number'])
print(wyckoffs)


s_name = spglibdata['international']
#str(spacegroup) + '_' + '_'.join(sorted(wyckoffs))

std_cell = spglibdata['std_lattice']
positions = spglibdata['std_positions']
numbers = spglibdata['std_types']

atoms = Atoms(numbers=numbers,
              cell=std_cell,
              pbc=True)

atoms.set_scaled_positions(positions)
atoms.wrap()
images += [atoms]

view(images)
new_name = name.rstrip('.cif') + '_conventional_' + str(spacegroup) + '.cif'
#print(new_name)
atoms.write(new_name)



