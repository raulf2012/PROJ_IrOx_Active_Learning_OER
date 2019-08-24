#!/usr/bin/env python

"""Run VASP job on slab.

Author(s): Michal Badich, Chris Paolucci, Raul A. Flores
"""

#| - Import Modules
import os
import ase.calculators.vasp as vasp_calculator
import subprocess
from ase import io
import numpy as np

from ase_modules.ase_methods import clean_up_dft
#__|

#| - Read Atoms Object
if os.path.isfile("init.cif"):
    atoms = io.read('init.cif')
elif os.path.isfile("init.traj"):
    atoms = io.read('init.traj')
#__|

subprocess.call('cp -rf OUTCAR OUTCAR_$(date +%s)', shell=True)
subprocess.call('cp -rf moments.traj moments.traj_$(date +%s)', shell=True)

#| - k-points
unitcell = atoms.get_cell()
print(unitcell)

kpoints = []
for i in [0, 1, 2]:
    l_var = np.sqrt(unitcell[i][0]**2 + unitcell[i][1]**2 + unitcell[i][2]**2)
    k = int(4 * 5 / l_var)
    if(k > 0):
        kpoints.append(k)
    else:
        kpoints.append(1)

kpar = 10
if kpoints[0] * kpoints[1] * kpoints[2] < 10:
    if(kpoints[0] * kpoints[1] * kpoints[2] % 2 == 0):
        kpar = 2
    else:
        kpar = 1

print(kpoints)
#__|

#| - Calculator

calc_params = dict(
    potim=0.03,
    encut=600,
    xc='PBE',
    #setups={'O': '_s', 'C': '_s'},
    gga='PE',
    #kpts  = (2,2,1),
    kpts=(kpoints[0], kpoints[1], kpoints[2]),
    kpar=5,
    npar=6,
    gamma=True,  # Gamma-centered (defaults to Monkhorst-Pack)
    ismear=0,
    inimix=0,
    amix=0.2,
    bmix=0.0001,
    amix_mag=0.1,
    bmix_mag=0.0001,
    #nupdown=0,
    nelm=220,
    sigma=0.05,
    #algo = 'normal',
    algo='fast',
    ibrion=2,
    isif=2,
    #isif=2,
    ediffg=1e-3,  # forces
    ediff=1e-6,  # energy conv.
    #nedos=2001,
    prec='High',
    nsw=150,  # Don't use the VASP internal relaxation, only use ASE
    lvtot=False,
    ispin=2,
    ldau=False,
    ldautype=2,
    #lreal  = 'False',#more accuratre volume
    lreal='auto',
    lasph=True,
    ldau_luj={
        'Ni': {'L': 2, 'U': 6.45, 'J': 0.0},
        'Co': {'L': 2, 'U': 3.32, 'J': 0.0},
        'Cr': {'L': 2, 'U': 3.5, 'J': 0.0},
        'Fe': {'L': 2, 'U': 5.3, 'J': 0.0},
        'Ce': {'L': 3, 'U': 4.50, 'J': 0.0},
        'V': {'L': 2, 'U': 3.25, 'J': 0.0},
        'Mn': {'L': 2, 'U': 3.75, 'J': 0.0},
        'Ti': {'L': 2, 'U': 3.00, 'J': 0.0},
        'W': {'L': -1, 'U': 0.0, 'J': 0.0},
        'O': {'L': -1, 'U': 0.0, 'J': 0.0},
        'C': {'L': -1, 'U': 0.0, 'J': 0.0},
        'Au': {'L': -1, 'U': 0.0, 'J': 0.0},
        'Cu': {'L': -1, 'U': 0.0, 'J': 0.0},
        'H': {'L': -1, 'U': 0.0, 'J': 0.0},
        },
    ldauprint=2,
    #lmaxmix=4,
    lmaxmix=6,
    lorbit=11,
    #idipol=3,
    #dipol=(0, 0, 0.5),
    #ldipol=True
    )

# Reading VASP parameters from file and merging with params in script
from ase_modules.dft_params import VASP_Params
VP = VASP_Params(load_defaults=False)
VP.load_params()
calc_params.update(VP.params)

calc = vasp_calculator.Vasp(**calc_params)
#__|

atoms.set_calculator(calc)
atoms.get_potential_energy()

io.write('out.traj', atoms)
io.write('out.traj', atoms)

clean_up_dft()
