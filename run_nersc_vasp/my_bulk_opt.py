#!/usr/bin/env python

"""Run VASP job on slab.

Author(s): Michal Badich, Raul A. Flores
"""

# | - Import Modules
import os
import subprocess
import ase.calculators.vasp as vasp_calculator
from ase import io

# My Modules
from ase_modules.ase_methods import clean_up_dft
#__|

# | - Script Inputs
dipole_corr = False
#__|

# | - Read Atoms Object
if os.path.isfile("init.cif"):
    atoms = io.read('init.cif')
elif os.path.isfile("init.traj"):
    atoms = io.read('init.traj')
#__|

# atoms.center(vacuum=7.5, axis=2)

# | - Copy Previous OUTCAR and moments.traj
subprocess.call('cp -rf OUTCAR OUTCAR_$(date +%s)', shell=True)
subprocess.call('cp -rf moments.traj moments.traj_$(date +%s)', shell=True)
#__|

# | - Calculator
calc = vasp_calculator.Vasp(
    encut=500,
    xc='PBE',
    #setups={'O': '_s', 'C': '_s'},
    gga='PE',
    #kpts  = (2,2,1),
    # kpts=(6, 6, 1),
    # For typical 2x2 cell IrO2 110 with 2 cusp sites and 2 bridge sites
    kpts=(4, 4, 4),
    kpar=10,
    npar=4,
    gamma=True,  # Gamma-centered (defaults to Monkhorst-Pack)
    ismear=0,
    inimix=0,

    # amix=0.2,
    # bmix=0.0001,
    # amix_mag=0.1,
    # bmix_mag=0.0001,

    # Conservative Mixing Paramters
    amix=0.05,
    bmix=0.00001,
    amix_mag=0.05,
    bmix_mag=0.00001,


    #nupdown= 0,
    nelm=250,
    sigma=0.05,
    algo='normal',
    ibrion=2,
    isif=2,
    ediffg=-0.02,  # forces
    ediff=1e-5,  # energy conv.
    #nedos=2001,
    prec='Normal',
    nsw=150,  # don't use the VASP internal relaxation, only use ASE
    lvtot=False,
    ispin=2,
    ldau=False,
    ldautype=2,
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

    # lmaxmix=4,
    lmaxmix=6,

    lorbit=11,

    idipol=3,
    dipol=(0, 0, 0.5),
    ldipol=dipole_corr,

    # ldipol=True,

    # addgrid=True,
    # isym=0,
    )

atoms.set_calculator(calc)
#__|

atoms.get_potential_energy()

io.write('out.cif', atoms)
clean_up_dft()
