# #########################################################
#
# #########################################################

#| - Import Modules
import os
import sys

import numpy as np

from ase import io
import ase.calculators.vasp as vasp_calculator
#__|

#| - Script Inputs
name = 'vasp_run2'
#__|

#| - Read Atoms
try:
    new_atoms = io.read("init.traj")
except:
    new_atoms = io.read("init.cif")

new_atoms.write(name + "_init.traj")
new_atoms.write(name + "_init.cif")
#__|

# | - k-points
atoms = new_atoms
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

#| - VASP Calculator
shared_calc_settings = dict(
    encut=600,
    xc='PBE',
    gga='PE',
    kpts=kpoints,
    npar=6,
    gamma=True, # Gamma-centered (defaults to Monkhorst-Pack)
    #ICHARG=1 ,
    #ISTART=2,
    #isym=0,

    ismear=0,
    inimix=0,
    amix=0.1,
    bmix=0.00001,
    amix_mag=0.1,
    bmix_mag=0.00001,
    #NUPDOWN=0,
    nelm=250,
    sigma=0.05,
    algo='fast',
    #algo='fast',
    ibrion=2,
    isif=2,
    ediffg=-0.02,  # forces
    ediff=1e-5,  #energy conv.
    nedos=2001,
    # prec='Normal',
    prec='Accurate',
    #nsw=130, # don't use the VASP internal relaxation, only use ASE
    nsw=0, # don't use the VASP internal relaxation, only use ASE
    lvtot=False,
    ispin=2,
    ldau=False,
    ldautype=2,
    lreal=False,
    lasph=True,
    #ldipol=True,
    #idipol=3,
    ldau_luj={'Ni':{'L':2, 'U':6.45, 'J':0.0},},
    #  ldau_luj={
    #      'Ni':{'L':2,  'U':8, 'J':0.0},
    #      'Co':{'L':2,  'U':3.32, 'J':0.0},
    #      'Fe':{'L':2,  'U':5.3, 'J':0.0},
    #      'Ce':{'L':3,  'U':5.5, 'J':1.0},
    #      'Sm':{'L':3,  'U':4.5, 'J':0.0},
    #      'W':{'L':-1,  'U':0.0, 'J':0.0},
    #      'O':{'L':-1, 'U':0.0, 'J':0.0},
    #      'C':{'L':-1, 'U':0.0, 'J':0.0},
    #      'Au':{'L':-1, 'U':0.0, 'J':0.0},
    #      'Cu':{'L':-1, 'U':0.0, 'J':0.0},
    #      'H':{'L':-1, 'U':0.0, 'J':0.0},
    #      },
    ldauprint=2,
    lmaxmix=6,
    lorbit=11,
    )

calc = vasp_calculator.Vasp(
    kpar=10,
    lwave=True,
    **shared_calc_settings)

new_atoms.set_calculator(calc)
new_atoms.get_potential_energy()

calc = vasp_calculator.Vasp(
    kpar=5,
    lwave=False,
    **shared_calc_settings)
#__|


#  assert False

#| - Phonon Calculation
from ase.phonons import Phonons
ph = Phonons(new_atoms, calc, supercell=(1, 1, 1))
ph.run()
ph.read(method='frederiksen', acoustic=True)

#dyn = QuasiNewton(atoms, logfile=name+'.log', trajectory=name+'.traj')
#dyn.run(fmax=0.05)

io.write('final.traj', new_atoms)

#subprocess.call('rm -f WAVECAR', shell=True)
#write(name+'opt.traj',atoms)
#__|
