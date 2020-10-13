# #########################################################
#
# #########################################################

#| - Import Modules
import os

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
calc = vasp_calculator.Vasp(
    encut=500,
    xc='PBE',
    gga='PS',
    kpts  = (3,3,3),
    kpar = 10,
    npar = 4,
    gamma = True, # Gamma-centered (defaults to Monkhorst-Pack)
    #ICHARG =1 ,
    #ISTART = 2,
    #isym = 0,
    ismear=0,
    inimix = 0,
    amix   = 0.1,
    bmix   = 0.00001,
    amix_mag= 0.1,
    bmix_mag= 0.00001,
    #NUPDOWN= 0,
    nelm=250,
    sigma = 0.05,
    algo = 'normal',
    #algo = 'fast',
    ibrion=2,
    isif=2,
    ediffg=-0.02,  # forces
    ediff=1e-5,  #energy conv.
    nedos=2001,
    prec='Normal',
    #prec='Accurate',
    #nsw=130, # don't use the VASP internal relaxation, only use ASE
    nsw=0, # don't use the VASP internal relaxation, only use ASE
    lvtot=False,
    ispin=1,
    ldau=False,
    ldautype=2,
    lreal  = False,
    lwave= True,
    lasph  = True,
    #ldipol  = True,
    #idipol  = 3,
    ldau_luj={'Ni':{'L':2,  'U':6.45, 'J':0.0},
    #ldau_luj={'Ni':{'L':2,  'U':8, 'J':0.0},
    'Co':{'L':2,  'U':3.32, 'J':0.0},
    'Fe':{'L':2,  'U':5.3, 'J':0.0},
    'Ce':{'L':3,  'U':5.5, 'J':1.0},
    'Sm':{'L':3,  'U':4.5, 'J':0.0},
    'W':{'L':-1,  'U':0.0, 'J':0.0},
    'O':{'L':-1, 'U':0.0, 'J':0.0},
    'C':{'L':-1, 'U':0.0, 'J':0.0},
    'Au':{'L':-1, 'U':0.0, 'J':0.0},
    'Cu':{'L':-1, 'U':0.0, 'J':0.0},
    'H':{'L':-1, 'U':0.0, 'J':0.0}},
    ldauprint=2,
    lmaxmix=6,
    lorbit=11,
    )



calc = vasp_calculator.Vasp(
    encut=500,
    xc='PBE',
    gga='PS',
    kpts  = (3,3,3),
    kpar = 5,
    npar = 4,
    gamma = True, # Gamma-centered (defaults to Monkhorst-Pack)
    #ICHARG =1 ,
    #ISTART = 2,
    #isym = 0,
    ismear=0,
    inimix = 0,
    amix   = 0.1,
    bmix   = 0.00001,
    amix_mag= 0.1,
    bmix_mag= 0.00001,
    #NUPDOWN= 0,
    nelm=250,
    sigma = 0.05,
    algo = 'normal',
    #algo = 'fast',
    ibrion=2,
    isif=2,
    ediffg=-0.02,  # forces
    ediff=1e-5,  #energy conv.
    nedos=2001,
    prec='Normal',
    #prec='Accurate',
    #nsw=130, # don't use the VASP internal relaxation, only use ASE
    nsw=0, # don't use the VASP internal relaxation, only use ASE
    lvtot=False,
    ispin=1,
    ldau=False,
    ldautype=2,
    lreal  = False,
    lwave= False,
    lasph  = True,
    #ldipol  = True,
    #idipol  = 3,
    ldau_luj={'Ni':{'L':2,  'U':6.45, 'J':0.0},
    #ldau_luj={'Ni':{'L':2,  'U':8, 'J':0.0},
    'Co':{'L':2,  'U':3.32, 'J':0.0},
    'Fe':{'L':2,  'U':5.3, 'J':0.0},
    'Ce':{'L':3,  'U':5.5, 'J':1.0},
    'Sm':{'L':3,  'U':4.5, 'J':0.0},
    'W':{'L':-1,  'U':0.0, 'J':0.0},
    'O':{'L':-1, 'U':0.0, 'J':0.0},
    'C':{'L':-1, 'U':0.0, 'J':0.0},
    'Au':{'L':-1, 'U':0.0, 'J':0.0},
    'Cu':{'L':-1, 'U':0.0, 'J':0.0},
    'H':{'L':-1, 'U':0.0, 'J':0.0}},
    ldauprint=2,
    lmaxmix=6,
    lorbit=11,
    )
#__|

#| - Phonon Calculation
from ase.phonons import Phonons
ph = Phonons(new_atoms, calc, supercell=(1, 1, 1))
#ph.run()
ph.read(method='frederiksen', acoustic=True)

phonon_energies, phonon_DOS = ph.dos(kpts=(40, 40, 40), npts=3000,
                                     delta=5e-4)

# Calculate the Helmholtz free energy
potentialenergy=0.0
from ase.thermochemistry import CrystalThermo
thermo = CrystalThermo(
    phonon_energies=phonon_energies,
    phonon_DOS=phonon_DOS,
    potentialenergy=potentialenergy,
    formula_units=1)
F = thermo.get_helmholtz_energy(temperature=298.15)


#dyn = QuasiNewton(atoms, logfile=name+'.log', trajectory=name+'.traj')
#dyn.run(fmax=0.05)

io.write('final.traj', new_atoms)
#__|




#| - __old__


# from ase.lattice.spacegroup import crystal
# from ase.io import *
#
# from sys import path
# #path.insert(0,'/nfs/slac/g/suncatfs/vossj/espressopre2')
#
# #if not (os.path.exists('vdw_kernel.bindat')):
# #		os.symlink('/nfs/slac/g/suncatfs/sw/vasp/vdw_kernel.bindat','vdw_kernel.bindat')
#
# #os.environ['ESP_PSP_PATH']='/nfs/slac/g/suncatfs/bajdich/XAS/gold111/customPPs'
# from ase.calculators.vasp import Vasp
#
# from ase.lattice.spacegroup import crystal
# from ase.io import *
# from ase.io.trajectory import PickleTrajectory
# from ase.visualize import view
# from ase.io import write
# from ase import Atoms
# from ase.io import read
# from ase import Atoms
# from ase import Atom
# from ase.constraints import FixAtoms
# from ase.optimize import QuasiNewton
# from ase.optimize import BFGS
# from ase.visualize import *
# from ase.io import *
# import numpy as np
#
# import os
# import subprocess



# #new_atoms=read('../OUTCAR')
# #new_atoms=read('C.json')
# #  new_atoms=read('optimized.json')
# new_atoms=read("init.traj")
# #from ase.build import bulk
# #new_atoms= bulk('Pd', 'fcc', a=4.0, cubic=True)
#
# #moms=new_atoms.get_magnetic_moments()
# #new_atoms.set_initial_magnetic_moments(moms)
#
# #new_atoms.center(vacuum=7.5, axis=2)
# new_atoms.write(name+'_init'+'.traj')
# new_atoms.write(name+'_init'+'.cif')
#
# #print cell


# import os
# import subprocess
# #subprocess.call('rm -f WAVECAR', shell=True)
# #write(name+'opt.traj',atoms)


##LSF -q suncat-test -n 8 -o out.log -e error.txt
##LSF -q suncat3 -n 16 -W 4:00 -o out.log -e error.txt
##LSF -q suncat2 -n 24 -W 4:00 -o out.log -e error.txt
##LSF -q suncat -n 64 -W 48:00  -o out.log -e error.txt
##LSF -q suncat -n 64 -W 48:00 -R "select[hname!=suncat0212 && hname!= suncat0232 && hname!= suncat0184]" -o out.log -e error.txt
##LSF -q psfehidleq  -n 64 -W 48:00 -o out.log -e error.txt
#LSF -q psnehidleq  -n 64 -W 48:00 -o out.log -e error.txt
##LSF -q suncat2 -n 24 -W 48:00 -o out.log -e error.txt
##LSF -q suncat-long -n 32 -o out.log -e error.txt
##LSF -q suncat-long -n 40 -o out.log -e error.txt
##LSF -q suncat3-long -n 32 -o out.log -e error.txt
##LSF -q psanaidleq -n 48 -W 48:00 -o out.log -e error.txt
##LSF -q psfehidleq -n 64 -W 48:00 -o out.log -e error.txt
##LSF -q psnehidleq -n 64 -W 48:00 -o out.log -e error.txt
##LSF -q psanaidleq -n 48 -W 48:00 -o out.log -e error.txt
#! /usr/bin/envpython


#__|
