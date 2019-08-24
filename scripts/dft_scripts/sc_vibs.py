## #! /usr/bin/envpython

#| - Import Modules
from ase.structure import bulk
from ase.lattice.spacegroup import crystal
from ase.io import *

from sys import path
#path.insert(0,'/nfs/slac/g/suncatfs/vossj/espressopre2')

import os
#if not (os.path.exists('vdw_kernel.bindat')):
#        os.symlink('/nfs/slac/g/suncatfs/sw/vasp/vdw_kernel.bindat','vdw_kernel.bindat')

#os.environ['ESP_PSP_PATH']='/nfs/slac/g/suncatfs/bajdich/XAS/gold111/customPPs'
from ase.calculators.vasp import Vasp
import ase.calculators.vasp as vasp_calculator

from ase.structure import bulk
from ase.lattice.spacegroup import crystal
from ase.io import *
from ase.io.trajectory import PickleTrajectory
from ase.visualize import view
from ase.io import write
from ase import Atoms
from ase.lattice import bulk
from ase.io import read
from ase import Atoms
from ase import Atom
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.optimize import BFGS
from ase.visualize import *
from ase.io import *
import numpy as np

import os
import subprocess
#__|

name = 'vasp_run2'
#new_atoms=read('../OUTCAR')
new_atoms=read('restart.traj')
#new_atoms=read('restart.traj')
#new_atoms=read('CONTCAR')

#moms=new_atoms.get_magnetic_moments()
#for i, m in enumerate(moms):
#        new_atoms[i].magmom=m

#new_atoms.center(vacuum=5.5, axis=2)
new_atoms.write(name+'_init'+'.traj')
new_atoms.write(name+'_init'+'.cif')

#print cell
#| - Calculator
vasp_calculator.int_keys.append("nedos")
vasp_calculator.int_keys.append("ICHARG")
vasp_calculator.float_keys.append("AMIX")
vasp_calculator.float_keys.append("AMIX_MAG")
vasp_calculator.float_keys.append("BMIX")
vasp_calculator.float_keys.append("BMIX_MAG")
vasp_calculator.int_keys.append("INIMIX")
vasp_calculator.bool_keys.append("LASPH")
vasp_calculator.bool_keys.append("LWAVE")
vasp_calculator.bool_keys.append("LDIPOL")
vasp_calculator.string_keys.append("LREAL")
vasp_calculator.int_keys.append("kpar")
vasp_calculator.int_keys.append("ncore")
vasp_calculator.int_keys.append("NUPDOWN")
vasp_calculator.int_keys.append("IDIPOL")

calc = vasp_calculator.Vasp(
    encut=500,
    xc='PBE',
    #setups={'O': '_s', 'C': '_s'},
    gga='PE',
    kpts  = (5,5,1),
    kpar = 5,
    npar = 8,
    gamma = True, # Gamma-centered (defaults to Monkhorst-Pack)
    ismear=0,
    INIMIX = 0,
    AMIX   = 0.2,
    BMIX   = 0.0001,
    AMIX_MAG= 0.4,
    BMIX_MAG= 0.0001,
    #NUPDOWN= 0,
    nelm=250,
    sigma = 0.05,
    algo = 'normal',
    ibrion=2,
    ediffg=-0.02,  # forces
    ediff=1e-6,  #energy conv.
    nedos=2001,
    prec='Normal',
    #prec='Accurate',
    nsw=0, # don't use the VASP internal relaxation, only use ASE
    #nsw=200, # don't use the VASP internal relaxation, only use ASE
    lvtot=False,
    ispin=2,
    ldau=True,
    ldautype=2,
    LREAL  = 'auto',
    LASPH  = True,
    LWAVE  = False,
    LDIPOL  = True,
    IDIPOL  = 3,
    dipol=(0, 0, 0),
    #ldau_luj={'Ni':{'L':2,  'U':6.45, 'J':0.0},
    ldau_luj={'Ni':{'L':-1,  'U':0.0, 'J':0.0},
    'Co':{'L':2,  'U':3.32, 'J':0.0},
    'Fe':{'L':2,  'U':5.3, 'J':0.0},
    'Ce':{'L':3,  'U':5.5, 'J':1.0},
    'W':{'L':-1,  'U':0.0, 'J':0.0},
    'O':{'L':-1, 'U':0.0, 'J':0.0},
    'C':{'L':-1, 'U':0.0, 'J':0.0},
    'Au':{'L':-1, 'U':0.0, 'J':0.0},
    'Cu':{'L':-1, 'U':0.0, 'J':0.0},
    'H':{'L':-1, 'U':0.0, 'J':0.0}},
    ldauprint=2,
    lmaxmix=6,
    lorbit=11)
#__|

new_atoms.set_calculator(calc)
#new_atoms.get_potential_energy()

from ase.vibrations import Vibrations

electronicenergy = 0.0# new_atoms.get_potential_energy()
write('final.traj',new_atoms)
#write(name+'opt.traj',atoms)
vibstr = 'vib_'+name
#vib = Vibrations(new_atoms, name=vibstr, indices=[87,88,85])
vib = Vibrations(new_atoms, name=vibstr, indices=[34],delta=0.005)
#vib = Vibrations(new_atoms, name=vibstr, indices=[88])
#vib = Vibrations(new_atoms, name=vibstr, indices=[87])

vib.run()
vib.summary(method='frederiksen')

# Make trajectory files to visualize normal modes:
for mode in range(len(vib.modes)):
    vib.write_mode(mode)


##vib.PrintHessianMatrix(mass=1,precision=3,suppress_small=1)
#vib.PrintFrequencies()
#print 'Zero-point energy = %1.2f eV' % vib.GetZeroPointEnergy()

#from ase.thermochemistry import IdealGasThermo
from ase.thermochemistry import HarmonicThermo
vib_energies = vib.get_energies()

print vib_energies

thermo = HarmonicThermo(vib_energies=vib_energies,
                        electronicenergy=electronicenergy,
                        )
G = thermo.get_internal_energy(temperature=298.15)

#import os
#import subprocess
#subprocess.call('rm -f WAVECAR', shell=True)
#write(name+'opt.traj',atoms)
