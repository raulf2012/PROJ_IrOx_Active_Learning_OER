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
from ase.lattice.spacegroup import crystal
from ase.io import *
	
from sys import path
#path.insert(0,'/nfs/slac/g/suncatfs/vossj/espressopre2')

import os
#if not (os.path.exists('vdw_kernel.bindat')):
#		os.symlink('/nfs/slac/g/suncatfs/sw/vasp/vdw_kernel.bindat','vdw_kernel.bindat')

#os.environ['ESP_PSP_PATH']='/nfs/slac/g/suncatfs/bajdich/XAS/gold111/customPPs'
from ase.calculators.vasp import Vasp
import ase.calculators.vasp as vasp_calculator

from ase.lattice.spacegroup import crystal
from ase.io import *
from ase.io.trajectory import PickleTrajectory
from ase.visualize import view
from ase.io import write
from ase import Atoms
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

name = 'vasp_run2'
#new_atoms=read('../OUTCAR')
#new_atoms=read('C.json')
new_atoms=read('restart.json')
#new_atoms=read('CeO2_4x4_ortho_w_epoxy.json')

#moms=new_atoms.get_magnetic_moments()
#new_atoms.set_initial_magnetic_moments(moms)


#new_atoms.center(vacuum=6.5, axis=2)
new_atoms.write(name+'_init'+'.traj')
new_atoms.write(name+'_init'+'.cif')

#print cell



calc = vasp_calculator.Vasp(encut=500,
			xc='PBE',
                        #setups={'O': '_s', 'C': '_s'},
                        gga='PE',
			kpts  = (3,3,1),
                        kpar = 5,
			npar = 6,
			gamma = True, # Gamma-centered (defaults to Monkhorst-Pack)
                        #ICHARG =1 ,
                        #ISTART = 2,
			ismear=0,
			inimix = 0,
			amix   = 0.05,
                        bmix   = 0.00001,
                        amix_mag= 0.05,
                        bmix_mag= 0.00001,
                        #NUPDOWN= 0,
                        nelm=250,	
            		#sigma = 0.05,
            		sigma = 0.10,
            		algo = 'normal',
			ibrion=2,
			ediffg=-0.02,  # forces
			ediff=1e-5,  #energy conv.
                        nedos=2001,
			prec='Normal',
			#prec='Accurate',
			#nsw=150, # don't use the VASP internal relaxation, only use ASE
			nsw=250, # don't use the VASP internal relaxation, only use ASE
			lvtot=False,
			ispin=2,
                        ldau=True, 
                        ldautype=2,
                        lreal  = 'auto',
			lasph  = True, 
			#ldipol  = True, 
			#idipol  = 3, 
             		ldau_luj={'Ni':{'L':2,  'U':6.45, 'J':0.0},
             		#ldau_luj={'Ni':{'L':2,  'U':8, 'J':0.0},
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

new_atoms.set_calculator(calc)
new_atoms.get_potential_energy()

#dyn = QuasiNewton(atoms, logfile=name+'.log', trajectory=name+'.traj')
#dyn.run(fmax=0.05)

write('final.traj',new_atoms)



import os
import subprocess
#subprocess.call('rm -f WAVECAR', shell=True)
subprocess.call('ase convert OUTCAR full_relax.json', shell=True)
subprocess.call('/global/homes/b/bajdich/bin/get_restart3', shell=True)
#write(name+'opt.traj',atoms)

