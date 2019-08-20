from ase.calculators.vasp import Vasp
import ase.calculators.vasp as vasp_calculator
import os
import subprocess


from ase.io import *
from ase.io.trajectory import PickleTrajectory
from ase.visualize import view
from ase.io import write
from ase import Atoms
#from ase.lattice import bulk
from ase.io import read
from ase import Atoms
from ase import Atom
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.optimize import BFGS
from ase.visualize import *
from ase.io import *
import numpy as np
from ase.io.trajectory import Trajectory

name = 'vasp_run2'
#atoms=read('start.cif')
#atoms=read('start.traj')
atoms=read('OUTCAR',  index=-1)
#atoms=read('moments.traj')
#watoms=read('CONTCAR')

#moms=atoms.get_magnetic_moments()
#atoms.set_initial_magnetic_moments(moms)

#atoms[56].magmom=1.0
#atoms[57].magmom=1.0
#atoms[60].magmom=1.0
#atoms[61].magmom=1.0

subprocess.call('cp -rf OUTCAR OUTCAR_$(date +%s)', shell=True)
subprocess.call('cp -rf moments.traj moments.traj_$(date +%s)', shell=True)


#atoms[1].magmom=-3.0

#print cell

#for a in atoms:
#	if(a.symbol=='Fe'):
#		a.magmom=5.0
#	if(a.symbol=='Cr'):
#		a.magmom=3.0
#	print a.magmom
#atoms.center(vacuum=7.5,axis=2)

#del atoms[44]

atoms.write(name+'_init'+'.traj')
atoms.write(name+'_init'+'.cif')

#vasp_calculator.int_keys.append("nedos")
#vasp_calculator.int_keys.append("ICHARG")
#vasp_calculator.float_keys.append("AMIX")
#vasp_calculator.float_keys.append("AMIX_MAG")
#vasp_calculator.float_keys.append("BMIX")
#vasp_calculator.float_keys.append("BMIX_MAG")
#vasp_calculator.int_keys.append("INIMIX")
#vasp_calculator.bool_keys.append("LASPH")
#vasp_calculator.string_keys.append("LREAL")
#vasp_calculator.int_keys.append("kpar")
#vasp_calculator.int_keys.append("ncore")
#vasp_calculator.int_keys.append("NUPDOWN")

unitcell=atoms.get_cell()
print unitcell
kpoints=[]
for i in [0,1,2]:
	l=np.sqrt(unitcell[i][0]**2+unitcell[i][1]**2+unitcell[i][2]**2)
	k=int(4*5/l)
	if(k>0):
		kpoints.append(k)
	else:
		kpoints.append(1)

print kpoints
kpar=10
if kpoints[0]*kpoints[1]*kpoints[2]<10 :
        if(kpoints[0]*kpoints[1]*kpoints[2] % 2 == 0):
                kpar=2
        else:
                kpar=1

print kpar

calc = vasp_calculator.Vasp(encut=600,
			xc='PBE',
                        #setups={'O': '_s', 'C': '_s'},
                        gga='PE',
			#kpts  = (2,2,1),
			kpts  = (kpoints[0],kpoints[1],kpoints[2]),
                        kpar = kpar,
			npar= 4,
			gamma = True, # Gamma-centered (defaults to Monkhorst-Pack)
			ismear=0,
			inimix = 0,
			amix  = 0.2,
			bmix   = 0.0001,
			amix_mag= 0.1,
			bmix_mag= 0.0001,
                        #nupdown= 0,
                        nelm=120,	
            		sigma = 0.05,
            		algo = 'normal',
			ibrion=2,
                        isif=3,
			ediffg=-0.02,  # forces
			ediff=1e-5,  #energy conv.
                        #nedos=2001,
			prec='Normal',
			nsw=150, # don't use the VASP internal relaxation, only use ASE
			lvtot=False,
			ispin=2,
                        ldau=False, 
                        ldautype=2,
                        #lreal  = 'False',#more accuratre volume
                        lreal  = 'auto',
			lasph =True, 
             		ldau_luj={'Ni':{'L':2,  'U':6.45, 'J':0.0},
				'Co':{'L':2,  'U':3.32, 'J':0.0},
                                'Cr':{'L':2,  'U':3.5,  'J':0.0},   
				'Fe':{'L':2,  'U':5.3, 'J':0.0},
                                'Ce':{'L':3,'U':4.50,'J':0.0},
				'V':{'L':2,  'U':3.25, 'J':0.0},
				'Mn':{'L':2,  'U':3.75, 'J':0.0},
				'Ti':{'L':2,  'U':3.00, 'J':0.0},
				'W':{'L':-1,  'U':0.0, 'J':0.0},
                        	'O':{'L':-1, 'U':0.0, 'J':0.0},
                        	'C':{'L':-1, 'U':0.0, 'J':0.0},
				'Au':{'L':-1, 'U':0.0, 'J':0.0},
				'Cu':{'L':-1, 'U':0.0, 'J':0.0},
				'H':{'L':-1, 'U':0.0, 'J':0.0}},
             		ldauprint=2,
                        #lmaxmix=4,
                        lmaxmix=6,
                        lorbit=11,
			#idipol=3,
                        #dipol=(0, 0, 0.5),
                        #ldipol=True
			)

atoms.set_calculator(calc)
atoms.get_potential_energy()

#dyn = QuasiNewton(atoms, logfile=name+'.log', trajectory=name+'.traj')
#dyn.run(fmax=0.05)

write('final.traj',atoms)



import os
import subprocess
#subprocess.call('rm -f WAVECAR', shell=True)
#write(name+'opt.traj',atoms)

