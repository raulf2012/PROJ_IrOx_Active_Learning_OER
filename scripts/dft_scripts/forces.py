#!/global/common/edison/software/python/2.7-anaconda-4.4/bin/python
##!/usr/common/software/python/2.7.9/bin/python

# | - Import Modules
from math import pow
#from pylab import *
#import matplotlib
#import matplotlib.pyplot as plt
import numpy as np
#from matplotlib.path import Path
#from matplotlib.patches import PathPatch
#from matplotlib.ticker import *
#from scipy import *
import os
import sys
import subprocess
from array import array
import numpy as np
import matplotlib.colors as colors
#from ase.structure import bulk
#from ase.lattice.spacegroup import crystal
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
#__|

if(len(sys.argv)>1):
	file=sys.argv[1]
else:
	file='OUTCAR'

if(len(sys.argv)>2):
	element=sys.argv[2]

atoms=read(file)

energy=subprocess.check_output("grep py= OUTCAR | tail -1 | gawk '{print $7}' ", shell=True)
energy=float(energy)
print(energy)

traj2=Trajectory('moments.traj',  'w')
traj2.write(atoms, energy=energy)
write('optimized.cif',atoms)
moments=[]
moments2=[]
moms=atoms.get_magnetic_moments()
forces=atoms.get_forces()
sum=0.0
largest=0.0
for a in range(len(atoms)):
#	#if(atoms[a].symbol=='Co'):
#	#	moments.append(moms[a])
	if(np.abs(moms[a])>0.1):
		moments.append([atoms[a].symbol,a,moms[a]])
	if(len(sys.argv)>2):
		if(atoms[a].symbol==element):
			moments2.append(moms[a])
	force=np.sqrt(forces[a][0]**2+forces[a][1]**2+forces[a][2]**2)
	sum+=force
	if(force>largest):
		largest=force



if(len(sys.argv)>2):
	print(moments2)
else:
	print('Moments from file %s' %(file))

	for i in range(len(moments)):
		print(moments[i],',')

	print('Forces')
	print(forces)
	print('largest force '+str(largest)+'  '+str(sum) )
