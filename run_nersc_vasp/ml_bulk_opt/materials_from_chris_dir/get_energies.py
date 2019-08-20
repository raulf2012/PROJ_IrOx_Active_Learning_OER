#!/global/common/cori/software/python/2.7-anaconda/bin/python  
from ase.io import *
	
from sys import path
#path.insert(0,'/nfs/slac/g/suncatfs/vossj/espressopre2')

import os
#if not (os.path.exists('vdw_kernel.bindat')):
#		os.symlink('/nfs/slac/g/suncatfs/sw/vasp/vdw_kernel.bindat','vdw_kernel.bindat')

#os.environ['ESP_PSP_PATH']='/nfs/slac/g/suncatfs/bajdich/XAS/gold111/customPPs'
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
import shutil

from itertools import combinations_with_replacement
from itertools import product

import os
alphabets=os.listdir("../all_cif")
print alphabets

complete_result=[]
iter =0

from operator import itemgetter
indices, L_sorted = zip(*sorted(enumerate((int(a.partition("_")[0])) for a in alphabets), key=itemgetter(1)))
print list(indices)


for i in list(indices):
	a=alphabets[i]
        main_string=a
	print int(a.partition("_")[0]) 
        os.chdir(main_string)
        if os.path.exists("./OUTCAR"):
        	energy=float(subprocess.check_output("grep py= OUTCAR | tail -1 |  gawk 'BEGIN{i=0;a[10]}{a[i]=$7; i=i+1;}END{print a[0]}'", shell=True))
        	energy_before_last=float(subprocess.check_output("grep py= OUTCAR | tail -2 |  gawk 'BEGIN{i=0;a[10]}{a[i]=$7; i=i+1;}END{print a[0]}'", shell=True))
        	natoms=int(subprocess.check_output("grep NIONS OUTCAR | awk '{print $12}'", shell=True))
        	force=float(subprocess.check_output("/project/projectdirs/m2997/bin/forces | grep 'largest force' | gawk '{print $3}'   ", shell=True))
        	pressure=float(subprocess.check_output("grep 'external pressure ' OUTCAR | tail -1 |  gawk 'BEGIN{i=0;a[10]}{a[i]=$4; i=i+1;}END{print a[0]}'", shell=True))
		#if (float(force)< 1 and abs(float(pressure)) < 1):
        	subprocess.check_output("rm -rf WAVECAR", shell=True)
		print 'Structure: ', a , '\n E/atom= ', energy/natoms, '\n E= ',energy, '\n E(N-1)= ',energy_before_last,'\n #atoms= ',natoms,'\n Max Force= ',force,'\n Pressure= ',pressure
		resubmit=False
		#if(abs(energy_before_last-energy)>0.2):
		#	print 'RESUBMIT (EN-1<EN)'
		#	resubmit=True
		#if(force> 3 ):
		#	print 'RESUBMIT (force too large)'
		#	resubmit=True
		#if(abs(pressure)>100):
		#	print 'RESUBMIT (pressure too large)'
		#	resubmit=True
		if(resubmit):
			print 'REDO'
        		#subprocess.call("sbatch ../submit_cori_fresh.pbs", shell=True)
		else:
			complete_result.append([int(a.partition("_")[0]), a, energy/natoms, energy, natoms,force, pressure])
        		#subprocess.call("sbatch ../submit_cori_n-1.pbs", shell=True)
			subprocess.call("cp optimized.cif ../opt_cifs/"+a, shell=True)
			print 'DONE'

		
        os.chdir('../')
	
        iter=iter+1
     
f=open('DFT_energies.dat', 'w')
f.write( '# number, cif file, energy/atom, max force, pressure n #atoms \n')
for k in sorted(complete_result, key=lambda student: student[0]):
	f.write(' %s %s %.6f %.6f %.6f  %i \n' % (k[0], k[1], k[2], k[5],k[6],k[4]))

f.close()
