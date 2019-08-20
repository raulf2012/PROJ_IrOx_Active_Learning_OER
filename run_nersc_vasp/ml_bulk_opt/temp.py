#!/global/common/cori/software/python/2.7-anaconda/bin/python

#| - Import  Modules
from ase.io import *

from sys import path
#path.insert(0,'/nfs/slac/g/suncatfs/vossj/espressopre2')

import os
#if not (os.path.exists('vdw_kernel.bindat')):
#        os.symlink('/nfs/slac/g/suncatfs/sw/vasp/vdw_kernel.bindat','vdw_kernel.bindat')

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
# __|

alphabets=os.listdir("../all_cif")
print alphabets

iter =0


from operator import itemgetter
indices, L_sorted = zip(*sorted(enumerate((int(a.partition("_")[0])) for a in alphabets), key=itemgetter(1)))
print list(indices)


for i in list(indices):
    a=alphabets[i]
    main_string=a
    print iter+1,'',a
    if not os.path.exists(main_string):
        os.mkdir(main_string)
        os.chdir(main_string)
        new_atoms=read('../../all_cif/'+a)
        write('start.traj',new_atoms)
        subprocess.call("sbatch ../submit_edison_fresh.pbs", shell=True)
        os.chdir('../')
    iter=iter+1
