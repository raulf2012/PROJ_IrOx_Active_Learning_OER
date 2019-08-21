import numpy as np
from ase import io
import re
import pandas as pd
import os


""" 
    rapiDOS v.0.5.1 
    -- Feb. - 19 - 2018 -- 
    Jose A. Garrido Torres & Michal Bajdich. 
    jagt@stanford.edu  &  bajdich@slac.stanford.edu 
"""

##############################################################################
# Get general information of the system:
##############################################################################

atoms = io.read('CONTCAR') # Open CONTCAR file
range_of_atoms = range(len(atoms)) # Number of atoms.
atomic_numbers = atoms.get_atomic_numbers() # Atomic number of the atoms.
atomic_symbols = atoms.get_chemical_symbols() # Chemical symbol of the atoms.
outcar_file = io.read('OUTCAR', format='vasp-out')

##############################################################################


##############################################################################
# Check spin:
##############################################################################

incar_file = open("INCAR", "r")
ispin = 1  # Non spin polarised calculations.
for line in incar_file:
    if re.match("(.*)ISPIN(.*)2", line):
        ispin = 2 # For spin polarised calculations.

##############################################################################


##############################################################################
# Check scale:
##############################################################################

# KBLOCK scales DOSCAR (See VASP Manual).
for line in open('OUTCAR'):
    if line.find('KBLOCK') != -1:
        ckblock = line.split()
        for i in range(0, len(ckblock)):
            if ckblock[i] == 'KBLOCK': kblock = float(ckblock[i+2])

##############################################################################


##############################################################################
# Code adapted from split_dos.py which belongs to Jonsson and Henkelman groups.
# URL: http://theory.cm.utexas.edu/vtsttools/scripts.html
##############################################################################

### READ DOSCAR ###
def read_dosfile():
    f = open("DOSCAR", 'r')
    lines = f.readlines()
    f.close()
    index = 0
    natoms = int(lines[index].strip().split()[0])
    index = 5
    nedos = int(lines[index].strip().split()[2])
    efermi = float(lines[index].strip().split()[3])
    print(natoms, nedos, efermi)

    return lines, index, natoms, nedos, efermi

### READ POSCAR or CONTCAR and save pos
def read_posfile():
    from ase.io import read

    try:
        atoms = read('POSCAR')
    except IOError:
        print("[__main__]: Couldn't open input file POSCAR, atomic positions will not be written...\n")
        atoms = []

    return atoms

### WRITE DOS0 CONTAINING TOTAL DOS ###
def write_dos0(lines, index, nedos, efermi):

    fdos = open("DOS0", 'w')
    line = lines[index+1].strip().split()
    ncols = int(len(line))
    # fdos.write('# %d \n' % (ncols)) #GH not sure why this is here

    for n in range(nedos):
        index +=1
        e = float(lines[index].strip().split()[0])
        e_f = e-efermi
        fdos.write('%15.8f ' % (e_f))

        for col in range(1, ncols):
            dos = float(lines[index].strip().split()[col])
            fdos.write('%15.8f ' % (dos))
        fdos.write('\n')
    return index

### LOOP OVER SETS OF DOS, NATOMS ###
def write_nospin(lines, index, nedos, natoms, ncols, efermi):
    
    atoms = read_posfile()
    if len(atoms) < natoms:
    	pos = np.zeros((natoms, 3))
    else:
        pos = atoms.get_positions()

    for i in range(1,natoms+1):
        si = str(i)

    ## OPEN DOSi FOR WRITING ##
        fdos = open("DOS"+si, 'w')
        index += 1
        ia = i-1
        # fdos.write('# %d \n' % (ncols))
        fdos.write('# %15.8f %15.8f %15.8f \n' % (pos[ia,0], pos[ia,1], pos[ia,2]))

    ### LOOP OVER NEDOS ###
        for n in range(nedos):
            index += 1
            e = float(lines[index].strip().split()[0])
            e_f = e-efermi
            fdos.write('%15.8f ' % (e_f))

            for col in range(1, ncols):
                dos = float(lines[index].strip().split()[col])
                fdos.write('%15.8f ' % (dos))
            fdos.write('\n')
    fdos.close()

def write_spin(lines, index, nedos, natoms, ncols, efermi):
    #pos=[]
    atoms = read_posfile()
    if len(atoms) < natoms:
        pos = np.zeros((natoms, 3))
    else:
        pos = atoms.get_positions()

    nsites = (ncols -1)/2

    for i in range(1,natoms+1):
        si = str(i)
    ## OPEN DOSi FOR WRITING ##
        fdos = open("DOS"+si, 'w')
        index += 1
        ia = i-1
        fdos.write('# %d \n' % (ncols))
        fdos.write('# %15.8f %15.8f %15.8f \n' % (pos[ia,0], pos[ia,1], pos[ia,2]))

    ### LOOP OVER NEDOS ###
        for n in range(nedos):
            index +=1   
            e = float(lines[index].strip().split()[0])
            e_f = e-efermi
            fdos.write('%15.8f ' % (e_f))

            for site in range(int(nsites)):
                dos_up = float(lines[index].strip().split()[site*2+1])
                dos_down = float(lines[index].strip().split()[site*2+2])*-1
                fdos.write('%15.8f %15.8f ' % (dos_up, dos_down))
            fdos.write('\n')
        fdos.close()

#
if __name__ == '__main__':
    import sys
    import os
    import datetime
    import time
    import optparse

    lines, index, natoms, nedos, efermi = read_dosfile()
    index = write_dos0(lines, index, nedos, efermi)
    ## Test if a spin polarized calculation was performed ##
    line = lines[index+2].strip().split()
    ncols = int(len(line)) 
    if ncols==7 or ncols==19 or ncols==9 or ncols==33:
        write_spin(lines, index, nedos, natoms, ncols, efermi)
        is_spin=True
    else: 
        write_nospin(lines, index, nedos, natoms, ncols, efermi)
        is_spin=False
    print("Spin unrestricted calculation: ", is_spin)
    print(ncols)


##############################################################################


##############################################################################
# Define the columns of the database for the PDOS:
##############################################################################

if ispin == 2:
    if ncols < 19:
        pdos_columns = ['Energy (E-Ef)', 's_up','s_down', 'py_up','py_down',
        'pz_up', 'pz_down', 'px_up', 'px_down', 'Element', 'Atom Number',
        'Atom Label']
    if ncols == 19: 
        pdos_columns = ['Energy (E-Ef)', 's_up', 's_down', 'py_up',
        'py_down', 'pz_up', 'pz_down', 'px_up', 'px_down', 'dxy_up',
        'dxy_down', 'dyz_up', 'dyz_down', 'dz2_up', 'dz2_down', 'dxz_up',
        'dxz_down', 'dx2_up', 'dx2_down', 'Element', 'Atom Number',
        'Atom Label']
    if ncols > 19: 
        pdos_columns = ['Energy (E-Ef)', 's_up', 's_down', 'py_up', 'py_down',
        'pz_up', 'pz_down', 'px_up', 'px_down', 'dxy_up', 'dxy_down',
        'dyz_up', 'dyz_down', 'dz2_up', 'dz2_down', 'dxz_up', 'dxz_down',
        'dx2_up', 'dx2_down', 'f1_up', 'f1_down', 'f2_up', 'f2_down', 'f3_up',
        'f3_down', 'f4_up', 'f4_down', 'f5_up', 'f5_down', 'f6_up',
        'f6_down', 'Element', 'Atom Number', 'Atom Label']

if ispin != 2:
    if ncols < 10 : 
        pdos_columns = ['Energy (E-Ef)', 's_up', 'py_up', 'pz_up', 'px_up', 'Atom Label']
    if ncols == 10:
        pdos_columns = ['Energy (E-Ef)', 's_up', 'py_up', 'pz_up', 'px_up', 'dxy_up', 'dyz_up',
        'dz2_up', 'dxz_up', 'dx2_up', 'Element', 'Atom Number', 'Atom Label']
    if ncols > 10: 
        pdos_columns = ['Energy (E-Ef)', 's_up', 'py_up', 'pz_up', 'px_up', 'dxy_up', 'dyz_up',
        'dz2_up', 'dxz_up', 'dx2_up', 'f1_up', 'f2_up', 'f3_up', 'f4_up', 'f5_up', 'f6_up', 'Element',
        'Atom Number', 'Atom Label']

##############################################################################


##############################################################################
# Build pandas dataframe for total DOS:
##############################################################################

total_dos = np.loadtxt('DOS0', skiprows=0) # Load file total DOS from DOS0.
total_dos[:,1:] = total_dos[:,1:]* kblock # Scale (See kblock VASP).


if ispin == 2:
    total_dos_columns = ['Energy (E-Ef)', 'Total DOS Spin Up', 'Total DOS '
    'Spin '
    'Down', 'Integrated DOS Spin Up', 'Integrated DOS Spin Down']
if ispin != 2:
    total_dos_columns = ['Energy (E-Ef)', 'Total DOS Spin Up', 'Integrated DOS']

total_dos_df = pd.DataFrame(total_dos, columns=total_dos_columns)

##############################################################################


##############################################################################
# Build pandas dataframe for PDOS:
##############################################################################

pdos_df = pd.DataFrame([], columns=pdos_columns)

for i in range_of_atoms:
    dos_atom_i = np.loadtxt('DOS'+str(i+1), skiprows=0) # Load the DOSxx files.
    dos_atom_i[:,1:] = dos_atom_i[:,1:] # Do NOT scale PDOS (See KBLOCK VASP).
    index_i = str(atomic_symbols[i]) # Element
    index_i2 = str(i+1) # Atom Number.
    index_i3 = str(atomic_symbols[i]+str(i+1)) # Element index.
    column_index_i = np.repeat(index_i, len(dos_atom_i))
    column_index_i2 = np.repeat(index_i2, len(dos_atom_i))
    column_index_i3 = np.repeat(index_i3, len(dos_atom_i))
    column_index_i = np.reshape(column_index_i, (len(column_index_i),1))
    column_index_i2 = np.reshape(column_index_i2, (len(column_index_i2),1))
    column_index_i3 = np.reshape(column_index_i3, (len(column_index_i3),1))
    dos_atom_i = np.append(dos_atom_i, column_index_i,axis=1) # Append Element.
    dos_atom_i = np.append(dos_atom_i, column_index_i2,axis=1) # Append Atom N.
    dos_atom_i = np.append(dos_atom_i, column_index_i3,axis=1) # Append Label.
    pdos_df_i = pd.DataFrame(dos_atom_i, columns=pdos_columns)
    pdos_df = pd.DataFrame.append(pdos_df_i, pdos_df, ignore_index=True)

##############################################################################


##############################################################################
# Get band gap:
##############################################################################

def get_bandgap():
    arg_fermi = np.where(total_dos[:, 0] > 0)[0][0]
    arg_fermi_upper = arg_fermi
    arg_fermi_lower = arg_fermi
    band_tolerance = 5e-2

    converged = False
    while not converged:
        occup_lower = total_dos[:, 1][arg_fermi_lower]
        occup_upper = total_dos[:, 1][arg_fermi_upper]

        if occup_lower > band_tolerance:
            arg_fermi_lower = arg_fermi_lower + 1
        if occup_upper > band_tolerance:
            arg_fermi_upper = arg_fermi_upper - 1

        if occup_lower > band_tolerance and occup_upper > band_tolerance:
            converged = True
            e_lower = total_dos[arg_fermi_lower][0]
            e_upper = total_dos[arg_fermi_upper][0]
            band_gap = e_lower - e_upper
        arg_fermi_lower = arg_fermi_lower - 1
        arg_fermi_upper = arg_fermi_upper + 1

    print("Approx. band gap: ", np.abs(band_gap), "eVv")
    return [e_lower, e_upper, band_gap]

band_gap = get_bandgap()
band_gap_columns = ['Lower Band Gap', 'Upper Band Gap', 'Band Gap']
band_gap_df = pd.DataFrame([band_gap], columns=band_gap_columns)

##############################################################################


##############################################################################
# Save files for Jupyter notebook
##############################################################################

total_dos_df.to_csv('TotalDOS.csv')
pdos_df.to_csv('PDOS.csv')
band_gap_df.to_csv('BandGap.csv')

##############################################################################

# Open example file (rapiDOS_analysis.ipynb) using Jupyter Notebook.

os.system('jupyter notebook rapiDOS_analysis.ipynb')
#use the second to run it without a browser
#os.system('jupyter nbconvert --to notebook --execute rapiDOS_analysis.ipynb')
