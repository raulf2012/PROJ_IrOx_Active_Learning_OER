import os
from ase import io

from os import listdir
from os.path import isfile, join


# | - Script Inputs
out_dir = "out_cifs"

#__|

# for root, dirs, files in os.walk("."):
#     tmp = 42

files = [f for f in listdir(".") if isfile(join(".", f))]

if not os.path.exists(out_dir):
    os.mkdir(out_dir)

cif_files = []
for file_i in files:
    if ".cif" in file_i:
        cif_files.append(file_i)


for cif_i in cif_files:
    atoms_i = io.read(cif_i)
    atoms_i.center(vacuum=7.5, axis=2)

    out_dir_i = os.path.join(
        out_dir,
        cif_i,
        )

    print(out_dir_i)

    atoms_i.write(out_dir_i)
