# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Import Modules

# +
import os
import sys
import pickle

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import bulk_dft_data_path

# #############################################################################
from ase_modules.ase_methods import view_in_vesta

# +
from ase_modules.ase_methods import view_in_vesta

import pandas as pd

# +
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling",
    "ccf_similarity_analysis/out_data",
    "all_ids_to_elim_1.pickle")
with open(path_i, "rb") as fle:
    all_ids_to_elim = pickle.load(fle)

with open(bulk_dft_data_path, "rb") as fle:
    df_bulk_dft = pickle.load(fle)

# +
df_bulk_dft = df_bulk_dft[
    (df_bulk_dft["source"] != "oqmd") & \
    (df_bulk_dft["source"] != "raul_oer") & \
    (df_bulk_dft["source"] != "chris") & \
    [True for i in range(len(df_bulk_dft))]
    ]

df_bulk_dft = df_bulk_dft.drop(all_ids_to_elim)

# + active=""
#
#
#
#
#
#
#
#

# +
df_bulk_dft = df_bulk_dft.sort_values("energy_pa")

# row_i = df_bulk_dft.iloc[2]
row_i = df_bulk_dft.loc["cqbrnhbacg"]


atoms = row_i["atoms"]
atoms

# view_in_vesta(atoms)

# +
from pymatgen.analysis import local_env
from pymatgen.io.ase import AseAtomsAdaptor

structure = AseAtomsAdaptor.get_structure(atoms)


CrysNN = local_env.CrystalNN(
    weighted_cn=False,
    cation_anion=False,
    distance_cutoffs=(0.5, 1),
    x_diff_weight=3.0,
    porous_adjustment=True,
    search_cutoff=7,
    fingerprint_length=None)


coord_data_dict = {
    # "": ,
    }

data_master = []
for i_cnt, site_i in enumerate(structure.sites):
    site_elem_i = site_i.species_string

    data_dict_i = dict()

    data_dict_i["element"] = site_elem_i
    data_dict_i["structure_index"] = i_cnt

    nn_info_i = CrysNN.get_nn_info(structure, i_cnt)

    if site_elem_i == "Ir":
        if len(nn_info_i) != 6:
            print("IOPSJDIFIDSJFID")

    if site_elem_i == "O":
        tmp = nn_info_i

    neighbor_list = []
    for neighbor_j in nn_info_i:
        neigh_elem_j = neighbor_j["site"].species_string
        neighbor_list.append(neigh_elem_j)

    neighbor_count_dict = dict()
    for i in neighbor_list:
        neighbor_count_dict[i] = neighbor_count_dict.get(i, 0) + 1
    
    data_dict_i["neighbor_count"] = neighbor_count_dict
    
    data_master.append(data_dict_i)

df_tmp = pd.DataFrame(data_master)
# -

# # Finding out if octahedra are corner sharing edge sharing

# +
nn_info_i

def number_of_neighbors(nn_info):
    num_neighbors = len(nn_info)
    return(num_neighbors)



def atom_is_in_central_octahedra(nn_info, verbose=True):
    is_octahedra = False

    correct_number_of_ligands = False
    correct_ligand_type = False
    

    num_nn = number_of_neighbors(nn_info)

    if num_nn == 6:
        correct_number_of_ligands = True
        if verbose:
            print("6 nearest neighbors!")

    nn_elems_unique = list(set([i["site"].species_string for i in nn_info_i]))
    if (len(nn_elems_unique) == 1) & (nn_elems_unique[0] == "O"):
        correct_ligand_type = True
        if verbose:
            print("Only 1 type of NN and it's oxygen")
    
    if correct_ligand_type and correct_number_of_ligands:
        is_octahedra = True

    return(is_octahedra)
# -

atom_is_in_central_octahedra(nn_info_i, verbose=False)

# +
ir_indices_list = []
for i_cnt, site_i in enumerate(structure.sites):
    if site_i.species_string == "Ir":
        ir_indices_list.append(i_cnt)


for index_i in ir_indices_list:
    nn_info_i = CrysNN.get_nn_info(structure, index_i)

    neighbor_list_indics_i = [i["site_index"] for i in nn_info_i]

    for index_j in ir_indices_list:
        
        if index_i == index_j:
            continue

        print(index_i, index_j)

        nn_info_j = CrysNN.get_nn_info(structure, index_j)
        neighbor_list_indics_j = [i["site_index"] for i in nn_info_j]


        num_shared_neigh = len(list(
            set(neighbor_list_indics_j) & set(neighbor_list_indics_i)))

# + active=""
#
#
#
#

# + jupyter={}
# neigh_elem_j
# neighbor_list

# neighbor_count_dict


# + jupyter={}
# dir(tmp[0]["site"])

# # tmp[0]["site"].specie
# # tmp[0]["site"].species

# tmp[0]["site"].species_string


# # species
# # species_and_occu
# # species_string
