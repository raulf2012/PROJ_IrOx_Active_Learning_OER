# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_IrOx_Active_Learning_OER]
#     language: python
#     name: conda-env-PROJ_IrOx_Active_Learning_OER-py
# ---

# # Import Modules

# +
import os
import sys
import pickle

import time
t0 = time.time()
sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import bulk_dft_data_path

# #############################################################################
from ase_modules.ase_methods import view_in_vesta

import pandas as pd

from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.analysis.local_env import NearNeighbors, VoronoiNN, site_is_of_motif_type

from methods import site_is_of_motif_type

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
print("df_bulk_dft.shape:", df_bulk_dft.shape)

df_bulk_dft = df_bulk_dft[
    (df_bulk_dft["source"] != "oqmd") & \
#     (df_bulk_dft["source"] != "raul_oer") & \
    (df_bulk_dft["source"] != "chris") & \
    [True for i in range(len(df_bulk_dft))]
    ]

print("df_bulk_dft.shape:", df_bulk_dft.shape)
# df_bulk_dft = df_bulk_dft.drop(all_ids_to_elim)
print("df_bulk_dft.shape:", df_bulk_dft.shape)

# +
df_bulk_dft = df_bulk_dft.sort_values("energy_pa")

# df_bulk_dft = df_bulk_dft.loc[["zk9q9yn3b2"]]
# df_bulk_dft = df_bulk_dft.loc[["c48lx363be"]]
# df_bulk_dft = df_bulk_dft.loc[["9ochme8294"]]

# df_bulk_dft = df_bulk_dft.iloc[0:50]
# -

# # METHOD | get_motiff_distribution

# +
thresh_dict = {
    "qtet": 0.5,
    "qoct": 0.5,
    # "qoct": 0.4,
    "qbcc": 0.5,
    "q6": 0.4,
    "qtribipyr": 0.8,
    "qsqpyr": 0.8,
    }

for key, val in thresh_dict.items():
    tmp = 42
    
    thresh_dict[key] = val / 1.5
# -

# # TEMP | Using NN from CrystalNN and my custom motiff method

path_i = os.path.join(
    "out_data",
    "coord_data_dict.pickle")
with open(path_i, "rb") as fle:
    coord_data_dict = pickle.load(fle)

# +
coord_data_dict
df_bulk_dft.head()

row_i = df_bulk_dft.loc["cubqbpzd7k"]
id_unique = row_i.name

coord_data_i = coord_data_dict.get(id_unique, None)

metal_index_j = 70
coord_data_j = coord_data_i[coord_data_i["structure_index"] == metal_index_j]

nn_info_list = coord_data_j["nn_info"].iloc[0]

nn_info_list = [i["site"] for i in nn_info_list][0]


# -

def get_motiff_distribution(atoms, nn_info):
    """
    """
    atoms_i = atoms

    struct_i = AseAtomsAdaptor.get_structure(atoms_i)

    metal_species_index_list = []
    for j_cnt, site_j in enumerate(struct_i):
        if site_j.species_string == "Ir":
            metal_species_index_list.append(j_cnt)

    coord_data_i = nn_info

    motiff_list_i = []
    for metal_index_j in metal_species_index_list:

        coord_data_j = coord_data_i[coord_data_i["structure_index"] == metal_index_j]
        nn_info_list = coord_data_j["nn_info"].iloc[0]
        nn_info_list = [i["site"] for i in nn_info_list]
        print(nn_info_list)

        motiff_j = site_is_of_motif_type(
            struct_i,
            metal_index_j,
            neighbors_list=nn_info_list,
            # "min_dist", "voronoi", "min_OKeeffe", "min_VIRE"
            approach="min_dist",
            delta=0.1,
            # delta=0.3,
            cutoff=10.0,
            thresh=thresh_dict,
            )

        motiff_list_i.append(motiff_j)


    motiff_count_dict = {}
    list_i = motiff_list_i
    for items in list(set(list_i)):
        motiff_count_dict[items] = list_i.count(items)

    num_motiffs = len(motiff_list_i)

    motiff_frac_dict = {}
    for key, val in motiff_count_dict.items():
        motiff_frac_i = val / num_motiffs
        motiff_frac_dict[key] = motiff_frac_i
        
    return(motiff_frac_dict)


# + jupyter={"outputs_hidden": true}
def method(row_i, get_motiff_distribution):
    # row_i = df_bulk_dft.loc["cubqbpzd7k"]
    atoms_i = row_i["atoms"]
    id_unique = row_i.name

    nn_info_i = coord_data_dict.get(id_unique, None)

    motiff_distr_i = get_motiff_distribution(atoms_i, nn_info=nn_info_i)
    return(motiff_distr_i)

coord_motiffs_series = df_bulk_dft.apply(
    method,
    axis=1, args=(get_motiff_distribution, ))

df_coord = pd.DataFrame(
    coord_motiffs_series,
    columns=["coord_motiff_distr"])


# +
def method(row_i):
    index_i = row_i.name

    coord_motiff_distr_i = row_i["coord_motiff_distr"]

    df_i = pd.DataFrame(coord_motiff_distr_i, index=["col"]).T
    df_i = df_i.sort_values("col", ascending=False)
    tmp = df_i.iloc[0]["col"]
    df_dominant_motiffs_i = df_i[df_i["col"] == tmp]

    out_string = "_".join(df_dominant_motiffs_i.index.tolist())

    return(out_string)

# #############################################################################
df_coord["major_motiffs"] = df_coord.apply(
    method,
    axis=1)

df_coord.head()
# -

df_coord

# +
unrec_list = []
for i_cnt, row_i in df_coord.iterrows():
    if "unrecognized" in row_i["major_motiffs"]:
        unrec_list.append(row_i)

len(unrec_list)

# + active=""
# 264
# 255 | Voronoi
# -

# # Save Data

# +
# # Pickling data ######################################################
# import os; import pickle
# directory = "out_data"
# if not os.path.exists(directory): os.makedirs(directory)
# with open(os.path.join(directory, "df_coord_motiff.pickle"), "wb") as fle:
#     pickle.dump(df_coord, fle)
# # #####################################################################

# +
# df_coord["coord_motiff_distr"].tolist()

# df_coord.head()
# -

print("Time to execute notebook: ", time.time() - t0, "(s)")
print("os.getcwd():", os.getcwd())

# + active=""
#
#
#
#

# +
coord_motiff_keys = []
for i in df_coord["coord_motiff_distr"].tolist():
    for key, val in i.items():
#         print(key)
        coord_motiff_keys.append(key)

set(coord_motiff_keys)
# -

df_coord["major_motiffs"].unique().tolist()

# +
# ['octahedral',
#  'bcc',
#  'unrecognized',
#  'trigonal bipyramidal',
#  'tetrahedral',
#  'cp',
#  'square pyramidal_bcc'
# ]

# + jupyter={}
# # df_bulk_dft.sort_values("energy_pa")

# row_i = df_bulk_dft.loc["n36axdbw65"]

# from ase.visualize import view
# from ase_modules.ase_methods import view_in_vesta

# # view_in_vesta(atoms_i)
# # view(atoms_i)

# atoms_i.write("~/temp.cif")

# VNN.site_is_of_motif_type(
#     struct_i,
#     0,
#     approach='min_dist',
#     delta=0.1,
#     cutoff=10.0,
#     thresh=None,
#     )

# structure_from_cif = Structure.from_file()
# y.get_bonded_structure(structure_from_cif)
# y.get_local_order_parameters(structure_from_cif, 0)

# site_is_of_motif_type(struct, n, approach="min_dist", delta=0.1, \          
# 1697                           cutoff=10.0, thresh=None):                            
# 1698     """                        

# import pymatgen

# print(pymatgen.__version__)

# print(pymatgen)

# metal_species_index_list = []
# for j_cnt, site_j in enumerate(struct_i):
#     if site_j.species_string == "Ir":
#         metal_species_index_list.append(
#             j_cnt)

# metal_species_index_list

# row_i = df_bulk_dft.loc["n36axdbw65"]

# # row_i = df_bulk_dft.iloc[0]
# atoms = row_i["atoms"]

# motiff_distr_i = get_motiff_distribution(atoms)
