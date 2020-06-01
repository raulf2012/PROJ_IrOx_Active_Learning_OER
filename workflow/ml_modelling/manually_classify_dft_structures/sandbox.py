# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_IrOx_Active_Learning_OER]
#     language: python
#     name: conda-env-PROJ_IrOx_Active_Learning_OER-py
# ---

# + jupyter={}
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

# + jupyter={}
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling",
    "ccf_similarity_analysis/out_data",
    "all_ids_to_elim_1.pickle")
with open(path_i, "rb") as fle:
    all_ids_to_elim = pickle.load(fle)

with open(bulk_dft_data_path, "rb") as fle:
    df_bulk_dft = pickle.load(fle)

# + jupyter={}
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
row_i = df_bulk_dft.iloc[0]

atoms_i = row_i["atoms"]
# -

row_i

# +
# atoms_i = atoms

struct_i = AseAtomsAdaptor.get_structure(atoms_i)

metal_species_index_list = []
for j_cnt, site_j in enumerate(struct_i):
    if site_j.species_string == "Ir":
        metal_species_index_list.append(j_cnt)

motiff_list_i = []
for metal_index_j in metal_species_index_list:
    tmp = 41

# +
# struct_i,
# metal_index_j,

# # "min_dist", "voronoi", "min_OKeeffe", "min_VIRE"
# approach="min_dist",
# delta=0.1,
# # delta=0.3,
# cutoff=10.0,
# thresh=thresh_dict,
# -

struct = struct_i
n = metal_index_j
approach = "min_dist"
delta = 0.1
cutoff = 10.0
thresh = None

from pymatgen.analysis.local_env import LocalStructOrderParams, get_neighbors_of_site_with_index

# +
if thresh is None:
    thresh = {
        "qtet": 0.5, "qoct": 0.5, "qbcc": 0.5, "q6": 0.4,
        "qtribipyr": 0.8, "qsqpyr": 0.8}

ops = LocalStructOrderParams([
    "cn", "tet", "oct", "bcc", "q6", "sq_pyr", "tri_bipyr"])

neighs_cent = get_neighbors_of_site_with_index(
    struct, n, approach=approach, delta=delta, cutoff=cutoff)

neighs_cent.append(struct.sites[n])
opvals = ops.get_order_parameters(
    neighs_cent, len(neighs_cent) - 1, indices_neighs=[
        i for i in range(len(neighs_cent) - 1)])
cn = int(opvals[0] + 0.5)
motif_type = "unrecognized"
nmotif = 0

if cn == 4 and opvals[1] > thresh["qtet"]:
    motif_type = "tetrahedral"
    nmotif += 1
if cn == 5 and opvals[5] > thresh["qsqpyr"]:
    motif_type = "square pyramidal"
    nmotif += 1
if cn == 5 and opvals[6] > thresh["qtribipyr"]:
    motif_type = "trigonal bipyramidal"
    nmotif += 1
if cn == 6 and opvals[2] > thresh["qoct"]:
    motif_type = "octahedral"
    nmotif += 1
if cn == 8 and (opvals[3] > thresh["qbcc"] and opvals[1] < thresh["qtet"]):
    motif_type = "bcc"
    nmotif += 1
if cn == 12 and (opvals[4] > thresh["q6"] and opvals[1] < thresh["q6"] and
                 opvals[2] < thresh["q6"] and opvals[3] < thresh["q6"]):
    motif_type = "cp"
    nmotif += 1

if nmotif > 1:
    motif_type = "multiple assignments"


# +
opvals


# -

motif_type
