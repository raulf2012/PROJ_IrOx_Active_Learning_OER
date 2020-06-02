# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_irox] *
#     language: python
#     name: conda-env-PROJ_irox-py
# ---

# # Compute CCF for Each Structure
# ---
#
#
# CCF is then used to compute simalirty for all pair-wise comparisons of all structures
#
# NOTE: This notebook can be set to read df_ccf from PROJ_DATA directory instead of actually running this script fully/correctly. This is becauses it takes almost an hour to run this notebook. To activate this mode set the `read_from_PROJ_DATA` to `True`

read_from_PROJ_DATA = False
read_from_PROJ_DATA = True

# # Import Modules

# +
import os
print(os.getcwd())
import sys

import pickle

import numpy as np

import pandas as pd

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import (
    static_irox_structures_path,
    bulk_dft_data_path,
    unique_ids_path,
    )

# from StructurePrototypeAnalysisPackage.ccf import struc2ccf
from spap.ccf import struc2ccf
# -

# # Script Inputs

# +
r_cut_off = 10.
r_vector = np.arange(1, 10, 0.02)

mean_density = 0.08407356
# -

directory = "out_data"
if not os.path.exists(directory):
    os.makedirs(directory)

# # Read Data

# +
with open(static_irox_structures_path, "rb") as fle:
    df_static_irox = pickle.load(fle)

with open(bulk_dft_data_path, "rb") as fle:
    df_bulk_dft = pickle.load(fle)

df_ids = pd.read_csv(unique_ids_path)

# try:
if read_from_PROJ_DATA:
    path_i = os.path.join(
        os.environ["PROJ_DATA"],
        "04_IrOx_surfaces_OER/PROJECT_COMPUTED_OUT_DATA/PROJ_IrOx_Active_Learning_OER",
        "workflow/ml_modelling/ccf_similarity_analysis/compute_ccf_and_dij_matrix/",
        "out_data/df_ccf.pickle"
        )
    with open(path_i, "rb") as fle:
        df_ccf_prev = pickle.load(fle)

else:
    with open("out_data/df_ccf.pickle", "rb") as fle:
        df_ccf_prev = pickle.load(fle)
# except:
#     # df_ccf_prev
#     df_ccf_prev = pd.DataFrame()

# +
# df_ccf_prev.shape

# +
# assert False
# -

# # Filtering df_bulk_dft

# +
sources_to_keep = [
    'raul',
    'raul_oer',
    # 'chris', 'oqmd',
    ]

df_bulk_dft = df_bulk_dft[df_bulk_dft["source"].isin(sources_to_keep)]
# -

# # Combining different datasets to process uniformly

df_static_irox = df_static_irox.set_index("static_id", drop=False)

# +
df_m = pd.concat([
    df_bulk_dft,
    df_static_irox,
    ], sort=False,)

print("df_tmp.index.shape:", df_m.index.shape)
print("df_tmp.index.unique().shape:", df_m.index.unique().shape)


# -

# # Creating scaled atoms with equal atomic density

# +
def get_atomic_density(row_i):
    atoms = row_i["atoms"]
    volume = atoms.get_volume()
    num_atoms = atoms.get_number_of_atoms()
    density = num_atoms / volume
    return(density)

df_bulk_dft_scaled = pd.DataFrame()
# df_bulk_dft_scaled["density_init"] = df_bulk_dft.apply(
df_bulk_dft_scaled["density_init"] = df_m.apply(
    get_atomic_density,
    axis=1)

# mean_density = df_bulk_dft_scaled["density_init"].mean()
# print("mean_density:", mean_density, "atoms/A3")
# assert False

# #############################################################################
# #############################################################################


data_list = []
for i_cnt, (name_i, row_i) in enumerate(df_m.iterrows()):
    atoms_i = row_i["atoms"]
    row_scaled = df_bulk_dft_scaled.loc[name_i]
    dens_init = row_scaled["density_init"]
    scale_fact = (mean_density / dens_init) ** (1 / 3)
    new_cell = atoms_i.cell / scale_fact
    atoms_i.set_cell(new_cell, scale_atoms=True)
    dens_final = atoms_i.get_number_of_atoms() / atoms_i.get_volume()

    out_dict = {
        "atoms_scaled": atoms_i,
        "index": name_i,
        "cell_scale_factor": scale_fact,
        "density_final": dens_final}
    data_list.append(out_dict)

df_scaled_atoms = pd.DataFrame(data_list).set_index("index")
df_bulk_dft_scaled = pd.concat([df_scaled_atoms, df_bulk_dft_scaled], axis=1)

num_unique_ids = df_bulk_dft_scaled.index.unique().shape[0]
num_ids = df_bulk_dft_scaled.index.shape[0]
assert num_unique_ids == num_ids, "JISFIDSIFJ"

# +
# df_scaled_atoms.head()

# df_bulk_dft_scaled.head()

# + active=""
#
#
#
#
#
#
#
#
#
# -

# # Calculate CCF for DFT Calculated IrO2 and IrO3 Systems

# +
# df_ccf.index

"bsv4nex29l" in df_bulk_dft.index

# df_bulk_dft

# +
indices_to_process = [i for i in df_bulk_dft_scaled.index if i not in df_ccf_prev.index]
print("len(indices_to_process):", len(indices_to_process))
index_before_splitting = df_bulk_dft_scaled.index

df_bulk_dft_scaled_not_processed = df_bulk_dft_scaled.loc[indices_to_process]
df_not_proc = df_bulk_dft_scaled_not_processed

if len(indices_to_process) == 0:
    print("No systems to process, exiting")

    with open("out_data/df_ccf.pickle", "wb") as fle:
        pickle.dump(df_ccf_prev, fle)

    print(20 * "# # ")
    print("All done!")
    assert False


# -
def method(row_i, argument_0, atoms_key="atoms"):
    """
    """
    atoms_i = row_i[atoms_key]
    print(20 * "*")
    ccf_i = struc2ccf(atoms_i, r_cut_off, r_vector)
    return(ccf_i)


# +
# assert False

# +
# TEMP
# df_not_proc = df_not_proc.iloc[0:10]

# +
df_i = df_not_proc
df_i["ccf"] = df_i.apply(
    method,
    axis=1,
    args=("TEMP", ),
    atoms_key="atoms_scaled"
    )
df_bulk_dft_scaled = df_i

df_ccf = df_bulk_dft_scaled["ccf"]
df_ccf = pd.DataFrame(df_ccf)

with open("out_data/df_ccf.pickle", "wb") as fle:
    pickle.dump(df_ccf, fle)

# +
df_ccf_new = pd.concat([
    df_ccf,
    df_ccf_prev,
    ], sort=False,
    )
df_ccf_new.shape

with open("out_data/df_ccf.pickle", "wb") as fle:
    pickle.dump(df_ccf_new, fle)

# +
with open("out_data/df_ccf.pickle", "rb") as fle:
    df_ccf_tmp = pickle.load(fle)

# with open("out_data/df_ccf.pickle", "rb") as fle:
#     df_ccf = pickle.load(fle)
# -

print("df_bulk_dft_scaled.shape:", df_bulk_dft_scaled.shape)
print("df_ccf_prev.shape:", df_ccf_prev.shape)
print("df_ccf_new.shape:", df_ccf_new.shape)

print(20 * "# # ")
print("All done!")
assert False

# + active=""
#
#
#
#
#

# + jupyter={}
# df_tmp.reindex(index_before_splitting)

# len(index_before_splitting)

# df_ccf_tmp.shape

# print(df_ccf_prev.index.unique().shape)


# assert df_bulk_dft_scaled.index.unique().shape[0] == df_bulk_dft_scaled.index.shape[0]

# indices_to_process = [i for i in df_bulk_dft_scaled.index if i not in df_ccf_prev.index]

# df_bulk_dft_scaled_not_processed = df_bulk_dft_scaled.loc[indices_to_process]
# df_not_proc = df_bulk_dft_scaled_not_processed

# df_i = df_bulk_dft_scaled
# df_i["ccf"] = df_i.apply(
#     method,
#     axis=1,
#     args=("TEMP", ),
#     atoms_key="atoms_scaled"
#     )
# df_bulk_dft_scaled = df_i

# df_ccf = df_bulk_dft_scaled["ccf"]
# df_ccf = pd.DataFrame(df_ccf)

# with open("out_data/df_ccf_test.pickle", "wb") as fle:
#     pickle.dump(df_ccf, fle)

# # TEMP | Reducing datasize for testing

# # TEMP
# # df_static_irox = df_static_irox.iloc[0:10]
# # df_bulk_dft = df_bulk_dft.iloc[0:10]


# ids_to_keep = [
# 'zy9dzknhnj',
# '6r716sxr9t',
# 'zlbp9jm19g',
# 'm16tvk8h7g',
# 'z39g648rnl',
# 'x18q7u8s9q',
# 'mp659s8i8q',
# 'ck638t75z3',
# '8u94muxqzw',
# 'mkbrzh8kv5',
# 'bwmrbuztvs',
# '82mwbwbgbl',
# '9knk63mjcw',
# '8495vjvsc3',
# 'xonhb1cgnq',
# '8aza84xlvs',
# '9wzqvtmw6a',
# 'n4xq8uv1cu',
# 'c48lx363be',
# 'bovtbhmf6p',
# 'vt658pbemi',
# '8hzf7ecw84',
# 'cf9lbhvjmq',
# 'bdzo7ixszl',
# 'brbizonjmy',
# 'zhci7ans9k',
# 'c3mp6jmgzq',
# '6fzy7kcrxy',
# 'cln3xecgbt',
# 'ckbg7d9u6u',
# 'mpbd6sz1z5',
# 'cant6qca6o',
# '6avov5cy64',
# '946gc3z2n3',
# '95zh7p8kno',
# 'cdm48gxynp',
# 'zunqxh6g9j',
# '9y85x1cp6e',
# 'xgmkzuxf8g',
# 'zuzwxhvuxe',
# 'cfcivdxrc2',
# 'm29j648g6i',
# '9g7f8evqzt',
# 'vhck6u84bj',
# 'mlva81bg7p',
# 'vj8l8uba9y',
# 'nt6wv1n27r',
# 'm4cacfvjvh',
# '9uxemw7rv2',
# 'bt8rz3babw',
       
# 'IrO3_rutile-like',
# 'IrO3',
# 'IrO3_battery',
# 'IrO2',
#       ]

# # df_bulk_dft[0:50].index.tolist()
# df_bulk_dft[df_bulk_dft["source"] == "raul_oer"].index.tolist()
# # df_bulk_dft = df_bulk_dft.loc[ids_to_keep]

# df_bulk_dft[df_bulk_dft.stoich == "AB3"].sort_values("dH").iloc[0:120]
# # ~0.3 eV/atom

# df_bulk_dft[df_bulk_dft.stoich == "AB2"].sort_values("dH").iloc
# # ~0.5 eV/atom

# ab2_min = df_bulk_dft[df_bulk_dft.stoich == "AB2"].dH.min()
# ab3_min = df_bulk_dft[df_bulk_dft.stoich == "AB3"].dH.min()

# from proj_data_irox import calc_dH

# print(calc_dH(-6.542, stoich="AB2"))

# print(calc_dH(-6.163, stoich="AB3"))

# ab3_min - -0.3438547784081729

# ab2_min - -0.33285956787756277

# -0.3068474600000002
# -0.5070619763541669

# with open(bulk_dft_data_path, "rb") as fle:
#     df_bulk_dft = pickle.load(fle)

# df_tmp = df_bulk_dft[df_bulk_dft.source == "chris"]

# -0.838172  # From chris | 9yz2mt8hbh
# -0.839922  # Mine | cubqbpzd7k

# -0.838286  # Mine | 9yz2mt8hbh

# df_tmp[df_tmp.stoich == "AB2"].sort_values("dH")
