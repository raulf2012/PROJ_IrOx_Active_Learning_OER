# ---
# jupyter:
#   jupytext:
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

# # Analysing Similarity Matrix for IrOx Systems Post-DFT
# ---
#
# Systems that have the same d but different energies
#
# list_0 = [
#  '8481z1n1na',
#  'zr9ic2zaz5',
#  '8h9snabqca',
#  '7f8pm5mhnu',
#  'cgx3mkzhmd',
#  'vwxfn3blxi',
#  '9obw8dbrvy',
#  'bpvynr7p9w',
#  '8gnovr727t',
#
#
#  '9pb4c1927h',
#  '8i63m2b5ve',
#
#
#  'vlxp9abd6h',
#  'z2nh817ene',
#  'xu6ivyvfvf',
#  ]

# # Import Modules

# +
import os
import sys

import pickle
import pandas as pd

# #############################################################################
import plotly.graph_objs as go

# #############################################################################
sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import (
    static_irox_structures_path,
    bulk_dft_data_path,
    df_dij_path)

from methods import plot_dij_matrix_heatmap
from plotting.my_plotly import my_plotly_plot
# -

# # Script Inputs

# +
# d_thresh = 0.075

# d_thresh = 0.01
d_thresh = 0.02
# d_thresh = 0.03
# d_thresh = 0.04
# d_thresh = 0.05
# d_thresh = 0.06
# d_thresh = 0.07
# d_thresh = 0.08
# d_thresh = 0.09
# d_thresh = 0.10
# d_thresh = 0.20
# d_thresh = 0.30
# d_thresh = 0.40
# d_thresh = 0.70


e_thresh = 0.01

create_plot = True
# -

# # Read Data

# +
# df_dij_path_tmp = df_dij_path[0:-18] + "df_d_ij_all_temp.pickle"
with open(df_dij_path, "rb") as fle:
# with open(df_dij_path_tmp, "rb") as fle:
    df_dij_dft = pickle.load(fle)
    print("df_dij_dft.shape:", df_dij_dft.shape)

with open(static_irox_structures_path, "rb") as fle:
    df_static_irox = pickle.load(fle)

with open(bulk_dft_data_path, "rb") as fle:
    df_bulk_dft = pickle.load(fle)

path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling",
    "ccf_similarity_analysis/out_data",
    "all_ids_to_elim_1.pickle")
with open(path_i, "rb") as fle:
    ids_to_drop_prev = pickle.load(fle)

ids_to_drop_prev = ids_to_drop_prev["AB2"] + ids_to_drop_prev["AB3"]

# sys.path.insert(0, "../04_final_ml_plots")

# +
# # df_dij_dft.loc["8p8evt9pcg", "9lmkmh8s8r"]


# df_dij_dft.loc[

#     "64cg6j9any",
#     "b46enqnq8e",
#     "9yz2mt8hbh",

# #     "6avov5cy64"
    
# #     "clc2b1mavs",
#     ]
# -

# # Dropping Static Structure from D_ij

# +
static_ids = df_static_irox["static_id"].tolist()
static_ids_in_dij = [i for i in static_ids if i in df_dij_dft.index]

df_dij_dft = df_dij_dft.drop(labels=static_ids_in_dij, axis=0)
df_dij_dft = df_dij_dft.drop(labels=static_ids_in_dij, axis=1)
# -

# # Filtering data to needed systems

# +
df_bulk_dft = df_bulk_dft[
    (df_bulk_dft["source"] != "chris") &
    (df_bulk_dft["source"] != "oqmd") &
    [True for i in range(len(df_bulk_dft))]
    ]

print("df_bulk_dft.shape:", "\n", df_bulk_dft.shape)
print("df_bulk_dft.index.unique().shape:", "\n",
    df_bulk_dft.index.unique().shape)
# -

# # Reorder index by Stoicheomtry first and then by energy

# +
ab2_indices = df_bulk_dft[df_bulk_dft["stoich"] == "AB2"].sort_values(
    "energy_pa").index.tolist()
ab3_indices = df_bulk_dft[df_bulk_dft["stoich"] == "AB3"].sort_values(
    "energy_pa").index.tolist()

ab2_indices_not_in_dij = [i for i in ab2_indices if i not in df_dij_dft.index]

new_ind_order = ab2_indices + ab3_indices
new_index_order_filtered = [i for i in new_ind_order if i in df_dij_dft.index]

df_dij_dft = df_dij_dft.reindex(new_index_order_filtered)
df_dij_dft = df_dij_dft[new_index_order_filtered]
# -

print("len(ab2_indices):", len(ab2_indices))
print("len(ab3_indices):", len(ab3_indices))
print("")
print("df_dij_dft.shape:", df_dij_dft.shape)

# # Reorder index to put OER bulk systems first

# +
oer_sys_ids = ['IrO3_rutile-like', 'IrO3', 'IrO3_battery', 'IrO2']

non_oer_ids = df_dij_dft.index.drop(oer_sys_ids)
new_index_order = oer_sys_ids + non_oer_ids.tolist()

df_dij_dft = df_dij_dft.reindex(new_index_order)
df_dij_dft = df_dij_dft[new_index_order]
# -

# # Drop ids that were identified to be redundant

# +
# df_dij_dft = df_dij_dft.drop(labels=ids_to_drop_prev, axis=0)
# df_dij_dft = df_dij_dft.drop(labels=ids_to_drop_prev, axis=1)
# -

df_dij_dft.loc["IrO3_rutile-like"][df_dij_dft.loc["IrO3_rutile-like"] < 0.01]

# # Create D_ij Matrix Plot

if create_plot:
    data = plot_dij_matrix_heatmap(
        df_dij_dft,
        d_thresh,
        e_thresh)

    layout = go.Layout(width=1100, height=1100)
    fig = go.Figure(data=data, layout=layout)

    fig = my_plotly_plot(
        figure=fig,
        plot_name='irox_dij_heatmap',
        # write_pdf_svg=True,
        write_html=True,
        write_png=True,
        write_pdf=False,
        write_svg=False,
        )

# +
# fig
# -

# # Analyzing systems that are duplicates

df_dij_ab2 = df_dij_dft.loc[ab2_indices, ab2_indices]
df_dij_ab3 = df_dij_dft.loc[ab3_indices, ab3_indices]


def ids_to_elim(df_dij):
    """
    """
    index_to_eliminate = []
    for i_cnt, (name_i, row_i) in enumerate(df_dij.iterrows()):
        cols_below_thresh = row_i[row_i < d_thresh]
        if cols_below_thresh.shape[0] > 1:
            df_i = df_bulk_dft.loc[cols_below_thresh.index]
            index_to_eliminate += df_i.iloc[1:].index.tolist()

    index_to_eliminate = list(set(index_to_eliminate))

    return(index_to_eliminate)


# +
ids_to_elim_ab2 = ids_to_elim(df_dij_ab2)
ids_to_elim_ab3 = ids_to_elim(df_dij_ab3)

all_ids_to_elim = {
    "AB2": ids_to_elim_ab2,
    "AB3": ids_to_elim_ab3,
    }


print("len(ids_to_elim_ab2):", len(ids_to_elim_ab2))
print("len(ids_to_elim_ab3):", len(ids_to_elim_ab3))

# + active=""
# len(ids_to_elim_ab2): 43
# len(ids_to_elim_ab3): 69
#
# len(ids_to_elim_ab2): 92
# len(ids_to_elim_ab3): 103
# -

# ## Saving ids of systmes that are duplicates

# +
# Pickling data ######################################################
import os; import pickle
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)

# TODO | Don't create this one anymore
with open(os.path.join(directory, "all_ids_to_elim_1.pickle"), "wb") as fle:
    pickle.dump(all_ids_to_elim, fle)

with open(os.path.join(directory, "all_ids_to_elim.pickle"), "wb") as fle:
    pickle.dump(all_ids_to_elim, fle)
# #####################################################################
# -

df_bulk_dft[df_bulk_dft.stoich == "AB2"].sort_values("dH")

# +
# df_dij_ab3.loc["8p8evt9pcg", "zimixdvdxd"]

df_dij_ab3.loc["xw9y6rbkxr", "zimixdvdxd"]
# -

assert False

# # TEST | TEST | TEST

# +
# df_dij_dft.loc["7h7yns937p"]
df_dij_dft.shape

"7h7yns937p" in df_dij_dft.index

# +

ids_to_elim_ab3 = all_ids_to_elim["AB3"]

print(len(ab3_indices))

unique_ids_ab3 = [i for i in ab3_indices if i not in ids_to_elim_ab3]

data_dict_list = []
for id_i in unique_ids_ab3:
    if id_i in df_dij_dft.index:
        num_duplicates = len(df_dij_dft.loc[id_i][df_dij_dft.loc[id_i] < d_thresh]) - 1,
        dict_i = {
            "id_unique": id_i,
            "num_duplicates": num_duplicates[0],
            }
        data_dict_list.append(dict_i)
    else:
        pass


df_tmp = pd.DataFrame(data_dict_list)

df_tmp.sort_values("num_duplicates", ascending=False)
# data_dict_list
# -



# +
# TEMP
# df_dij_dft = df_dij_dft.loc[all_ids_to_elim, all_ids_to_elim]


ids_dict_master = {}
for i_cnt, (name_i, row_i) in enumerate(df_dij_dft.iterrows()):
    # tmp = row_i[row_i < d_thresh]
    tmp = row_i[row_i < d_thresh].drop(name_i)
    df_i = df_bulk_dft.loc[tmp.index]

    # ids_dict_list_i = {i_cnt: df_i.index.sort_values().tolist()}
    # ids_dict_lists.append(ids_dict_list_i)
    if len(df_i) > 0:
        equiv_ids_list = df_i.index.sort_values().tolist()
        id_joined_str = "_".join(equiv_ids_list)

        # ids_dict_master[i_cnt] = df_i.index.sort_values().tolist()
        ids_dict_master[name_i] = {
            "id_joined_str": id_joined_str,
            "equiv_ids_list": equiv_ids_list,
            }

# #############################################################################
# df_i.index.sort_values().tolist()

df_test = pd.DataFrame(ids_dict_master,
#     index=ids_dict_master.keys()
#     index=["id_str_joined"],
    ).T


# df_test["id_str_joined"].unique().shape

df_test

# + active=""
#
#
#
#
# -

df_dij_dft.loc["8p8evt9pcg", "zimixdvdxd"]

# +
# ids_dict_master
# tmp_list = []
# for key_i, val_i in ids_dict_master.items():
#     for key_j, val_j in ids_dict_master.items():

#         if key_i == key_j:
#             continue

#         print(key_i, key_j)

#         if val_j == val_i:
#         else:
#             tmp_list.append(key_i)

# ids_dict_master

# if val_j == val_i:

# len(all_ids_to_elim)


# print(len([i for i in all_ids_to_elim if i in ab2_indices]))
# print(len([i for i in all_ids_to_elim if i in ab3_indices]))

# # np.fill_diagonal(df_dij_dft.values, np.nan)
# # e_thresh = 0.01
# use_energy_simil = False

# trouble_ids_list = []

# unique_id_list = []
# all_ids_to_elim = []
# for i_cnt, (name_i, row_i) in enumerate(df_dij_dft.iterrows()):
#     tmp = row_i[row_i < d_thresh]

#     # if len(tmp) > 1:
#     #     break

#     if tmp.shape[0] == 1:
#         mess = "No other structures close to this one"
#         # print(mess)
#         unique_id_list.append(tmp.index[0])
#     else:
#         df_i = df_bulk_dft.loc[tmp.index]

#         # if "8k7expx2bp" in df_i.index.tolist():
#         # if "6s648e8s6p" in df_i.index.tolist():
#         # if "b5cgvsb16w" in df_i.index.tolist():
#         #     display(df_i)

#         e_range = abs(df_i["energy_pa"].min() - df_i["energy_pa"].max())
#         e_thresh_u = df_i.loc[name_i]["energy_pa"] + e_thresh
#         e_thresh_l = df_i.loc[name_i]["energy_pa"] - e_thresh

#         # Using enery similarity criteria
#         if use_energy_simil:
#             df_j = df_i[
#                 (df_i["energy_pa"] < e_thresh_u) &
#                 (df_i["energy_pa"] > e_thresh_l)]
#             index_to_keep_i = df_j.sort_values("energy_pa").iloc[0].name
#             index_to_eliminate = df_j.iloc[1:].index.tolist()
#         else:
#             index_to_keep_i = df_i.sort_values("energy_pa").iloc[0].name
#             index_to_eliminate = df_i.iloc[1:].index.tolist()


#         all_ids_to_elim += index_to_eliminate

#         if e_range > e_thresh:
#             # display(df_i)
#             df_i_tmp = df_i

#             ids_tmp = df_i_tmp.index.tolist()
#             trouble_ids_list += ids_tmp

#             # print("Energies span greater range than 'e_thresh'")
#             # print(e_range)
#             # print("")

# all_ids_to_elim = list(set(all_ids_to_elim))

# trouble_ids_list = list(set(trouble_ids_list))

# #############################################################################
# Drop ab2 stoicheomtry

# df_dij_dft = df_dij_dft.drop(labels=ab2_indices, axis=0)
# df_dij_dft = df_dij_dft.drop(labels=ab2_indices, axis=1)

# dft_indices = ab2_indices + ab3_indices
# non_dft_indices = [i for i in df_dij_dft.index if i not in dft_indices]
