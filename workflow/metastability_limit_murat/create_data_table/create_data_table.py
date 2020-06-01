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

# +
import os
print(os.getcwd())
import sys

import pandas as pd

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import calc_dH

# +
iro2_dft_energies = [
    -662.273,
    -665.361,
    -659.370,
    -663.256,
    -665.901,
    -667.284,
    ]


iro3_dft_energies = [
    -612.420,
    -611.473,
    -615.891,
    -609.275,
    -608.229,
    -616.290,
    ]


# +
# stoich_i = "AB2"

# if stoich_i == "AB2":
#     raw_energies = iro2_dft_energies
# elif stoich_i == "AB3":
#     raw_energies = iro3_dft_energies
# -

def process_df(raw_energies, stoich_i=None):

    stoich_conv_dict = {
        "AB2": "$IrO_{2}$",
        "AB3": "$IrO_{3}$",
        }
    
    col_conv_dict = {
        "raw_dft": "$E_{DFT}$",
        "ev_atom": "$E_{DFT}$",
        "dH_atom": "$\Delta H_{f}$",
        "e_above_hull": "$\Delta E_{hull}$",
        }
    unit_conv_dict = {
        "raw_dft": "(eV)",
        "ev_atom": "(eV/atom)",
        "dH_atom": "(eV/atom)",
        "e_above_hull": "(eV/atom)",
        }

    num_atoms_dict = {
        "AB2": 102,
        "AB3": 100,
        }

    df = pd.DataFrame()

    df["raw_dft"] = raw_energies
    df["ev_atom"] = df.raw_dft / num_atoms_dict[stoich_i]
    df["dH_atom"] = [calc_dH(i, stoich=stoich_i) for i in df.ev_atom.tolist()]


    hull_e_per_atom_dict = {
        "AB2": -7.04751560624999,
        "AB3": -6.46984746,
        }

    df["e_above_hull"] = df.ev_atom - hull_e_per_atom_dict[stoich_i]

    # Sorting by energy
    df = df.sort_values("raw_dft")

    # Dropping raw DFT column
    # df = df.drop(columns=[df.columns[0]])
    df = df.drop(columns=["raw_dft"])

    [stoich_conv_dict.get(stoich_i, "TEMP") for i in df.columns.values],
    
    stoich_col_list = []
    for i_cnt, col_i in enumerate(df.columns.values):
        if i_cnt == 0:
            stoich_col_list.append(
                stoich_conv_dict.get(stoich_i, "TEMP")
                )
        else:

            import random
            import numpy as np

            rand_float = np.round(random.random(), decimals=3)

            stoich_col_list.append(
                "\phantom{" + str(i_cnt) + str(rand_float) + "}"
                )

    # print("stoich_col_list:", stoich_col_list)
        

    tuples = list(zip(*[
        stoich_col_list,
        [col_conv_dict.get(i, "TEMP") for i in df.columns.values],
        [unit_conv_dict.get(i, "TEMP") for i in df.columns.values],
        # ["a", "b", "c", ],
        ]))


    index = pd.MultiIndex.from_tuples(tuples, names=["Stoich.", "", "Units"])
    df.columns = index

    return(df)

# +
df_ab2 = process_df(iro2_dft_energies, stoich_i="AB2")
df_ab3 = process_df(iro3_dft_energies, stoich_i="AB3")


df = pd.concat([df_ab2, df_ab3], axis=1)
# -

# # Round values

# +
df = df.round(decimals=3)



# +
last_row_data = []
for i in df.iloc[-1].values:
    i = "\textbf{" + str(i) + "}"
    last_row_data.append(i)
    

# df.iloc[-1].values = last_row_data
# df.loc[5] = last_row_data
df.loc[5] = tuple(last_row_data)

# tmp = df.columns.values[0]
# df.loc[5][tmp] = "TEMP"

df

# + active=""
#
#
#
#
#
# -

# # Write Latex Table to File

# + jupyter={}
alignment_list = ["c" for i in range(len(df.columns))]

alignment_str = ""
for i in alignment_list:
    alignment_str += i

# #########################################################
shared_props = dict(
    # buf="oer_table.tex",
    columns=None,
    col_space=None,

    # #####################################################
    header=True,
    # #####################################################
    index=False,

    # #####################################################
    na_rep='-',
    formatters=None,

    # #####################################################
    # float_format="{:0.2f}",
    # float_format=["c" for i in range(len(df.columns))],
    # float_format=alignment_str,
    # float_format="%%.2f",

    sparsify=None,
    index_names=True,
    bold_rows=False,
    
    # #####################################################
    column_format=alignment_str,
    
    longtable=None,
    
    # #####################################################
    escape=False,
    encoding=None,
    decimal='.',
    multicolumn=None,
    multicolumn_format=None,
    multirow=None,
    )

# + jupyter={}
df.to_latex(
    buf="TEMP.tex",
    **shared_props)

path_i = os.path.join(
    os.environ["PROJ_irox_paper"],
    "04_data_tables/amorphous_limit_data",
    "amorphous_lim_table.tex")
df.to_latex(
    buf=path_i,
    **shared_props)

# +
# # df.to_latex?
# -

df.columns
# 
df.columns.value_counts()



# + active=""
#
#
#
#

# + jupyter={}

# # ab2_meta_lim = 
# calc_dH(i, stoich="AB2")

# + jupyter={}
# ab2_meta_lim = calc_dH(
#     raw_dft_most_stable_amorph["AB2"],
#     stoich="AB2")

# ab3_meta_lim = calc_dH(
#     raw_dft_most_stable_amorph["AB3"],
#     stoich="AB3")

# # calc_dH(stoich="AB2")

# + jupyter={}
# print("AB2:", ab2_meta_lim)
# print("AB3:", ab3_meta_lim)

# + jupyter={}
# raw_dft_most_stable_amorph = dict(
#     AB2=-6.542,
#     AB3=-6.163,
#     )
