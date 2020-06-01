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

# # Import Modules

# +
import os
print(os.getcwd())
import sys

import pandas as pd

# #########################################################
sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import calc_dH

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "workflow/ml_modelling"))
from ml_methods import get_ml_dataframes
# -

# # Script Inputs

# +
TdS_IrO2 = -0.562803124838058
TdS_IrO3 = -0.789056096134837

mu_H2O = -2.4583

# #########################################################
dH_r_iro2 = -2.515126703632689
dH_b_iro3 = 4 * -0.596831474814422

# #########################################################
dG_ir_aq = -2.038081775

# +
all_systems = [
    "ir",
    "r-iro2",
    "a-iro3",
    "r-iro3",
    "b-iro3",
    "ir_aq",
    ]
    
TdS_dict = {
    "AB2": TdS_IrO2,
    "AB3": TdS_IrO3,
    }
# -

# # Read Data

# ## Read DFT data for a-IrO3 and R-IrO3

# +
DF_dict = get_ml_dataframes()

df_dft = DF_dict.get("df_dft_final_final")

df_i = df_dft[
    (df_dft.stoich == "AB3")
    ].sort_values("dH")

alpha_iro3_row = df_i.iloc[0]
rutile_iro3_row = df_i.loc["b5cgvsb16w"]
# -

# # Put together Enthalpy dict

# +
dH_dict = {
    "ir": 0.,
    "r-iro2": dH_r_iro2,
    "a-iro3": 4 * alpha_iro3_row.dH,
    "r-iro3": 4 * rutile_iro3_row.dH,
    "b-iro3": dH_b_iro3,

    # "ir_aq": None,
    }

# dH_dict
# -

# # TEMP

# +
dG_dict = dict()

for sys_i, dH_i in dH_dict.items():

    if "iro2" in sys_i:
        stoich_i = "AB2"
    elif "iro3" in sys_i:
        stoich_i = "AB3"
    else:
        continue

    # Getting Entropy correction
    TdS_i = TdS_dict[stoich_i]


    dG = dH_i - TdS_i
    dG_dict[sys_i] = dG


# TdS_IrO2 = -0.562803124838058
# TdS_IrO3 = -0.789056096134837

dG_dict["ir_aq"] = dG_ir_aq
dG_dict["ir"] = 0.

# +
data_dict_list = []
for sys_i in all_systems:
    # print(sys_i)

    data_dict_i = dict()
    data_dict_i["system"] = sys_i
    
    
    data_dict_i["dH"] = dH_dict.get(sys_i, None)
    data_dict_i["dG"] = dG_dict.get(sys_i, None)
    
    data_dict_list.append(data_dict_i)

df = pd.DataFrame(data_dict_list)
df

# + active=""
#
#
#
#
#
#
#

# +
num_dec = 3

df = df.round({
    "dH": num_dec,
    "dG": num_dec,
    })

# +
units_dict = {
    "system": "",
    "dH": "(eV/f.u.)",
    "dG": "(eV/f.u.)",
    }

units_list = []
for col_i in df.columns.values:
    unit_i = units_dict[col_i]
    units_list.append(unit_i)

tuples = list(zip(*[
    df.columns.values,
    units_list,
    ]))

index = pd.MultiIndex.from_tuples(
    tuples,
    # names=["Header", "Units"],
    )
df.columns = index

# +
column_rename_dict = {
    "dH": "$\Delta H_{f}$",
    "dG": "$\Delta G_{f}$",
    }

df = df.rename(columns=column_rename_dict)


df = df.replace(
    to_replace="ir",
    value="Ir(s)")

df = df.replace(
    to_replace="r-iro2",
    value="$R$-$IrO_{2}(s)$")

df = df.replace(
    to_replace="a-iro3",
    value="$\\alpha$-$IrO_{3}(s)$")
df = df.replace(
    to_replace="r-iro3",
    value="$R$-$IrO_{3}(s)$")
df = df.replace(
    to_replace="b-iro3",
    value="$\\beta$-$IrO_{3}(s)$")

df = df.replace(
    to_replace="ir_aq",
    value="$IrO_{4}^{-}(aq)$")
# -

alignment_str = "lcc"

# + active=""
#
#
#
#
#
# -

# # Write latex table to file

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

# +
df.to_latex(
    buf="bulk_pourb_table.tex",
    **shared_props)

path_i = os.path.join(
    os.environ["PROJ_irox_paper"],
    # "PAPER_IrOx_Active_Learning_OER",
    "04_data_tables/bulk_pourb_energy",
    "bulk_pourb_table.tex")

df.to_latex(
    buf=path_i,
    **shared_props)
# -

df

# + active=""
#
#
#
#

# + jupyter={}
# # /mnt/f/Dropbox/01_norskov/00_git_repos
# "PAPER_IrOx_Active_Learning_OER"
# "04_data_tables/bulk_pourb_energy"

# + jupyter={}
# df_tmp = df.copy()
# df_tmp = df_tmp.set_index(('system', ''))

# + jupyter={}
# # df_tmp = df


# df = df.set_index(
#     'system',
#     drop=True,
#     append=False,
#     inplace=False,
#     verify_integrity=False,
#     )

# + jupyter={}

# df_tmp.columns.tolist()

# df_tmp

# + jupyter={}
# assert False
