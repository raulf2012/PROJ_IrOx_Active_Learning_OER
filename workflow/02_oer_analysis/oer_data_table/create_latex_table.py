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

import os
print(os.getcwd())
import sys

# +
# %%capture

import pandas as pd
import numpy as np

# #########################################################
from IPython.display import display

# #########################################################
sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import data_dir
from proj_data_irox import corrections_dict

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "workflow"))
from an_data_processing import load_df

# #########################################################
from oxr_reaction.oxr_rxn import ORR_Free_E_Plot
from misc_modules.pandas_methods import drop_columns, reorder_df_columns
# -

pd.options.display.max_colwidth = 200

# # Read Data

# +
# path_i = "/mnt/f/Dropbox/01_norskov/PROJECT_DATA/04_IrOx_surfaces_OER/oer_slabs_results/190321_new_job_df/job_dataframe.pickle"

# # #########################################################
# import pickle; import os
# # with open(path_i, "rb") as fle:
# with open(path_i, "rb") as fle:
#     data = pickle.load(fle)
# # #########################################################

# +
# # %%capture

df_pourbaix, df_ads, df_surf = load_df(
    from_file=False,
    root_dir=data_dir,
    data_dir=data_dir,
    file_name="df_master.pickle",
    process_df=True)

df = df_ads
# -

df = drop_columns(df=df, keep_or_drop="keep",
    columns=[

        "bulk_system",
        "facet",
        "adsorbate",
        "coverage_type",
        # "ooh_direction",
        "ads_e",
        # "elec_energy",
        # "total_magmom",
        # "abs_magmom",
        # "path_short",
        # "name_i",
        # "max_force",
        # "sum_force",
        # "elem_num_dict",
        # "incar_parsed",
        # "init_atoms",
        # "atoms_object",
        # "N_atoms",
        # "dipole_correction",
        # "path",
        # "name_i_2",
        # "name_i_3",
        # "priority",
        "surface_type",

        ],
    )


# +
data_dict_list = []
grouped = df.groupby(["bulk_system", "facet", "coverage_type", "surface_type"])
for name, group in grouped:
    new_row_i = dict()
    # display(group)

    # #####################################################
    ORR_PLT = ORR_Free_E_Plot(
        free_energy_df=None,
        state_title="adsorbate",
        free_e_title="ads_e",
        # smart_format=smart_format_dict,
        color_list=None,
        rxn_type="OER")

    ORR_PLT.add_series(
        group,
        plot_mode="all",
        overpotential_type="OER",
        # property_key_list=prop_name_list,
        add_overpot=False,
        name_i="TEMP")

    oxr_series = ORR_PLT.series_list[0]
    overpot, lim_step = oxr_series.calc_overpotential_OER()

    # #####################################################
    # arrow_str = " \rightarrow "
    # if lim_step == ["oh", "o"]:
    #     lim_step = "$OHT"   + arrow_str + "OTT$"
    # elif lim_step == ["ooh", "bulk"]:
    #     lim_step = "$OOH"   + arrow_str + "*TTT$"
    # elif lim_step == ["o", "ooh"]:
    #     lim_step = "$OTT"   + arrow_str + "OOH$"
    # elif lim_step == ["bulk", "oh"]:
    #     lim_step = "$*TTT"   + arrow_str + "OHT$"
    # else:
    #     tmp = 42

    phan = "\phantom{T}"
    arrow_str = " \rightarrow "
    if lim_step == ["oh", "o"]:
        lim_step = "$*OH \phantom{T}"   + arrow_str + "*O \phantom{T} \phantom{T} $"
    elif lim_step == ["ooh", "bulk"]:
        lim_step = "$*OOH"   + arrow_str + "* \phantom{T} \phantom{T} \phantom{T} $"
    elif lim_step == ["o", "ooh"]:
        lim_step = "$*O \phantom{T} \phantom{T} "   + arrow_str + "*OOH$"
    elif lim_step == ["bulk", "oh"]:
        lim_step = "$* \phantom{T} \phantom{T} \phantom{T} "   + arrow_str + "*OH \phantom{T} $"
    else:
        tmp = 42

    # #####################################################
    new_row_i["overpot"] = overpot
    new_row_i["lim_step"] = lim_step

    # #####################################################
    ads_e_dict = dict(zip(
        group.adsorbate,
        group.ads_e,
        ))
    new_row_i.update(ads_e_dict)

    group = drop_columns(
        df=group,
        columns=["adsorbate", "ads_e"],
        keep_or_drop="drop",
        )

    other_props = dict()
    for column in group.columns:
        num_unique = group[column].unique().shape[0]
        if num_unique == 1:
            other_props[column] = group[column].iloc[0]

    new_row_i.update(other_props)

    data_dict_list.append(new_row_i)
    
df_new = pd.DataFrame(data_dict_list)

# + active=""
#
#
#
#
#
# -

# # Drop bare adsorption energy column

# +


df_new = drop_columns(
    df=df_new,
    columns=["bare"],
    keep_or_drop="drop",
    )
# -

# # Add new composite columns

# +
df_new["g_o-g_oh"] = df_new.o - df_new.oh

df_new["lim_pot"] = 1.23 + df_new["overpot"]

# #########################################################
df_new["e_oh"] = df_new["oh"] - corrections_dict["oh"]
df_new["e_o"] = df_new["o"] - corrections_dict["o"]
df_new["e_ooh"] = df_new["ooh"] - corrections_dict["ooh"]
# -

# # Sort column order

df_new = reorder_df_columns(
    [
        "bulk_system",
        "facet",
        "coverage_type",
        "surface_type",
        "e_oh", "e_o", "e_ooh",
        "oh", "o", "ooh",
        "g_o-g_oh",
        "lim_pot",
        "overpot",
        "lim_step",
        ],
    df_new,
    )

# # Round float columns

# +
num_dec = 3

df_new = df_new.round({
    "oh": num_dec, "o": num_dec, "ooh": num_dec,
    "e_oh": num_dec, "e_o": num_dec, "e_ooh": num_dec,
    "g_o-g_oh": num_dec,
    "overpot": num_dec, "lim_pot": num_dec,
    })
# -

# # Format facet column and change b-IrO3 entries

# +
df_new.facet = "(" + df_new.facet + ")"

ind_a = df_new[
    (df_new.bulk_system == "IrO3_battery") & \
    (df_new.facet == "(010)") & \
    (df_new.surface_type == "a")
    ].index[0]

ind_b = df_new[
    (df_new.bulk_system == "IrO3_battery") & \
    (df_new.facet == "(010)") & \
    (df_new.surface_type == "b")
    ].index[0]

df_new.at[ind_a, "facet"] = "(010)-A"
df_new.at[ind_b, "facet"] = "(010)-B"

df_new = drop_columns(df=df_new, columns="surface_type", keep_or_drop="drop")
# -

# # Sort based on bulk_system, facet, and coverage_type

df_new = df_new.sort_values(["bulk_system", "facet", "coverage_type"])

units_dict = {
    "bulk_system": "-",
    "facet": "-",
    "coverage_type": "-",
    "e_oh": "(eV)",
    "e_o": "(eV)",
    "e_ooh": "(eV)",
    "oh": "(eV)",
    "o": "(eV)",
    "ooh": "(eV)",
    "g_o-g_oh": "(eV)",
    "lim_pot": "(V)",
    "overpot": "(V)",
    "lim_step": "-",
    }


# +
# assert False
# -

# # Rename columns and entries

# +
# #########################################################
# #########################################################
column_rename_dict = {
    "bulk_system": "Bulk Sys.",
    "facet": "Facet",
    "coverage_type": "Coverage",
    "surface_type": "SurfaceTMP",

    "g_o-g_oh": "$\\Delta G_{O}-\\Delta G_{OH}$",

    "oh": "$\\Delta G_{OH}$",
    "o": "$\\Delta G_{O}$",
    "ooh": "$\\Delta G_{OOH}$",

    
    "e_oh": "$\\Delta E_{OH}$",
    "e_o": "$\\Delta E_{O}$",
    "e_ooh": "$\\Delta E_{OOH}$",
    
    "overpot": "$\\eta$",
    "lim_pot": "Lim. Pot.",
    "lim_step": "RDS",
    }
column_rename_dict_inv = {v: k for k, v in column_rename_dict.items()}

df_new = df_new.rename(columns=column_rename_dict)

# #########################################################
df_new = df_new.replace(
    to_replace="IrO2",
    value="$R{\text -}IrO_{2}$")
    # value="$R-IrO_{2}$")

df_new = df_new.replace(
    to_replace="IrO3",
    value="$\\alpha{\text -}IrO_{3}$")

df_new = df_new.replace(
    to_replace="IrO3_battery",
    value="$\\beta{\text -}IrO_{3}$")
    # value="\beta-IrO_{3}")

df_new = df_new.replace(
    to_replace="IrO3_rutile-like",
    value="$R{\text -}IrO_{3}$")

df_new = df_new.replace(
    to_replace="NaN",
    value=np.nan)

# #########################################################
df_new = df_new.replace(
    to_replace="h_covered",
    value="*OH")
df_new = df_new.replace(
    to_replace="o_covered",
    value="*O")
df_new = df_new.replace(
    to_replace="o_covered_2",
    value="*O-partial")
# -

# # Adding units index

# +
units_list = []
for col_i in df_new.columns.values:
    unit_i = units_dict[
        column_rename_dict_inv[col_i]
        ]
    units_list.append(unit_i)


tuples = list(zip(*[
    df_new.columns.values,
    units_list,
    ]))


index = pd.MultiIndex.from_tuples(tuples, names=["Header", "Units"])
df_new.columns = index
# -

# # Column Alignment

# +
alignment_column_dict = {
    "Bulk Sys.": "l",
    "Facet": "l",
    "Coverage": "l",
    "$\\Delta E_{OH}$": "c",
    "$\\Delta E_{O}$": "c",
    "$\\Delta E_{OOH}$": "c",
    "$\\Delta G_{OH}$": "c",
    "$\\Delta G_{O}$": "c",
    "$\\Delta G_{OOH}$": "c",
    "$\\Delta G_{O}-\\Delta G_{OH}$": "c",
    "Lim. Pot.": "c",
    "$\\eta$": "c",
    "RDS": "c",
    }

alignment_list = [alignment_column_dict[i[0]] for i in df_new.columns.tolist()]

alignment_str = ""
for i in alignment_list:
    alignment_str += i
# -

# # Write Dataframe to Latex Table

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
df_new.to_latex(
    buf="oer_table.tex",
    **shared_props)

path_i = os.path.join(
    os.environ["PROJ_irox_paper"],
    "04_data_tables/oer_energetics",
    "oer_table.tex")
df_new.to_latex(
    buf=path_i,
    **shared_props)
# -

df_new
