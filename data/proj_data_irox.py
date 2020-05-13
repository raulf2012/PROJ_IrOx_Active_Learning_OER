#!/usr/bin/env python

"""Misc storage location for data.

Author: Raul A. Flores


Setting up paths to read from this file:

import os
import sys
sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import calc_dH
"""

# | - Import Modules
import os
import sys
#__|

tmp = 42

# | - VASP Gas-phase References
# Exp. heat of formation of H2O
DG_f_H2O = -2.4583  # eV   -4.9166 eV for 2 H2O molecules

h2_ref = -6.759300

# My Calculations
# /global/cscratch1/sd/flores12/IrOx_Project/06_gas_phase/01-att
h2o_ref = -14.21890762

# Recent gas-phase references from Michal
h2_ref = -6.77149190
h2o_ref = -14.23091949


# 500 eV, Normal
gas_dft_references_dict = {
    500: dict(
        h2o=-14.23091949,  # PBE, 500 eV
        h2=-6.77149190,  # PBE, 500 eV
        o2=-9.88557216,  # PBE, 500 eV
        ),

    600: dict(
        h2o=-14.23835553,  # PBE, 600 eV
        h2=-6.77402325,  # PBE, 600 eV
        o2=-9.89149843,  # PBE, 600 eV
        )
    }



# | - Free Energy Corrections
#Contributions to Gibbs Energies for gas molecules
# (VASP-PBE calculated by Max; T= 300K)

# zpeh2o=0.560    #exp. NIST 0.558
# zpeh2=0.268     #exp. NIST 0.273
# cvh2o=0.103      #From Colin at P = 0.035 bar
# cvh2=0.0905
# tsh2o=0.675     #From Colin at P = 0.035 bar
# tsh2=0.408     #at P = 1bar

# "Standard"
zpe_h2o = 0.56
cv_h2o = 0.103
ts_h2o = 0.675

zpe_h2 = 0.268
cv_h2 = 0.091
ts_h2 = 0.408

h2o_corr = zpe_h2o + cv_h2o - ts_h2o
h2_corr = zpe_h2 + cv_h2 - ts_h2

# h2o_corr = zpe_h2o - ts_h2o
# h2_corr = zpe_h2 - ts_h2
#__|

#__|


# | - Adsorbate Free Energy Corrections
# Adsorbate Vibrational Analysis --> Free Energy Contributions
# Obtained from DFT jobs in the following dir:
# /global/cscratch1/sd/flores12/IrOx_Project/04_ads_vib/IrO3/110/01_O_covered
ads_fe_dict = {
    "oh": {
        "zpe": 0.354,
        "cv": 0.053,
        "ts": 0.1,
        },

    "o": {
        "zpe": 0.093,
        "cv": 0.023,
        "ts": 0.035,
        },

    "ooh": {
        "zpe": 0.464,
        "cv": 0.058,
        "ts": 0.096,
        },
    }

# | - NEW | Attempt to automate the calculation of the adsorbate FE corr dict
# 0.42611961567987805


# def calc_ads_corr_i(
#     ads_spec,
#     ads_fe_dict,
#     h2o_corr,
#     h2_corr,
#     ):
#     """
#     """
#     # | - calc_ads_corr_i
#     if ads_spec == "oh":
#         num_H = 1
#         num_O = 1
#     elif ads_spec == "o":
#         num_H = 0
#         num_O = 1
#     elif ads_spec == "ooh":
#         num_H = 1
#         num_O = 2
#
#     corr_i = (0. +
#
#         + ( 0. +
#             + ads_fe_dict[ads_spec]["zpe"] +
#             + ads_fe_dict[ads_spec]["cv"] +
#             - ads_fe_dict[ads_spec]["ts"]
#             ) +
#
#         - (0. +
#             + (num_O) * h2o_corr +
#             + ((num_H - num_O * 2.) / 2.) * h2_corr
#             )
#
#         )
#
#     return(corr_i)
#     #__|
#
# corr_oh = calc_ads_corr_i("oh", ads_fe_dict, h2o_corr, h2_corr)
# corr_o = calc_ads_corr_i("o", ads_fe_dict, h2o_corr, h2_corr)
# corr_ooh = calc_ads_corr_i("ooh", ads_fe_dict, h2o_corr, h2_corr)
#
# corrections_dict_tmp = {
#     "ooh": corr_ooh,
#     "o": corr_o,
#     "oh": corr_oh,
#     "bare": 0.,
#     }
#__|

# My data
corrections_dict = {
    "ooh": 0.3765,
    "o": 0.044,  # I had this as negative before | RF - 181105
    "oh": 0.2945,
    "bare": 0.,
    }
#__|


# | - Bulk Calculations
# Bulk energies from typical bulk DFT calculation
# Energies are on a per atom basis
# Regressed bulk energies from surface energy extrapolations

# | - IrO2
# $PROJ_irox/workflow/an_bulk_systems/r_iro2/_1

# DFT bulk energy
IrO2_bulk_e_dft = -7.044164135

# Regressed bulk energy
IrO2_ave_regressed_bulk_e = -7.044138074749999


# Potential Energy:
# -42.26498481
# Ir2O4

# Potential Energy:
# -42.26498481
#
# Atoms.get_chemical_formula:
# Ir2O4
#__|

# | - IrO3 (aAlF3)
# DFT bulk energy
IrO3_bulk_e_dft = -6.46018589875

# Regressed bulk energy
IrO3_ave_regressed_bulk_e = -6.463877450346546

# Potential Energy:
# -51.68148719
#
# Atoms.get_chemical_formula:
# Ir2O6
#__|

# | - IrO3 (rutile)
# DFT bulk energy
IrO3_rutile_like_bulk_e_dft = -6.457128223125

# Regressed bulk energy
IrO3_rutile_like_ave_regressed_bulk_e = -6.42967135329762

# Potential Energy:
# -103.31405157
#
# Atoms.get_chemical_formula:
# Ir4O12
#__|

# | - IrO3 (battery)
# DFT bulk energy
IrO3_battery_bulk_e_dft = -6.38668709984375

# Potential Energy:
# -408.74797439
#
# Atoms.get_chemical_formula:
# Ir16O48

# Same DFT settings as other systems

# Potential Energy:
# -410.62404231 eV
# -6.41600066109375 eV/atom


# -410.61718419 eV
# -6.41589350296875 eV/atom

# -410.61653865 eV
# -6.41588341640625 eV/atom

# -410.50375054 eV
# -6.4141211021875 eV/atom
#__|

bulk_e_per_atom_dict = {
    "IrO2": IrO2_ave_regressed_bulk_e,
    "IrO3": IrO3_ave_regressed_bulk_e,
    "IrO3_rutile-like": IrO3_rutile_like_ave_regressed_bulk_e,
    "IrO3_battery": IrO3_battery_bulk_e_dft,
    }
#__|


#| - Experimental Thermochemical Data
# dh_iro2_exp = -242.672 / kjmol  # Barin
dh_iro2_exp = -2.515126703632689

# dg_iro2_exp = -188.386 / kjmol  # Barin
dg_iro2_exp = -1.9524900243561174

#__|


# | - Computing dH and dG for IrO2 and IrO3
def calc_dH(
    e_per_atom,
    stoich=None
    ):
    """
    Based on a E_DFT/atom of -7.047516 for rutile-IrO2
    
    See the following dir for derivation:
        PROJ_IrOx_Active_Learning_OER/workflow/energy_treatment_deriv/calc_references
    """
    # | - calc_dH
    o_ref = -4.64915959
    ir_metal_fit = -9.32910211636731

    # o_ref = -4.657947279999998
    # ir_metal_fit = -9.316164736367316
    # iro2: -21.147186
    # o_ref2: -4.657947279999998
    # ir_metal_fit: -9.316164736367316
    # IrO2 (exp) -2.515126703632689

    if stoich == "AB2":
        dH = (2 + 1) * e_per_atom - 2 * o_ref - ir_metal_fit
        dH_per_atom = dH / 3.
    elif stoich == "AB3":
        dH = (3 + 1) * e_per_atom - 3 * o_ref - ir_metal_fit
        dH_per_atom = dH / 4.

    return(dH_per_atom)
    #__|

#__|


# | - Bulk Pourbaix Transitions
# From my bulk Pourbaix plot | RF | 190104
bulk_pourb_trans_dict = {
    "ir_iro2": 0.74107,
    "iro2_iro3": 1.07148,
    "iro3_iro4-": 2.33329,

    "IrO3_battery_low": 1.225857,
    "IrO3_battery_high": 2.024,
    "IrO3_rutile-like_low": 1.13989,
    "IrO3_rutile-like_high": 2.196,
    }
#__|


# | - Metastability limits of IrO2 and IrO3
metastability_limits_raw_dft = dict(
    AB2=-6.542,
    AB3=-6.163,
    )

# Computed in
# PROJ_IrOx_Active_Learning_OER/workflow/metastability_limit_murat
# AB2: -0.33285956787756277
# AB3: -0.3438547784081729
#__|


# | - Color Palettes
irox_bulk_color_map = {
    "Ir": "#6799A3",
    "IrO2": "#7FC97F",
    "IrO3": "#BEAED4",
    "IrO3_a-AlF3": "#BEAED4",
    "IrO3_rutile-like": "#FDC086",
    # "IrO3_battery": "#FFFF99",
    "IrO3_battery": "#ff9dcd",
    # "HIrO3": "red",

    # "IrO4-": "pink",
    # "IrO4-": "#FFC3CB",
    # "IrO4-": "#7c7c7c",
    # "IrO4-": "#9c9c9c",
    "IrO4-": "#3a3a3a",

    }

pymatgen_to_my_naming_convention = {
    # "Ir": "Ir",
    # "IrO2": "IrO2",
    # "IrO3": "IrO3",
    # "IrO3_a-AlF3": "NEWNEW",
    "HIrO3(s)": "HIrO3",
    # "IrO4-": "IrO4-",
    }


irox_surface_e_color_map = {
    "IrO2_bare": "#7fc97f",
    "IrO2_h_covered": "#57b557",
    "IrO2_o_covered": "#3d933d",
    "IrO3_bare": "#beaed4",
    "IrO3_h_covered": "#9f87bf",
    "IrO3_o_covered": "#7f5fab",
    "IrO3_rutile-like_bare": "#fdc187",
    "IrO3_rutile-like_h_covered": "#fca24c",
    "IrO3_rutile-like_o_covered": "#fa7d06",
    "IrO3_battery_bare": "#ff9dcd",
    "IrO3_battery_h_covered": "#ff61ae",
    "IrO3_battery_o_covered": "#ff1b8c",
    "IrO3_battery_half_o_covered": "#ff1b8c",
    }
# __|


# | - Surface Energies
surface_energies = {
    "IrO2_100": 0.22334409161204058,
    "IrO2_110": 0.16552537099343662,

    "IrO3_100": 0.25584134125701796,
    "IrO3_110": 0.15491737507750267,
    "IrO3_111": 0.272534627288274,
    "IrO3_211": 0.2177463607996281,

    "IrO3_rutile-like_001": 0.170714237243587,
    "IrO3_rutile-like_100": 0.16401138601124265,
    "IrO3_rutile-like_110": 0.10704959112203366,
    # "IrO3_rutile-like_110": 1.10704959112203366,
    }

max_surf_e = max(list(surface_energies.values()))
min_surf_e = min(list(surface_energies.values()))
# __|


# | - Smart Format Dict
smart_format_dict = [

    [
        {"bulk_system": "IrO2"},
        {"color2": irox_bulk_color_map["IrO2"]}],
    [
        {"bulk_system": "IrO3"},
        {"color2": irox_bulk_color_map["IrO3"]}],
    [
        {"bulk_system": "IrO3_rutile-like"},
        {"color2": irox_bulk_color_map["IrO3_rutile-like"]}],
    [
        {"bulk_system": "IrO3_battery"},
        {"color2": irox_bulk_color_map["IrO3_battery"]}],

    [{"coverage_type": "o_covered"}, {"symbol": "circle"}],
    [{"coverage_type": "o_covered_2"}, {"symbol": "circle"}],
    [{"coverage_type": "h_covered"}, {"symbol": "triangle-up"}],

    ]

smart_format_dict_volcano = [
    [{"bulk_system": "IrO3"}, {"color2": "black"}],
    [{"bulk_system": "IrO2"}, {"color2": "blue"}],

    [{"coverage_type": "o_covered"}, {"symbol": "o"}],
    [{"coverage_type": "h_covered"}, {"symbol": "^"}],

    [{"coverage_type": "O-4_OH-0"}, {"symbol": "o"}],
    [{"coverage_type": "O-2_OH-0"}, {"symbol": "o"}],
    [{"coverage_type": "O-2_OH-2"}, {"symbol": "^"}],

    [{"facet": "110"}, {"color1": "red"}],
    [{"facet": "211"}, {"color1": "green"}],
    [{"facet": "100"}, {"color1": "black"}],
    [{"facet": "111"}, {"color1": "pink"}],

    ]

smart_format_dict_FED = [

    [
        {"bulk_system": "IrO3"},
        {"color": "red"}, ],

    [
        {"bulk_system": "IrO2"},
        {"color": "blue"},
        ],

    [
        {"coverage_type": "o_covered"},
        {"dash": "dot"},
        ],

    [
        {"coverage_type": "h_covered"},
        {"dash": None},
        ],
    ]

#__|


#| - OER Systems to plot
oer_systems_to_plot = [
    "IrO2_100_h_covered_NaN",
    "IrO2_100_o_covered_NaN",
    "IrO2_110_h_covered_NaN",
    "IrO2_110_o_covered_NaN",
    "IrO3_100_h_covered_NaN",
    "IrO3_100_o_covered_NaN",
    "IrO3_110_h_covered_NaN",
    "IrO3_110_o_covered_NaN",
    "IrO3_111_o_covered_NaN",
    #  "IrO3_211_h_covered_NaN",
    "IrO3_211_o_covered_NaN",
    "IrO3_battery_010_o_covered_a",
    "IrO3_battery_010_o_covered_b",
    "IrO3_rutile-like_100_h_covered_NaN",
    "IrO3_rutile-like_100_o_covered_NaN",
    "IrO3_rutile-like_100_o_covered_2_NaN",
    #  "IrO3_rutile-like_110_h_covered_NaN",
    "IrO3_rutile-like_110_o_covered_NaN",

    ]
#__|


# | - Scaling Relation Data
g_h2o = 0
g_o2 = 4.92

gas_molec_dict = {
    "h2": 0.,
    "o2": 4.92,
    "h2o": 0,
    }

# Ideal scaling
scaling_dict_ideal = {

    "ooh": {
        "m": 1.,
        "b": 3.2,
        },

    "o": {
        "m": 2.,
        "b": 0.,
        },

    "oh": {
        "m": 1.,
        "b": 0.,
        },

    }


scaling_dict_fitted = {

    "ooh": {
        "m": 0.976,
        "b": 3.09,
        },
    "o": {
        "m": 1.30,
        "b": 1.19,
        },
    "oh": {
        "m": 1.,
        "b": 0.,
        },
    }

#__|


# | - Experimental IrOx Data From Nature Paper (R14)
exp_irox_lim_pot_10mA = {
    # R14
    "SrIrO3": 1.51,
    "IrO2(110)": 1.94,

    # R95
    "IrO2(110)_R95": 1.8,
    }

exp_irox_lim_pot_1mA = {
    "SrIrO3": 1.45,
    "IrO2(110)": 1.8}


exp_irox_lim_pot = {
    "10_mA/cm2": exp_irox_lim_pot_10mA,
    "1_mA/cm2": exp_irox_lim_pot_1mA,
    "IrOx": 1.57,
    # "RuOx": ,
    }
#__|


# | - Figure Settings

# | - Fonts styling
axis_label_font_size = 10 * (4 / 3)
# axis_tick_labels_font_size = 9 * (4 / 3)
axis_tick_labels_font_size = 8 * (4 / 3)

font_family = "Arial"
base_font_color = "black"
#__|

# How we'll write out the voltage axis labels
voltage_name = "U<sub>RHE</sub> (V)"

#__|


# | - __misc__
proj_dir_name = "04_IrOx_surfaces_OER"

data_dir_name = "04_IrOx_surfaces_OER"

# main_systems = ["IrO2", "IrO3", "IrO3_battery", "IrO3_rutile-like"]
main_systems = ["IrO2", "IrO3", "IrO3_rutile-like", "IrO3_battery"]

data_dir = os.path.join(
    os.environ["PROJ_DATA"],
    data_dir_name,
    "oer_slabs_results",
    "190321_new_job_df",
    )

system_names_dict = {
    "IrO2": "R-IrO<sub>2</sub>",
    "IrO3": "α-IrO<sub>3</sub>",
    "IrO3_rutile-like": "R-IrO<sub>3</sub>",
    "IrO3_battery": "β-IrO<sub>3</sub>"

    #  "IrO2": "IrO<sub>2</sub> (rutile)",
    #  "IrO3": "IrO<sub>3</sub> (a-AlF<sub>3</sub>)",
    #  "IrO3_rutile-like": "IrO<sub>3</sub> (rutile)",
    #  "IrO3_battery": "IrO<sub>3</sub> (battery)"
    }

# The following properties split the OER ads data into sets of bare,O,OH,OOH
groupby_props = [
    "bulk_system",
    "facet",
    "coverage_type",
    "surface_type",
    ]
#__|


# | - Location of important data files
# import os; import sys
# sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
# from proj_data_irox import (
#     bulk_dft_data_path,
#     unique_ids_path,
#     prototypes_data_path,
#     static_irox_structures_path,
#     oqmd_irox_data_path,
#     )

# Processed bulk data (No OQMD data here)
bulk_dft_data_path = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling",
    "processing_bulk_dft/out_data",
    "df_bulk_dft.pickle")

unique_ids_path = os.path.join(
    os.environ["PROJ_irox"],
    "data/ml_irox_data",
    "unique_ids.csv")

# /home/raulf2012/Dropbox/01_norskov/00_git_repos/PROJ_IrOx_Active_Learning_OER

prototypes_data_path = os.path.join(
    os.environ["PROJ_irox"],

    "workflow/ml_modelling/processing_bulk_dft",
    "static_prototypes_structures/out_data",

    # "workflow/ml_modelling",
    # "static_prototypes_structures/out_data",

    "data_prototypes.pickle",

    # "chris_prototypes_structures",
    # "out_data",
    # "data_prototypes.pickle",

    )


# $dropbox/01_norskov/00_git_repos/PROJ_IrOx_Active_Learning_OER/workflow/ml_modelling/static_prototypes_structures
# workflow/ml_modelling/static_prototypes_structures
static_irox_structures_path = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/processing_bulk_dft",
    "static_prototypes_structures/out_data",
    "data_structures.pickle",

    # "workflow/ml_modelling",
    # "workflow/ml_modelling/processing_bulk_dft/static_prototypes_structures/out_data",
    # "static_prototypes_structures/out_data",
    )

static_irox_structures_kirsten_path = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/processing_bulk_dft",
    "static_prototypes_structures/out_data",

    # "workflow/ml_modelling",
    # "static_prototypes_structures/out_data",
    "data_structures_kirsten.pickle",
    )
oqmd_irox_data_path = os.path.join(
    os.environ["PROJ_irox"],
    # "workflow/ml_modelling/create_oqmd_data_df/out_data",
    "workflow/ml_modelling/processing_bulk_dft",
    "parse_oqmd_data/out_data",
    "df_oqmd_data.pickle")


# Fingerprints ################################################################
# TODO COMBAK
# /home/raulf2012/Dropbox/01_norskov/00_git_repos/PROJ_IrOx_Active_Learning_OER

fp_base_path = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/voronoi_featurize/out_data",

    # "workflow/ml_modelling/00_ml_workflow",
    # "190611_new_workflow/01_data_coll_feat/out_data",
    )

df_features_pre_opt_path = os.path.join(
    fp_base_path, "df_features_pre_opt.pickle")
df_features_pre_opt_kirsten_path = os.path.join(
    fp_base_path, "df_features_pre_opt_kirsten.pickle")

df_features_post_opt_path = os.path.join(
    fp_base_path, "df_features_post_opt.pickle")


# Small dataframe of OER bulk structures
oer_bulk_structures_path = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/processing_bulk_dft",
    "parse_my_oer_bulk_dft/out_data",
    "df_oer_bulk.pickle",
    )


# #############################################################################
# CCF Analysis ################################################################
df_ccf_path = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/ccf_similarity_analysis",
    "compute_ccf_and_dij_matrix/out_data",
    "df_ccf.pickle")

df_dij_path = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/ccf_similarity_analysis",
    "compute_ccf_and_dij_matrix/out_data",
    "df_d_ij_all.pickle")

# /mnt/c/Users/raulf2012/Dropbox/01_norskov/00_git_repos/PROJ_IrOx_Active_Learning_OER

ids_to_discard__too_many_atoms_path = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/processing_bulk_dft/static_prototypes_structures/out_data",
    # "workflow/ml_modelling/static_prototypes_structures/out_data",
    "ids_to_discard__too_many_atoms.pickle"
    )


ids_duplicates_path = os.path.join( 
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow/191102_new_workflow/get_duplicates_from_al",
    "out_data/duplicates.pickle")


# Prototype Classification Analysis

df_prototype_dft_path = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/processing_bulk_dft/symmetry_analysis_dft_static",
    "out_data",
    "df_prototype_dft.pickle")

df_prototype_static_path = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/processing_bulk_dft/symmetry_analysis_dft_static",
    "out_data",
    "df_prototype_static.pickle")

df_dft_final_final_path = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/processing_bulk_dft/creating_final_dataset_for_upload",
    "out_data/df_dft_final_no_dupl.pickle")

#  /mnt/f/Dropbox/01_norskov/00_git_repos/PROJ_IrOx_Active_Learning_OER


# | - __old__
# voronoi_features_data_path = os.path.join(
#     os.environ["PROJ_irox"],
#     "workflow/ml_modelling/00_ml_workflow",
#     "190611_new_workflow/01_data_coll_feat/out_data",
#     "df_features_pca.pickle",
#     )
#
# voronoi_features_all_data_path = os.path.join(
#     os.environ["PROJ_irox"],
#     "workflow/ml_modelling/00_ml_workflow",
#     "190611_new_workflow/01_data_coll_feat/out_data",
#     "df_features.pickle",
#     )
#__|

#__|


# | - METHODS
def get_relative_path_to_proj(path):
    relative_path_to_proj = path.replace(
        os.environ["PROJ_irox"], "")

    relative_path_to_proj
    if relative_path_to_proj[0] == "/":
        relative_path_to_proj = relative_path_to_proj[1:]

    return(relative_path_to_proj)

#__|


#| - AL Settings and Variables | OLD
sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/00_ml_workflow/191102_new_workflow"))
#__|

