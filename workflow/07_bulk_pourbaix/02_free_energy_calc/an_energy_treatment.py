# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
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

# # Bulk Energy Treatment for IrOx Systems
# ___
#
#
# TODO:
#   * Switch the bulk DFT energies to my numbers

# # Import Modules

# +
import os
import sys

sys.path.insert(0,
    os.path.join(
        os.environ["PROJ_irox"],
        "data",
        )
    )

# #############################################################################
from energetics.dft_energy import Element_Refs

# #############################################################################
# Gas phase molecules
from proj_data_irox import (
    h2o_ref,
    zpe_h2o,
    cv_h2o,
    ts_h2o,

    h2_ref,
    zpe_h2,
    cv_h2,
    ts_h2,

#     h2o_corr,
#     h2_corr,
    )

# Bulk DFT calculations
from proj_data_irox import (
    IrO2_ave_regressed_bulk_e,
    IrO3_ave_regressed_bulk_e,
    IrO3_rutile_like_ave_regressed_bulk_e,
    IrO3_battery_bulk_e_dft,    
    )
# -

# # Script Inputs

# +
kjmol = 96.485
T = 298.15

convert_dS = T / (1000 * kjmol)
# -

# # Gas Phase References

# +
e_h2o_r = h2o_ref
e_h2_r = h2_ref

h_h2o_r = e_h2o_r + zpe_h2o + cv_h2o
h_h2_r = e_h2_r + zpe_h2 + cv_h2

g_h2o_r = h_h2o_r + ts_h2o
g_h2_r = h_h2_r + ts_h2

# #############################################################################
Elem_Refs = Element_Refs(
    H2O_dict={
        "electronic_e": e_h2o_r, "enthalpy_e": h_h2o_r, "gibbs_e": g_h2o_r},
    H2_dict={
        "electronic_e": e_h2_r, "enthalpy_e": h_h2_r, "gibbs_e": g_h2_r},
    oxygen_ref="O2",
    hydrogen_ref="H2")

oxy_ref, hyd_ref = Elem_Refs.calc_ref_energies()
# -

print("O_ref =", oxy_ref.enthalpy_e)
print("H_ref =", hyd_ref.enthalpy_e)

# + [markdown] toc-hr-collapsed=false
# # IrOx Energy Treatment
# -

# ## IrOx Energy Quantities

# +
# IrO2_ave_regressed_bulk_e
# IrO3_ave_regressed_bulk_e
# IrO3_rutile_like_ave_regressed_bulk_e
# IrO3_battery_bulk_e_dft



# Ir  #########################################################################
e_ir_m = -8.860644725
s_ir_metal = 35.5  # exp.



# IrO2  #######################################################################
# e_iro2 = -7.045317 * 3  # calculated Rutile PBE energy, 600 eV
e_iro2 = IrO2_ave_regressed_bulk_e * 3
h_iro2 = e_iro2  # Assumping PV term is small and E ~= H

dh_iro2_exp = -242.672 / kjmol
dg_iro2_exp = -188.386 / kjmol

s_iro2_solid = 58.57  # exp.
s_iro2_solid_calc_phonon = 0.1485



# IrO3  #######################################################################
# e_iro3 = -6.442159 * 4  # calculated lowest PBE energy , 600 eV
# e_iro3 = -6.492707 * 4  # Taken from Chris's script
e_iro3 = IrO3_ave_regressed_bulk_e * 4
# e_iro3 = IrO3_rutile_like_ave_regressed_bulk_e * 4

h_iro3 = e_iro3  # Assumping PV term is small and E ~= H



# IrHO3  ######################################################################
e_irho3 = -6.123924 * 5  # calculated lowest PBE energy , 600 eV
h_irho3 = e_irho3  # Assumping PV term is small and E ~= H

# +
# print(S_iro2_solid * convert_dS, S_iro2_solid_calc_phonon)
# -

# ## Ir Metal Fit to Exp. IrO2 Formation E.

h_ir_m_fit = h_iro2 - (2 * -4.1811960) - dh_iro2_exp  # fit to exp dH
print(h_ir_m_fit)
# h_ir_m_fit = e_iro2 - (2 * -4.6579473) - dh_iro2_exp  # fit to exp dH

# + [markdown] toc-hr-collapsed=true
# ## Calculating dH_f of IrOx Species
# -

# ### IrO2

# +
dh_iro2 = h_iro2 - (2 * oxy_ref.enthalpy_e + h_ir_m_fit)

print(
    "\n",
    "IrO2: dH=", dh_iro2, " eV",
    "\n",
    "IrO2: dH=", dh_iro2 / 3, " eV/atom",
    "\n",
    "IrO2: dH=", dh_iro2 * kjmol, " kj/mol",
    "\n",
    "IrO2: dH_exp", dh_iro2_exp * kjmol, ' kj/mol',
    )
# -

# ### IrO3 (α-AlF3)

# +
dh_iro3 = h_iro3 - (3 * oxy_ref.enthalpy_e + h_ir_m_fit)

print(
    "\n",
    "IrO3_a-AlF3: dH=", dh_iro3, " eV",
    "\n",
    "IrO3_a-AlF3: dH=", dh_iro3 / 3, " eV/atom",
    "\n",
    "IrO3_a-AlF3: dH=", dh_iro3 * kjmol, " kj/mol",
    )
# -

# ### IrO3 (rutile-like)

# +
e_iro3_rutile = IrO3_rutile_like_ave_regressed_bulk_e * 4.
h_iro3_rutile = e_iro3_rutile

dh_iro3_rutile = h_iro3_rutile - (3 * oxy_ref.enthalpy_e + h_ir_m_fit)

print(
    "\n",
    "IrO3_rutile: dH=", dh_iro3_rutile, " eV",
    "\n",
    "IrO3_rutile: dH=", dh_iro3_rutile / 3, " eV/atom",
    "\n",
    "IrO3_rutile: dH=", dh_iro3_rutile * kjmol, " kj/mol",
    )
# -

# ### IrO3 (battery)

# +
e_iro3_battery = IrO3_battery_bulk_e_dft * 4.
h_iro3_battery = e_iro3_battery

dh_iro3_battery = h_iro3_battery - (3 * oxy_ref.enthalpy_e + h_ir_m_fit)

print(
    "\n",
    "IrO3_battery: dH=", dh_iro3_battery, " eV",
    "\n",
    "IrO3_battery: dH=", dh_iro3_battery / 3, " eV/atom",
    "\n",
    "IrO3_battery: dH=", dh_iro3_battery * kjmol, " kj/mol",
    )
# -

# ### IrHO3

# +
dh_irho3 = h_irho3 - (3 * oxy_ref.enthalpy_e + h_ir_m_fit + hyd_ref.enthalpy_e)

print(
    "\n",
    "IrHO3: dH=", dh_irho3, " eV",
    "\n",
    "IrHO3: dH=", dh_irho3 / 3, " eV/atom",
    "\n",
    "IrHO3: dH=", dh_irho3 * kjmol, " kj/mol",
    )

# + [markdown] toc-hr-collapsed=true
# ## Calculating dG_f of IrOx Species
# -

s_o2_gas = 205.2

# ### IrO2

# +
TdS_iro2 = (s_iro2_solid - s_ir_metal - s_o2_gas) * convert_dS
#print (S_iro2_solid-S_ir_metal-S_o2_gas)*T/(1000)

dg_iro2 = dh_iro2 - TdS_iro2
print(
    "\n",

    'IrO2: dG=',
    dg_iro2,
    ' eV',
    "\n",

    'IrO2: dG/3=',
    dg_iro2 / 3,
    ' eV/atom',
    "\n",

    'IrO2: dG/3=',
    dg_iro2 * kjmol,
    ' kj/mol',
    "\n",

    'IrO2_exp: dG/3=',
    dg_iro2_exp * kjmol,
    ' kj/mol',
    "\n",
    )
# -

# ### IrO3 (α-AlF3)

# +
factor = 3 / 2  # adjusted to reflect more O in IrO3 vs IrO2
TdS_iro3 = (s_iro2_solid * factor - s_ir_metal - 3 / 2 * s_o2_gas) * convert_dS
dg_iro3 = dh_iro3 - TdS_iro3

print(
    "\n",

    'IrO3: dG=',
    dg_iro3,
    ' eV',
    "\n",

    'IrO3: dG/3=',
    dg_iro3 / 3,
    ' eV/atom',
    "\n",

    'IrO3: dG=',
    dg_iro3 * kjmol,
    ' kj/mol',
    "\n",
    )
# -

# ### IrO3 (rutile-like)

# +
factor = 3 / 2  # adjusted to reflect more O in IrO3 vs IrO2
TdS_i = (s_iro2_solid * factor - s_ir_metal - 3 / 2 * s_o2_gas) * convert_dS
dg_iro3_rutile = dh_iro3_rutile - TdS_i

print(
    "\n",

    'IrO3: dG=',
    dg_iro3_rutile,
    ' eV',
    "\n",

    'IrO3: dG/3=',
    dg_iro3_rutile / 3,
    ' eV/atom',
    "\n",

    'IrO3: dG=',
    dg_iro3_rutile * kjmol,
    ' kj/mol',
    "\n",
    )
# -

# ### IrO3 (battery)

# +
factor = 3 / 2  # adjusted to reflect more O in IrO3 vs IrO2
TdS_i = (s_iro2_solid * factor - s_ir_metal - 3 / 2 * s_o2_gas) * convert_dS
dg_iro3_battery = dh_iro3_battery - TdS_i

print(
    "\n",

    'IrO3_battery: dG=',
    dg_iro3_battery,
    ' eV',
    "\n",

    'IrO3_battery: dG/3=',
    dg_iro3_battery / 3,
    ' eV/atom',
    "\n",

    'IrO3_battery: dG=',
    dg_iro3_battery * kjmol,
    ' kj/mol',
    "\n",
    )
# -

# ### IrHO3

# +
# # Adjusted to reflect more O H in IriHO3 vs. IrO2
# # Phonons needs to be calculated

# factor = 4 / 2
# TdS_irho3 = (S_iro2_solid * factor - S_ir_metal -
#     S_o2_gas * 3 / 2 - S_h2_gas / 2) * convert_dS
# #TdS_irho3=(-S_o2_gas*3/2-S_h2_gas/2)*convert_dS
# dg_irho3 = dh_irho3 - TdS_irho3
# print(
#     'IrHO3: dG= ',
#     dg_irho3,
#     ' eV ',
#     dg_irho3 / 5,
#     dg_irho3 * kjmol,
#     ' kjmol',
#     )
# -

# # Final Energies

# +
final_dict = {
    "Ir": 0.,
    "IrO2": dg_iro2,
    "IrO3": dg_iro3,
    "IrO3_rutile": dg_iro3_rutile,
    "IrO3_battery": dg_iro3_battery,
    }

print("ion_dict_solids_expt = " + str(final_dict))
