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

# # Bulk Energy Treatment for IrOx Systems
# ---

# # Import Modules

# +
import os
print(os.getcwd())
import sys

sys.path.insert(0,
    os.path.join(
        os.environ["PROJ_irox"],
        "data"))

# #############################################################################
from energetics.dft_energy import Element_Refs

# #############################################################################
# Gas phase molecules
from proj_data_irox import (
    zpe_h2o,
    cv_h2o,
    ts_h2o,

    zpe_h2,
    cv_h2,
    ts_h2)

# #############################################################################
from proj_data_irox import gas_dft_references_dict
gas_dft_refs_i = gas_dft_references_dict[500]

h2o_ref = gas_dft_refs_i["h2o"]
h2_ref = gas_dft_refs_i["h2"]
# -

# # Script Inputs

-7.047516 * 6

# +
T = 298.15

# DFT Quantities
dft_energy_dict = {
    # FROM AL calculations at 500 eV
    "iro2": -7.047516,
    # "iro2": -7.047426,

    "iro3": -6.469847,
    # "iro3": -6.467450,
    }


# #############################################################################
# IrO2 Experimental thermochemical data #######################################

# # dh_iro2_exp = -242.672 / kjmol  # Barin
# dh_iro2_exp = -2.515126703632689
# # dg_iro2_exp = -188.386 / kjmol  # Barin
# dg_iro2_exp = -1.9524900243561174

from proj_data_irox import (
    dg_iro2_exp,
    dh_iro2_exp)

# #############################################################################
# Entropies ###################################################################

# Reference for IrO2 thermochemical data
# https://onlinelibrary.wiley.com/doi/book/10.1002/9783527619825
# Thermochemical Data of Pure Substances, Third Edition

# rutile-IrO2
S_iro2_solid = 58.576  # J/(mol K)  # Barin
TS_iro2_solid = 0.18098818987407367  # eV


# Ir metal
S_ir_metal = 35.505  # J/(mol K)  # Barin
TS_ir_metal = 0.10969917603772607  # eV


# Oxygen gas entropy reference
# https://webbook.nist.gov/cgi/cbook.cgi?ID=C7782447&Mask=1
# CODATA key values for thermodynamics
# https://www.worldcat.org/title/codata-key-values-for-thermodynamics/oclc/18559968
S_o2_gas = 205.2  # J/(mol K)  # CODATA thermo book
TS_o2_gas = 0.6340921386744053  # eV
# -

# # Gas Phase References

print(h2o_ref)
h2_ref

# +
# # %%capture

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

# +
oxy_ref

# e_h2o_r
# h_h2o_r

# +
# assert False
# -

# # Fitting the Ir metal reference to the experimental formation energy of IrO2

# +
e_iro2_pa = dft_energy_dict["iro2"]
e_iro2 = 3 * e_iro2_pa

PV_i = 0.
h_iro2 = e_iro2 + PV_i

oxy_ref_h = oxy_ref.enthalpy_e

h_ir_m_fit = h_iro2 - (2 * oxy_ref_h) - dh_iro2_exp  # fit to exp dH

dh_iro2 = h_iro2 - (2 * oxy_ref.enthalpy_e + h_ir_m_fit)

# +
e_iro3_pa = dft_energy_dict["iro3"]
e_iro3 = 4 * e_iro3_pa

PV_i = 0.
h_iro3 = e_iro3 + PV_i

dh_iro3 = h_iro3 - (3 * oxy_ref.enthalpy_e + h_ir_m_fit)

# +
h_iro3

print(4 * dft_energy_dict["iro3"])

# -6.442159*4

# +
print(dft_energy_dict["iro3"])
# -0.650702238408173

-6.46984746
# -

# # Calculate Gibbs Free Energy

# +
TdS_iro2 = TS_iro2_solid - TS_ir_metal - TS_o2_gas

# Alternate calculate of TdS for IrO2
# TdS_iro2 = dh_iro2_exp - dg_iro2_exp

dg_iro2 = dh_iro2 - TdS_iro2
# -

# # Calculating Gibbs free energy of IrO3 polymorph

# +
# Entropy for WOx from Barrin tables
# WO3: 75.898
# WO2: 50.543

factor = 75.898 / 50.543

# adjusted to reflect more O in IrO3 vs IrO2
TdS_iro3 = factor * TS_iro2_solid - TS_ir_metal - 3 / 2 * TS_o2_gas

# +
e_iro3_pa = dft_energy_dict["iro3"]
e_iro3 = 4 * e_iro3_pa

PV_i = 0.
h_iro3 = e_iro3 + PV_i

dh_iro3 = h_iro3 - (3 * oxy_ref.enthalpy_e + h_ir_m_fit)
dg_iro3 = dh_iro3 - TdS_iro3

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

# +
print(40 * "-")
print("O ref:", oxy_ref.enthalpy_e)
print("Ir ref:", h_ir_m_fit)

print(40 * "-")
print("TdS_iro2:", TdS_iro2)
print("TdS_iro3:", TdS_iro3)

print(40 * "-")
print("dh_iro2:", dh_iro2)
print("dh_iro3:", dh_iro3)
print("")
print("dg_iro2:", dg_iro2)
print("dg_iro3:", dg_iro3)
print("")
print("dh_iro2 - dh_iro3:", dh_iro2 - dh_iro3)
print("dg_iro2 - dg_iro3:", dg_iro2 - dg_iro3)

print(40 * "-")

# + active=""
#
#
#

# +
TdS_iro3

-1.8137510174978528 / 4
# -2.6028071136326894 / 4

# -2.6028071136326894 / 4

(4 * -0.650702238408173) - TdS_iro3

-0.471408668408173 * 4 - TdS_iro3
# -

# # Collecting variables to save

TdS_iro2
# TdS_iro3

# +
# assert False

# +
out_dict = dict(
    o_ref=oxy_ref.enthalpy_e,
    ir_ref=h_ir_m_fit,
    TdS_iro2=TdS_iro2,
    TdS_iro3=TdS_iro3,

    # "": ,
    )

# #######################################################################
import os; import json
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)
with open(os.path.join(directory, "data.json"), "w") as outfile:
    json.dump(out_dict, outfile, indent=2)
# #######################################################################

# +
# /mnt/c/Users/raulf2012/Dropbox/01_norskov/00_git_repos/PROJ_IrOx_Active_Learning_OER
# workflow/energy_treatment_deriv/calc_references

os.path.join(
    os.environ["PROJ_irox"],
    "workflow/energy_treatment_deriv/calc_references",
    "out_data/data.json"
    )
# -

# #######################################################################
import json
data_path = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/energy_treatment_deriv/calc_references",
    "out_data/data.json",
    )
with open(data_path, "r") as fle:
    data = json.load(fle)
# #######################################################################

data
