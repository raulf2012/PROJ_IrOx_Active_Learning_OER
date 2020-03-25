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
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Derivation of Corrections for Surface Energy Calculations

# +
import os
import sys
import sympy

sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "data"))

from proj_data_irox import *

E_h2o = h2o_ref
E_h2 = h2_ref

# + active=""
#
#
# -

# # General Gibbs Free Energy Formulation
# ---

# In general, the Gibbs free energy $G$ of a state is composed of an enthalpy and an entropy term
# $$
# G = H - TS = % TEMP
# \Bigg(
# \lbrack
# E_{elec} + E_{ZPE} + PV +
# \int_{0}^{T} C_p dT
# \rbrack
# - TS
# \Bigg)_i
# $$
#
# Sweeping every term except the DFT electronic energy into one term ($\varphi$) yields
# $$
# G_i = E_{elec,i} + \varphi _i
# $$

# +
E_elec = sympy.symbols("E_elec")
phi = sympy.symbols("phi")

G = E_elec + phi
G
# -

# Where ($\varphi$) is defined as
# $$
# \varphi =
# E_{ZPE} + PV +
# \int_{0}^{T} C_p dT
# - TS
# $$

# +
# Defining the Phi function
ZPE, PV, Cp, TS = sympy.symbols(["ZPE", "PV", "Cp", "TS"])

phi = ZPE + PV + Cp - TS

phi
# -

# ---

# # Adsorption Energy Formulas

# With this formulasim, we can now write out an adsorption reaction energy as follows
#
# $$
# \Delta G_{rxn} = G_{*O_xH_y} - G_{*} - x G_O - y G_H
# $$
#
# For convience, we then separate the electronic energy and corrections terms as before
#
#
# $$
# \Delta G_{rxn} =
# \Delta E_{rxn} + \Delta \varphi_{rxn} =
# % Electronic energy terms
# \lbrack
# E_{*O_xH_y} - E_{*} - x \cdot E_O - y \cdot E_H
# \rbrack +
# % Correction terms
# \lbrack
# \varphi_{*O_xH_y} - \varphi_{*} - x \cdot \varphi_O - y \cdot \varphi_H
# \rbrack
# $$

# + active=""
#
#
# -

# # Gas Phase References (O and H Reference Energies)
# ---

# The reference $H$ energy is simply obtained from a DFT calculation of an $H_2$ molecule
#
# <center> |-------- Hydrogen Reference Energy --------| </center>
#
# $$
# G^{ref}_{H} = \frac{1}{2} G_{H_2}
# $$
#
# Using this definition for the H reference energy we can show that the energy of an $H_2$ molecule is $0$
#
# $$2 \cdot H^{Ref} = H_{2}$$
#
# $$G^{f}_{H_2} = G_{H_2} - 2 \cdot G^{Ref}_{H} =
# G_{H_2} - 2 \cdot \frac{1}{2} G_{H_2} = 0
# $$
#
# ---

# Because of the well known self interaction errors associated with the $O_2$ molecule, we calculate the oxygen reference energy ($G^{ref}_{O}$) by invoking the chemical equation of water formation in terms of oxygen and hydrogen gas
#
# $$
# H_2O = \frac{1}{2} O_2 + H_2
# $$
#
# Rearranging for the oxygen term, we can express the energy of oxygen in this $H_2O$ and $H_2$ reference state as follows
# $$
# O^{Ref} = H_2O - H_2
# $$
#
# $$
# G^{Ref}_{O} = G_{H_2O} - G_{H_2} =
# \lbrack
# E_{H_2O} - E_{H_2}
# \rbrack
# +
# \lbrack
# \varphi_{H_2O} - \varphi_{H_2}
# \rbrack
# $$
#

# Again we can show that the energy of a water molecule in this scheme is 0
# $$
# G^{f}_{H_2O} = G_{H_2O} - G^{Ref}_{O} - 2 \cdot G^{Ref}_{H} =
# G_{H_2O} - (G_{H_2O} - G_{H_2}) - 2 \cdot (\frac{1}{2} G_{H_2}) = 0
# $$

# # Energy of O2 Molucule in H2O/H2 Referencec State

# Under the water and hydrogen as reference state, the energy of a oxygen molecule can be computed with the water formation equation
# $$
# H_2O = \frac{1}{2} O_2 + H_2
# $$
#
# With a reaction equation as follows
#
# $$
# \Delta G_{formation,H_2O} =
# G_{H_2O} - \frac{1}{2} G_{O_2} - G_{H_2} =
# \frac{4.92}{2}eV
# $$

# With this, we can rearrange for the energy of an $O_2$ molecule
#
# $$
# \frac{1}{2} G_{O_2} = -\Delta G_{formation,H_2O} + G_{H_2O} - G_{H_2}
# $$
#
# $$
# G_{O_2} =
# 2 \cdot
# \lbrack
# -\Delta G_{formation,H_2O} + G_{H_2O} - G_{H_2}
# \rbrack
# $$
#
# Instead of plugging DFT valuese into these terms, we simply use $H_2O$ and $H_2$ reference state which sets their energy to  0

# $$
# G_{O_2} =
# 2 \cdot
# \lbrack
# -\Delta G_{formation,H_2O}
# \rbrack =
# 4.92 eV
# $$
#

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
# #############################################################################
# #### H2O Phi term ####
phi_h2o = phi.subs({
    "Cp": cv_h2o,
    "PV": 0.,
    "TS": ts_h2o,
    "ZPE": zpe_h2o,
    })
print("Phi_H2O", phi_h2o)

G_H2O = G.subs({
    "phi": phi_h2o,
    "E_elec": E_h2o,
    })
print("G_H2O", G_H2O)

print("")
# #############################################################################
# #### H2 Phi term ####
phi_h2 = phi.subs({
    "Cp": cv_h2,
    "PV": 0.,
    "TS": ts_h2,
    "ZPE": zpe_h2,
    })
phi_h2
print("Phi_H2: ", phi_h2)
# print("Phi_H: ", phi_h2 / 2)

G_H2 = G.subs({
    "phi": phi_h2,
    "E_elec": E_h2,
    })
print("G_H2", G_H2)
# print("G_H", G_H2 / 2)

# +
phi_O_ref = phi_h2o - phi_h2
print("phi_O_ref: ", phi_O_ref)

E_O_ref = E_h2o - E_h2
print("E_O_ref: ", E_O_ref)


phi_H_ref = phi_h2 / 2
print("phi_H_ref: ", phi_H_ref)

E_H_ref = E_h2 / 2
print("E_H_ref: ", E_H_ref)

# + active=""
#
#
#
#
#
#
#
#
# -

# # Adsorbate Energy Corrections (*OOH, *OH, *O)

# +
zpe_oh_ads = ads_fe_dict["oh"]["zpe"]
cv_oh_ads = ads_fe_dict["oh"]["cv"]
ts_oh_ads = ads_fe_dict["oh"]["ts"]

zpe_o_ads = ads_fe_dict["o"]["zpe"]
cv_o_ads = ads_fe_dict["o"]["cv"]
ts_o_ads = ads_fe_dict["o"]["ts"]

zpe_ooh_ads = ads_fe_dict["ooh"]["zpe"]
cv_ooh_ads = ads_fe_dict["ooh"]["cv"]
ts_ooh_ads = ads_fe_dict["ooh"]["ts"]

# +
# #############################################################################
# #### *OH Phi term ####
phi_oh_ads = phi.subs({
    "Cp": cv_oh_ads,
    "PV": 0.,
    "TS": ts_oh_ads,
    "ZPE": zpe_oh_ads,
    })
print("Phi_*OH", phi_oh_ads)


# #############################################################################
# #### *O Phi term ####
phi_o_ads = phi.subs({
    "Cp": cv_o_ads,
    "PV": 0.,
    "TS": ts_o_ads,
    "ZPE": zpe_o_ads,
    })
print("Phi_*O", phi_o_ads)


# #############################################################################
# #### *OH Phi term ####
phi_ooh_ads = phi.subs({
    "Cp": cv_ooh_ads,
    "PV": 0.,
    "TS": ts_ooh_ads,
    "ZPE": zpe_ooh_ads,
    })
print("Phi_*OOH", phi_ooh_ads)

# +
D_phi_ooh_ads = phi_ooh_ads - (2 * phi_O_ref + 1 * phi_H_ref)
D_phi_oh_ads = phi_oh_ads - (1 * phi_O_ref + 1 * phi_H_ref)
D_phi_o_ads = phi_o_ads - (1 * phi_O_ref + 0 * phi_H_ref)

print("D_phi_ooh_ads: ", D_phi_ooh_ads)
print("D_phi_oh_ads: ", D_phi_oh_ads)
print("D_phi_o_ads: ", D_phi_o_ads)
# + active=""
#
#
#
#
# -

# # Collecting Important Terms for Saving

# +
out_dict = dict()

out_dict["D_phi_ooh_ads"] = float(D_phi_ooh_ads)
out_dict["D_phi_o_ads"] = float(D_phi_o_ads)
out_dict["D_phi_oh_ads"] = float(D_phi_oh_ads)

out_dict["phi_ooh_ads"] = float(phi_ooh_ads)
out_dict["phi_o_ads"] = float(phi_o_ads)
out_dict["phi_oh_ads"] = float(phi_oh_ads)

out_dict["phi_H_ref"] = float(phi_H_ref)
out_dict["phi_O_ref"] = float(phi_O_ref)

out_dict["E_H_ref"] = float(E_H_ref)
out_dict["E_O_ref"] = float(E_O_ref)


import pickle
with open("out_data/data.pickle", "wb") as fle:
    pickle.dump(out_dict, fle)
# -

# # Output from Michal's Script

# + active=""
# corr_OH: 0.30225
# corr_O: -0.0145
# corr_H:  0.31675
#
#
# OH: 11.14742354
# O:  7.44492759
# H:  3.70249595
