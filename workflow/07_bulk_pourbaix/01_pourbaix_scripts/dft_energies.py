#!/usr/bin/env python

"""Bulk crystal Gibbs energy per formula unit

Obtained from following dir:
$PROJ_irox/workflow/energy_treatment_deriv/calc_references

Author: Raul A. Flores
"""

# | - Import Modules
import os
import sys

import json

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import calc_dH
#__|


# #######################################################################
energy_treatment_data_path = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/energy_treatment_deriv/calc_references",
    "out_data/data.json")
with open(energy_treatment_data_path, "r") as fle:
    energy_treatment_data_dict = json.load(fle)

TdS_iro3 = energy_treatment_data_dict["TdS_iro3"]
TdS_irho3 = energy_treatment_data_dict["TdS_irho3"]
# #######################################################################


kJmol = 96.485
kcalmol = 23.06035


# TdS_irho3 = -0.972
# TdS_irho3 = -1.263

# Newest effort, using energetics at 500eV (AL energies)
ion_dict_solids_expt = {

    # #########################################################################
    # Uncomment to get Hydroxide phase in diagram
    # "IrHO3_a-AlF3": (5 * calc_dH(-6.151177564999999, stoich="IrHO3")) - TdS_irho3,

    # Raw DFT values:
    # See following WF location:
    #   https://beta.workflowy.com/#/de6d404eb1f5

    # -6.125331012  # 03
    # -6.134396575  # 11
    # -6.156622886  # 20

    # -6.151177564999999  # a-IrO3 hydrated (9lmkmh8s8r reduced cell)


    # #########################################################################
    "Ir": 0.0,

    # #########################################################################
    "IrO2": -1.9524900243561174,  # Completely exp. (Barrin)


    # #########################################################################
    # "IrO3_a-AlF3": (4 * -0.650702238408173) - TdS_iro3,
    "IrO3_a-AlF3": (4 * calc_dH(-6.46984746, stoich="AB3")) - TdS_iro3,

    # #########################################################################
    "IrO3_rutile-like": (4 * calc_dH(-6.456962203125, stoich="AB3")) - TdS_iro3,  # b5cgvsb16w

    # #########################################################################
    "IrO3_battery": (4 * calc_dH(-6.41597669640625, stoich="AB3")) - TdS_iro3,  # Current calc, running


    # #########################################################################
    # SrIro3 type with Sr removed
    # "IrO3_cubic_perovskite": -1.0965785774978554,  # vwxfn3blxi
    # "IrO3_cubic_perovskite": -1.0965785774978554 - 0.3,  # TEST

    # #########################################################################
    # IrO3 polymorphs #########################################################

    #  "IrO3_TEMP00": (4 * -0.650702238) - TdS_iro3,  # 8p8evt9pcg
    #  "IrO3_TEMP01": (4 * -0.648304405) - TdS_iro3,  # zimixdvdxd
    #  "IrO3_TEMP02": (4 * -0.64247475)  - TdS_iro3,  # xw9y6rbkxr
    #  "IrO3_TEMP03": (4 * -0.637816982) - TdS_iro3,  # b5cgvsb16w
    #  "IrO3_TEMP04": (4 * -0.634008712) - TdS_iro3,  # mj7wbfb5nt
    #  "IrO3_TEMP05": (4 * -0.619444569) - TdS_iro3,  # 949rnem5z2
    #  "IrO3_TEMP06": (4 * -0.616537586) - TdS_iro3,  # 6pvt6r95ve
    #  "IrO3_TEMP07": (4 * -0.608641347) - TdS_iro3,  # 8l919k6s7p
    #  "IrO3_TEMP08": (4 * -0.60295535)  - TdS_iro3,  # 7ic1vt7pz4
    #  "IrO3_TEMP09": (4 * -0.601362395) - TdS_iro3,  # zwvqnhbk7f
    #  "IrO3_TEMP10": (4 * -0.599971464) - TdS_iro3,  # 8wxibl7lm4
    #  "IrO3_TEMP11": (4 * -0.599769987) - TdS_iro3,  # ngn4xec1mo
    #  "IrO3_TEMP12": (4 * -0.599425003) - TdS_iro3,  # xlziv2zr9g
    #  "IrO3_TEMP13": (4 * -0.598528052) - TdS_iro3,  # mp6lno9jzr
    #  "IrO3_TEMP14": (4 * -0.598300114) - TdS_iro3,  # zgxg9o7kny
    #  "IrO3_TEMP15": (4 * -0.580270245) - TdS_iro3,  # v2blxebixh
    #  "IrO3_TEMP16": (4 * -0.578769959) - TdS_iro3,  # 6384vt7pml
    #  "IrO3_TEMP17": (4 * -0.57347596)  - TdS_iro3,  # 9i6ixublcr
    #  "IrO3_TEMP18": (4 * -0.565705101) - TdS_iro3,  # xlbfb49wml
    #  "IrO3_TEMP19": (4 * -0.563072263) - TdS_iro3,  # cqbrnhbacg
    #  "IrO3_TEMP20": (4 * -0.56131933)  - TdS_iro3,  # bgcpc2vabf
    "IrO3_TEMP21": (4 * -0.55375463)  - TdS_iro3,  # nrml6dms9l

    }


#| - __old__

# calc_dH(-6.389916089375, stoich="AB3")


# | - Default Energetics
# ion_dict_solids_expt = {
#     "Ir": 0.,
#     # dG Thermochemical Data of Pure Substances good
#     "IrO2": -1.952490,
#     # calc by MB see IrOx_free_energy_calculator.py
#     "IrO3": -1.9039513,
#     "IrHO3": -3.1616628,
#     "Ir2O7": -1.1782971,
#     "IrO5H4": -5.998632,
#     }
#__|

# | - Chris's Version
# ion_dict_solids_expt = {
#     "Ir": 0.,
#     "IrO2": -1.952324,
#     "IrO3": -1.903951,
#     "IrHO3": -3.161663,
#     "Ir2O7": -1.178304,
#     "IrO6": -0.879916,
#     "IrO5H4": -5.960167,
#     }
#__|


#  tmp = calc_dH(-6.456962203125, stoich="AB3")
#  print("ISFIIDIFJDISJFIDSF")
#  print("R-IrO3", tmp)
#  print("ISFIIDIFJDISJFIDSF")
#  print("")
#  tmp = calc_dH(-6.41597669640625, stoich="AB3")
#  print("ISFIIDIFJDISJFIDSF")
#  print("B-IrO3", tmp)
#  print("ISFIIDIFJDISJFIDSF")


#__|
