#!/usr/bin/env python

"""Bulk crystal Gibbs energy per formula unit

Author: Raul A. Flores
"""

kJmol = 96.485
kcalmol = 23.06035

#| - Default Energetics
ion_dict_solids_expt = {
    "Ir": 0.,
    # dG Thermochemical Data of Pure Substances good
    "IrO2": -1.952490,
    # calc by MB see IrOx_free_energy_calculator.py
    "IrO3": -1.9039513,
    "IrHO3": -3.1616628,
    "Ir2O7": -1.1782971,
    "IrO5H4": -5.998632,
    }
#__|

#| - Chris's Version
ion_dict_solids_expt = {
    "Ir": 0.,
    "IrO2": -1.952324,
    "IrO3": -1.903951,
    "IrHO3": -3.161663,
    "Ir2O7": -1.178304,
    "IrO6": -0.879916,
    "IrO5H4": -5.960167,
    }
#__|

#| - 181231 | NEW | an_energy_treatment.ipynb
# ion_dict_solids_expt = {
#     'Ir': 0.,
#     'IrO2': -1.95232354,
#     'IrO3': -2.17726054,
#     'IrHO3': -3.1616628,
#     }
#
# ion_dict_solids_expt = {
#     'Ir': 0.,
#     'IrO2': -1.9523235387946347,
#     'IrO3': -2.267671121530652}

ion_dict_solids_expt = {
    'Ir': 0.0,
    'IrO2': -1.9523235387946347,

    # 'IrO3_a-AlF3': -2.267671121530652,
    # 'IrO3_rutile-like': -2.13084673333495,
    # 'IrO3_battery': -1.95890971951947,
    }
#__|







#| - __old__
# #############################################################################
# Default Energetics
# ion_dict_solids_expt = {
#     "Ir": 0,
#
#     # dG Thermochemical Data of Pure Substances good
#     "IrO2": -188.386 / kJmol,
#
#     # calc by MB see IrOx_free_energy_calculator.py
#     "IrO3": -183.702743834 / kJmol,
#     "IrHO3": -305.053034883 / kJmol,
#     "Ir2O7": -113.688 / kJmol,
#     "IrO5H4": -578.778019727 / kJmol,
#     }

# #############################################################################
# Version 2
# ion_dict_solids_expt = {
#     "Ir": 0,
#     "IrO2": -188.3699405 / kJmol,
#     "IrO3": -155.342319414 / kJmol,
#     "IrHO3": -311.972187654 / kJmol,
#     }

#__|
