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


# Newest effort, using energetics at 500eV (AL energies)
ion_dict_solids_expt = {
    'Ir': 0.0,
    'IrO2': -1.9524900243561174,  # Completely exp. (Barrin)

    'IrO3_a-AlF3': -1.8137510174978528,  # Real
    # 'IrO3_a-AlF3': -2.2,  # TEST

    "IrHO3": -3.161663,
    # "Ir2O7": -1.178304,
    # "IrO6": -0.879916,
    # "IrO5H4": -5.960167,

    }
