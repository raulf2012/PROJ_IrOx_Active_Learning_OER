# # !/usr/bin/env python

"""
"""


# | - Import Modules
import os
import sys

print("sys.executable:", sys.executable)

import pickle
import yaml

# Decomposes ion into composition object that can handle the charge string
from pymatgen import MPRester
from pymatgen.core.ion import Ion
from pymatgen.analysis.pourbaix_diagram import PourbaixEntry, IonEntry
from pymatgen.entries.computed_entries import ComputedEntry

# #########################################################
from dft_energies import ion_dict_solids_expt
#__|

# | - Script Inputs
# This initializes the REST adaptor. Put your own API key in.

path_i = os.path.join(os.environ["PROJ_irox"], "config", "config.yml")
with open(path_i) as file:
    config_dict = yaml.load(file, Loader=yaml.FullLoader)

api_key = config_dict["materials_project"]["api_key"]

mpr = MPRester(api_key)  # Raul

use_iro3 = True
# use_iro3 = False

stoich_i = "AB3"
#__|


def create_all_pourb_entries(
    stoich_i=None,

    ):
    if sys.argv[-1] == "AB3" or sys.argv[-1] == "AB2":
        stoich_i = sys.argv[-1]

    if stoich_i == "AB2":
        ion_dict_solids_expt.pop("IrO3_a-AlF3")
        ion_dict_solids_expt.pop("IrO3_rutile-like")
        ion_dict_solids_expt.pop("IrO3_battery")
        ion_dict_solids_expt.pop("IrO3_TEMP21")
    elif stoich_i == "AB3":
        pass

    # | - Methods
    #Used later to filter duplicate entries
    #If entry is already in entry_list, then return True
    def contains_entry(entry_list, entry):
        """
        """
        # | - contains_entry
        for e in entry_list:
            if e.entry_id == entry.entry_id or (abs(entry.energy_per_atom - e.energy_per_atom) < 1e-6 and entry.composition.reduced_formula == e.composition.reduced_formula):
                return True
        #__|

    # __|

    # | - Main Code ***************************************

    # | - Processing My Solid Entries
    # for key in ion_dict_solids_expt:
    pbx_solid_entries = []
    for key, value in ion_dict_solids_expt.items():

        split_key = key.split("_")
        formula_i = split_key[0]

        comp = Ion.from_formula(formula_i)
        energy = value
        pbx_entry_ion = PourbaixEntry(
            ComputedEntry(
                comp,
                energy,
                parameters={
                    "full_name": key,
                    },
                )
            )

        # AP pbx_entry_ion.name = key
        pbx_entry_ion.conc = 1
        pbx_solid_entries.append(pbx_entry_ion)

    # __|


    # | - Ion Entries

    #Dictionary of reference state:experimental formation energy
    ion_dict_Ir = mpr._make_request('/pourbaix_diagram/reference_data/Ir')

    for i, dummy in enumerate(ion_dict_Ir):
        print(ion_dict_Ir[i])
    print("")

    ion_dict = ion_dict_Ir


    # Calculate DFT reference energy for ions (See Persson et al, PRB (2012))
    pbx_ion_entries = []
    for id in ion_dict:
        # Ion name-> Ion comp name (ex. Fe[3+] -> Ion: Fe1 +3)
        comp = Ion.from_formula(id['Name'])

        energy = id['Energy']  # + ion_correction * factor
        print(id['Name'], comp, energy)
        # energy = -3.9840465215000007
        print(id['Name'], comp, energy)

        pbx_entry_ion = PourbaixEntry(IonEntry(comp, energy))
        # AP pbx_entry_ion.name = id['Name']
        pbx_entry_ion.conc = 1.0e-6
        if pbx_entry_ion.name not in ['jfdksl']:  # ['H2RuO2[2+]']: #["RuO4(aq)"]:
            pbx_ion_entries.append(pbx_entry_ion)
    #__|


    all_entries = pbx_solid_entries + pbx_ion_entries
    #__|

    # | - Save all_entries
    import os; import pickle
    directory = "out_data"
    if not os.path.exists(directory): os.makedirs(directory)
    path_i = os.path.join(directory, "all_entries_" + stoich_i + ".pickle")
    with open(path_i, "wb") as fle:
        pickle.dump(all_entries, fle)
    #__|


# print(20 * "# # ")
# print("All done!")
# assert False
