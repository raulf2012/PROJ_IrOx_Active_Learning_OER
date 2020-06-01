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

# + jupyter={}
# | - Import Modules
import os
print(os.getcwd())
import sys

import pickle
import yaml

from pymatgen import MPRester

# Decomposes ion into composition object that can handle the charge string
from pymatgen.core.ion import Ion

from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.analysis.pourbaix_diagram import (
    PourbaixEntry,
    IonEntry,
    )

#GGA/GGA+U Mixing Scheme
from pymatgen.entries.compatibility import MaterialsProjectAqueousCompatibility
from pymatgen.entries.computed_entries import ComputedEntry

import warnings
warnings.filterwarnings('ignore')

# #############################################################################
# #############################################################################
from dft_energies import ion_dict_solids_expt
#__|
# -

ion_dict_solids_expt

# +
# ion_dict_solids_expt = ion_dict_solids_expt["IrO2"]

# ion_dict_solids_expt = dict(IrO2=-1.9523235387946347)

# +
MU_H2O = -2.4583
# MU_H2O = -2.4577835

# IrO2
IrO2_G = -1.9523235387946347 - 2 * MU_H2O

# IrO3
# IrO3_G = -2.267671121530652 - 3 * MU_H2O
IrO3_G = -2.067671121530652 - 3 * MU_H2O


(IrO3_G - IrO2_G) / 2  # 1.0714762086319913
# -

1.0714762086319913 - 1.1214762086319912

# +
path_i = os.path.join(os.environ["PROJ_irox"], "config", "config.yml")
with open(path_i) as file:
    config_dict = yaml.load(file, Loader=yaml.FullLoader)

api_key = config_dict["materials_project"]["api_key"]

# +
# assert False

# +
# | - Script Inputs
#This initializes the REST adaptor. Put your own API key in.
mpr = MPRester(api_key)  # Raul

#__|

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

#__|

# | - Main Code ****************************************************************
#Entries are the basic unit for thermodynamic and other analyses in pymatgen.
entries = mpr.get_entries_in_chemsys(['O', 'H'])

# | - Ion Reference Data
#Dictionary of reference state:experimental formation energy

# ion_dict_Co = mpr._make_request('/pourbaix_diagram/reference_data/Ru')
# ion_dict_S = mpr._make_request('/pourbaix_diagram/reference_data/Y')

ion_dict_Co = mpr._make_request('/pourbaix_diagram/reference_data/Ir')
ion_dict_S = mpr._make_request('/pourbaix_diagram/reference_data/Sr')


for i, dummy in enumerate(ion_dict_Co):
    print(ion_dict_Co[i])
    print('hi')

# for i,dummy in enumerate(ion_dict_S):
# 	print ion_dict_S[i]['Name'], ion_dict_S[i]['Energy']

ion_dict = ion_dict_Co  # + ion_dict_S
#ion_dict = ion_dict_Co + ion_dict_S
#__|

#NOTE This line assumes that the every entry in the experimental ion energy
# has the same ref. st. solid
ref_state = str(ion_dict[0]['Reference Solid'])
ref_dict = {ref_state: ion_dict[0]['Reference solid energy']}

# Run aqueouscorrection on the entries
# Entries without applicable corrections will be discarded
# Implements the GGA/GGA+U mixing scheme
aqcompat = MaterialsProjectAqueousCompatibility()

entries_aqcorr = list()
for entry in entries:
    # Applies corrections to entry, if none applicable it gets rid of entry
    aq_corrected_entry = aqcompat.process_entry(entry)
    # If entry already in entries_aqcorr then don't add to list
    if not contains_entry(entries_aqcorr, aq_corrected_entry):
        entries_aqcorr.append(aq_corrected_entry)

# Generate a phase diagram to consider only solid entries stable in water.
pd = PhaseDiagram(entries_aqcorr)
stable_solids = pd.stable_entries
stable_solids_minus_h2o = [entry for entry in stable_solids if
    entry.composition.reduced_formula not in ["H2", "O2", "H2O", "H2O2"]]



# +
pbx_solid_entries = []
for entry in stable_solids_minus_h2o:
    pbx_entry = PourbaixEntry(entry)

    # Replace E with newly corrected E
    pbx_entry.g0_replace(pd.get_form_energy(entry))
    pbx_entry.reduced_entry()  # Applies reduction factor?????
    pbx_solid_entries.append(pbx_entry)

# | - Processing My Solid Entries
for key, value in ion_dict_solids_expt.items():
    print("IDJFIJDS")
    split_key = key.split("_")
    formula_i = split_key[0]

    comp = Ion.from_formula(formula_i)
    energy = value
    pbx_entry_ion = PourbaixEntry(
        ComputedEntry(
            comp,
            energy,
            attribute={
                "full_name": key}))

    # AP pbx_entry_ion.name = key
    pbx_entry_ion.conc = 1
    pbx_solid_entries.append(pbx_entry_ion)

for a in pbx_solid_entries:
    print(a, a.conc)
    print('hkjo')

#__|

# +
self = pbx_entry_ion


self.uncorrected_energy

self.conc_term

self.nH2O

# +
-2.369

-2.4577835

-2.46  # 1.23 * 2

# +
ion_dict_solids_expt

2.9643 - -1.9523235387946347
# -

-2.4583

1.23 * 2

# +
# PourbaixEntry?
# -

pbx_entry_ion.energy

assert False

# + jupyter={}
# | - Ion Entries
# Calculate DFT reference energy for ions (See Persson et al, PRB (2012))
pbx_ion_entries = []
for id in ion_dict:
    # Ion name-> Ion comp name (ex. Fe[3+] -> Ion: Fe1 +3)
    comp = Ion.from_formula(id['Name'])
    energy = id['Energy']  # + ion_correction * factor
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

# Pickling data ######################################################
# import os; import pickle
# directory = "out_data"
# if not os.path.exists(directory): os.makedirs(directory)
# with open(os.path.join(directory, "all_entries.pickle"), "wb") as fle:
#     pickle.dump(all_entries, fle)
# #####################################################################

# + active=""
#
#
#
#

# + jupyter={}
# from pymatgen.analysis.pourbaix_diagram import (
#     PourbaixDiagram,
#     PourbaixPlotter,
#     generate_entry_label,
#     )


# # #############################################################################
# import pickle; import os
# path_i = os.path.join(
#     "out_data",
#     "all_entries.pickle")
# with open(path_i, "rb") as fle:
#     all_entries = pickle.load(fle)
# # #############################################################################

# all_entries = all_entries[0:2]

# all_entries


# pourbaix = PourbaixDiagram(all_entries)  # comp_dict={'Ru':0.5,'Y':0.5})
