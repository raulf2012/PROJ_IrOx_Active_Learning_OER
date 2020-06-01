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

# +
# | - Import Modules
import os
import sys

import pickle
import yaml

# Decomposes ion into composition object that can handle the charge string
from pymatgen import MPRester
from pymatgen.core.ion import Ion
from pymatgen.analysis.pourbaix_diagram import PourbaixEntry, IonEntry
from pymatgen.entries.computed_entries import ComputedEntry

# #########################################################
# from dft_energies import ion_dict_solids_expt
#__|

# +
# | - Script Inputs
#This initializes the REST adaptor. Put your own API key in.

path_i = os.path.join(os.environ["PROJ_irox"], "config", "config.yml")
with open(path_i) as file:
    config_dict = yaml.load(file, Loader=yaml.FullLoader)

api_key = config_dict["materials_project"]["api_key"]

mpr = MPRester(api_key)  # Raul
#__|
# -

entries = mpr.get_pourbaix_entries(
    [
        "Ir",
        ]
    )

entry = entries[0]

entry.uncorrected_energy
# entry

# Pourbaix Entry : Ir1 O4 with energy = 7.4405, npH = -8.0, nPhi = -7.0, nH2O = 4.0, entry_id = None


assert False


# + jupyter={}
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


# + jupyter={}
# | - Main Code ****************************************************************

# | - Processing My Solid Entries
# for key in ion_dict_solids_expt:
pbx_solid_entries = []
for key, value in ion_dict_solids_expt.items():

    split_key = key.split("_")
    formula_i = split_key[0]

    # if len(split_key) > 1:
    #     # attribute_i = key.split("_")[-1]
    #     attribute_i = split_key[-1]
    # else:
    #     attribute_i = None

    comp = Ion.from_formula(formula_i)
    # comp = Ion.from_formula(key)
    energy = value
    # energy = ion_dict_solids_expt[key]  # + ion_correction * factor
    pbx_entry_ion = PourbaixEntry(
        ComputedEntry(
            comp,
            energy,
            # attribute={
            parameters={
                "full_name": key,
                },
            )
        )

    # AP pbx_entry_ion.name = key
    pbx_entry_ion.conc = 1
    pbx_solid_entries.append(pbx_entry_ion)

#__|


# | - Ion Entries


# Ion Reference Data

#Dictionary of reference state:experimental formation energy
ion_dict_Ir = mpr._make_request('/pourbaix_diagram/reference_data/Ir')
ion_dict_S = mpr._make_request('/pourbaix_diagram/reference_data/Sr')


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
    pbx_entry_ion = PourbaixEntry(IonEntry(comp, energy))
    # AP pbx_entry_ion.name = id['Name']
    pbx_entry_ion.conc = 1.0e-6
    if pbx_entry_ion.name not in ['jfdksl']:  # ['H2RuO2[2+]']: #["RuO4(aq)"]:
        pbx_ion_entries.append(pbx_entry_ion)
#__|


all_entries = pbx_solid_entries + pbx_ion_entries
#__|

# + active=""
#
#
#
#
#
#

# + jupyter={}
# def get_pourbaix_entries(self, chemsys):
#     """
#     A helper function to get all entries necessary to generate
#     a pourbaix diagram from the rest interface.

#     Args:
#         chemsys ([str]): A list of elements comprising the chemical
#             system, e.g. ['Li', 'Fe']
#     """
#     from pymatgen.analysis.pourbaix_diagram import PourbaixEntry, IonEntry
#     from pymatgen.analysis.phase_diagram import PhaseDiagram
#     from pymatgen.core.ion import Ion
#     from pymatgen.entries.compatibility import \
#         MaterialsProjectAqueousCompatibility

#     pbx_entries = []

#     # Get ion entries first, because certain ions have reference
#     # solids that aren't necessarily in the chemsys (Na2SO4)
#     url = '/pourbaix_diagram/reference_data/' + '-'.join(chemsys)
#     ion_data = self._make_request(url)
#     ion_ref_comps = [Composition(d['Reference Solid']) for d in ion_data]
#     ion_ref_elts = list(itertools.chain.from_iterable(
#         i.elements for i in ion_ref_comps))
#     ion_ref_entries = self.get_entries_in_chemsys(
#         list(set([str(e) for e in ion_ref_elts] + ['O', 'H'])),
#         property_data=['e_above_hull'], compatible_only=False)
#     compat = MaterialsProjectAqueousCompatibility("Advanced")
#     ion_ref_entries = compat.process_entries(ion_ref_entries)
#     ion_ref_pd = PhaseDiagram(ion_ref_entries)

#     # position the ion energies relative to most stable reference state
#     for n, i_d in enumerate(ion_data):
#         ion_entry = IonEntry(Ion.from_formula(i_d['Name']), i_d['Energy'])
#         refs = [e for e in ion_ref_entries
#                 if e.composition.reduced_formula == i_d['Reference Solid']]
#         if not refs:
#             raise ValueError("Reference solid not contained in entry list")
#         stable_ref = sorted(refs, key=lambda x: x.data['e_above_hull'])[0]
#         rf = stable_ref.composition.get_reduced_composition_and_factor()[1]
#         solid_diff = ion_ref_pd.get_form_energy(stable_ref) - i_d['Reference solid energy'] * rf
#         elt = i_d['Major_Elements'][0]
#         correction_factor = ion_entry.ion.composition[elt] / stable_ref.composition[elt]
#         ion_entry.energy += solid_diff * correction_factor
#         pbx_entries.append(PourbaixEntry(ion_entry, 'ion-{}'.format(n)))

#     # Construct the solid pourbaix entries from filtered ion_ref entries
#     extra_elts = set(ion_ref_elts) - {Element(s) for s in chemsys} \
#         - {Element('H'), Element('O')}
#     for entry in ion_ref_entries:
#         entry_elts = set(entry.composition.elements)
#         # Ensure no OH chemsys or extraneous elements from ion references
#         if not (entry_elts <= {Element('H'), Element('O')} or
#                 extra_elts.intersection(entry_elts)):
#             # replace energy with formation energy, use dict to
#             # avoid messing with the ion_ref_pd and to keep all old params
#             form_e = ion_ref_pd.get_form_energy(entry)
#             new_entry = deepcopy(entry)
#             new_entry.uncorrected_energy = form_e
#             new_entry.correction = 0.0
#             pbx_entry = PourbaixEntry(new_entry)
#             pbx_entries.append(pbx_entry)

#     return pbx_entries
