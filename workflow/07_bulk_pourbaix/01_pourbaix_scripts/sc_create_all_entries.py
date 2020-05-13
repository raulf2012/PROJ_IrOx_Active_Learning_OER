#!/usr/bin/env python

"""
"""

# | - Import Modules
import pickle

# Decomposes ion into composition object that can handle the charge string
from pymatgen import MPRester
from pymatgen.core.ion import Ion
from pymatgen.analysis.pourbaix_diagram import PourbaixEntry, IonEntry
from pymatgen.entries.computed_entries import ComputedEntry

# #########################################################
from dft_energies import ion_dict_solids_expt
#__|

# | - Script Inputs
#This initializes the REST adaptor. Put your own API key in.
mpr = MPRester('NJTXWGbreuLAq8O5')  # Raul
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
#  ion_dict_S = mpr._make_request('/pourbaix_diagram/reference_data/Sr')


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
with open(os.path.join(directory, "all_entries.pickle"), "wb") as fle:
    pickle.dump(all_entries, fle)
#__|





# | - __old__

#  from pymatgen.analysis.phase_diagram import PhaseDiagram
#  import warnings
#  warnings.filterwarnings('ignore')
#  #GGA/GGA+U Mixing Scheme
#  from pymatgen.entries.compatibility import MaterialsProjectAqueousCompatibility

#Entries are the basic unit for thermodynamic and other analyses in pymatgen.
#  entries = mpr.get_entries_in_chemsys(['O', 'H'])



#| - __old__

#NOTE This line assumes that the every entry in the experimental ion energy
# has the same ref. st. solid
#  ref_state = str(ion_dict[0]['Reference Solid'])
#  ref_dict = {ref_state: ion_dict[0]['Reference solid energy']}

# Run aqueouscorrection on the entries
# Entries without applicable corrections will be discarded
# Implements the GGA/GGA+U mixing scheme

#  aqcompat = MaterialsProjectAqueousCompatibility()
#
#  entries_aqcorr = list()
#  for entry in entries:
#      # Applies corrections to entry, if none applicable it gets rid of entry
#      aq_corrected_entry = aqcompat.process_entry(entry)
#      # If entry already in entries_aqcorr then don't add to list
#      if not contains_entry(entries_aqcorr, aq_corrected_entry):
#          entries_aqcorr.append(aq_corrected_entry)

# Generate a phase diagram to consider only solid entries stable in water.
#  pd = PhaseDiagram(entries_aqcorr)
#  stable_solids = pd.stable_entries
#  stable_solids_minus_h2o = [entry for entry in stable_solids if
#      entry.composition.reduced_formula not in ["H2", "O2", "H2O", "H2O2"]]


#  for entry in stable_solids_minus_h2o:
#      pbx_entry = PourbaixEntry(entry)
#
#      # Replace E with newly corrected E
#      pbx_entry.g0_replace(pd.get_form_energy(entry))
#      pbx_entry.reduced_entry()  # Applies reduction factor?????
#      pbx_solid_entries.append(pbx_entry)
#__|




# print("")
# print("entry attributes:")
# for i in all_entries:
#     if hasattr(i.entry, 'attribute'):
#         print(i.entry.attribute)
#

# | - __old__
# 'Ru': 0,
# 'RuO2': 1.0*(-252.681)/kJmol,
# 'RuO4H4': -691.0/kJmol, # AP - from expt. paper
# 'Ru2O7': -147.473202819/kJmol,
# 'RuO3': -235.183254948/kJmol,
# #'RuO4': -293.694662969/kJmol,
# #'H8RuO6': -1296.50182566/kJmol, ## AP - completely disconnected, from materials project for Pt
# 'H2Ru2O6': -795.867243923/kJmol,
# 'RuHO3': -418.32273716/kJmol,
# 'RuO3H3': -5.480340000000, # from: file:///Users/anjli/Downloads/229_1535_3_PB%20(2).pdf; https://hrcak.srce.hr/file/230620
# 'RuO4H5': -766.0/kJmol, # from: file:///Users/anjli/Downloads/229_1535_3_PB%20(2).pdf; https://hrcak.srce.hr/file/230620
# 'Y': 0,
# 'Y2Ru2O7': -2367.44953534/kJmol,
# 'HRu4O12': -1242.80340535/kJmol,
# 'HRu2O6': -714.064039347/kJmol,
# 'H3Ru4O12': -1544.3919267/kJmol
#__|

# replot(all_entries)


# #,comp_dict={'Ir':0.5}):
# def replot(all_entries):
#     """
#     """
#     # | - replot
#     pourbaix = PourbaixDiagram(all_entries)  # comp_dict={'Ru':0.5,'Y':0.5})
#     plotter = PourbaixPlotter(pourbaix)
#     print_entries(all_entries)
#     plt = plotter.get_pourbaix_plot(limits=[[0, 14], [0.0, 2.0]])
#     f = plt.gcf()
#     f.set_size_inches((21.5, 21.5))
#     plt.savefig(plotname + '.pdf')
#     plt.savefig(plotname + '.png')
#
#     # plt2=plotter.get_pourbaix_plot_colorfill_by_element(
#     #     limits=[[-2,16],[-3,3]])
#     #plt2.savefig(plotname+'colorfill.pdf')
#     #plt2.savefig(plotname+'colorfill.png')
#
#     plotter.show()
#     #__|
#
# def print_entries(all_entries):
#     """
#     """
#     # | - print_entries
#     for e in all_entries:
#         print(e)
#     #__|

#__|
