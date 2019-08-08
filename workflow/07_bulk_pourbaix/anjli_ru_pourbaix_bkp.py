##!/usr/bin/env python3


#| - Import Modules
import sys, os

import pymatgen
print(pymatgen.__file__)
from pymatgen import MPRester
a = MPRester("8EKmva1TerBtUdxz") #michal

from pymatgen.core.ion import Ion #Decomposes ion into composition object that can handle the charge string
from pymatgen import Element #Accesses properties of element inputed (as a string) using an enumeration
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.analysis.pourbaix_diagram import PourbaixEntry, IonEntry, MultiEntry
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram
from pymatgen.analysis.pourbaix_diagram import PourbaixPlotter
from pymatgen.entries.compatibility import MaterialsProjectAqueousCompatibility #GGA/GGA+U Mixing Scheme
from pymatgen.entries.computed_entries import ComputedEntry
import sys
import warnings
warnings.filterwarnings('ignore')
#__|

#| - Script Inputs
plotname='Volume_pourbaix_RuOx_only_pymatgen'

#| - Energetics Data
kJmol = 96.485
kcalmol = 23.06035
ion_dict_solids_expt={
'Ru': 0,
'RuO2': 1.0*(-252.681)/kJmol,
'RuO4H4': -691.0/kJmol, # AP - from expt. paper
'Ru2O7': -147.473202819/kJmol,
'RuO3': -235.183254948/kJmol,
#'RuO4': -293.694662969/kJmol,
#'H8RuO6': -1296.50182566/kJmol, ## AP - completely disconnected, from materials project for Pt
'H2Ru2O6': -795.867243923/kJmol,
'RuHO3': -418.32273716/kJmol,
'RuO3H3': -5.480340000000, # from: file:///Users/anjli/Downloads/229_1535_3_PB%20(2).pdf; https://hrcak.srce.hr/file/230620
'RuO4H5': -766.0/kJmol, # from: file:///Users/anjli/Downloads/229_1535_3_PB%20(2).pdf; https://hrcak.srce.hr/file/230620
'Y': 0,
'Y2Ru2O7': -2367.44953534/kJmol,
'HRu4O12': -1242.80340535/kJmol,
'HRu2O6': -714.064039347/kJmol,
'H3Ru4O12': -1544.3919267/kJmol
}

'''
'Ir': 0,
'IrO2' : -188.386/kJmol,#dG Thermochemical Data of Pure Substances good
'IrO3' :   -183.702743834/kJmol , # calc by MB  see IrOx_free_energy_calculator.py
'IrHO3' : -305.053034883/kJmol , # calc by MB see IrOx_free_energy_calculator.py
'Ir2O7' :   -113.688/kJmol , # calc by MB  see IrOx_free_energy_calculator.py
'IrO5H4' :   -578.778019727/kJmol , # calc by MB  see IrOx_free_energy_calculator.py
'''
#__|

#__|

#| - Methods
#Used later to filter duplicate entries
#If entry is already in entry_list, then return True
def contains_entry(entry_list, entry):
    """
    """
    #| - contains_entry
    for e in entry_list:
        if e.entry_id == entry.entry_id or (abs(entry.energy_per_atom - e.energy_per_atom) < 1e-6 and entry.composition.reduced_formula == e.composition.reduced_formula):
            return True
    #__|

def replot(all_entries): #,comp_dict={'Ir':0.5}):
    """
    """
    #| - replot
    pourbaix = PourbaixDiagram(all_entries) #, comp_dict={'Ru':0.5,'Y':0.5})
    plotter = PourbaixPlotter(pourbaix)
    print_entries(all_entries)
    plt=plotter.get_pourbaix_plot(limits=[[0, 14],[0.0, 2.0]])
    f = plt.gcf()
    f.set_size_inches((21.5, 21.5))
    plt.savefig(plotname+'.pdf')
    plt.savefig(plotname+'.png')
    #plt2=plotter.get_pourbaix_plot_colorfill_by_element(limits=[[-2, 16],[-3, 3]])
    #plt2.savefig(plotname+'colorfill.pdf')
    #plt2.savefig(plotname+'colorfill.png')
    plotter.show()
    #__|

def print_entries(all_entries):
    """
    """
    #| - print_entries
    for e in all_entries:
        print(e)
    #__|

#__|


#| - Main Code ****************************************************************
#This initializes the REST adaptor. Put your own API key in.
mpr = MPRester('8EKmva1TerBtUdxz')

#Entries are the basic unit for thermodynamic and other analyses in pymatgen.
#entries = mpr.get_entries_in_chemsys(['Co','S', 'O', 'H'])
#entries = mpr.get_entries_in_chemsys(['S', 'O', 'H'])
entries = mpr.get_entries_in_chemsys(['O', 'H'])
#print entries

#Dictionary of reference state:experimental formation energy
ion_dict_Co = mpr._make_request('/pourbaix_diagram/reference_data/Ru')
ion_dict_S = mpr._make_request('/pourbaix_diagram/reference_data/Y')

for i,dummy in enumerate(ion_dict_Co):
    print(ion_dict_Co[i])
    print('hi')

'''
for i,dummy in enumerate(ion_dict_S):
	print ion_dict_S[i]['Name'], ion_dict_S[i]['Energy']
'''

ion_dict = ion_dict_Co + ion_dict_S

#N# Note that this line assumes that the every entry in the experimental ion energy has the same ref. st. solid
ref_state=str(ion_dict[0]['Reference Solid'])
ref_dict = {ref_state: ion_dict[0]['Reference solid energy']}

# Run aqueouscorrection on the entries. Entries without applicable corrections will be discarded
aqcompat = MaterialsProjectAqueousCompatibility() #Implements the GGA/GGA+U mixing scheme,

entries_aqcorr = list()

for entry in entries:
    aq_corrected_entry = aqcompat.process_entry(entry) #Applies corrections to entry, if none applicable it gets rid of entry
    if not contains_entry(entries_aqcorr, aq_corrected_entry): #If entry already in entries_aqcorr then don't add to list
        entries_aqcorr.append(aq_corrected_entry)


### Generate a phase diagram to consider only solid entries stable in water.
pd = PhaseDiagram(entries_aqcorr)
stable_solids = pd.stable_entries
stable_solids_minus_h2o = [entry for entry in stable_solids if
			   entry.composition.reduced_formula not in ["H2", "O2", "H2O", "H2O2"]]
pbx_solid_entries = []


for entry in stable_solids_minus_h2o:
    pbx_entry = PourbaixEntry(entry)
    pbx_entry.g0_replace(pd.get_form_energy(entry)) #Replace E with newly corrected E
    pbx_entry.reduced_entry() #Applies reduction factor?????
    pbx_solid_entries.append(pbx_entry)


for key in ion_dict_solids_expt:
	comp = Ion.from_formula(key)
	energy = ion_dict_solids_expt[key] #+ ion_correction * factor
	pbx_entry_ion = PourbaixEntry(ComputedEntry(comp, energy))
	### AP pbx_entry_ion.name = key
	pbx_entry_ion.conc = 1
	pbx_solid_entries.append(pbx_entry_ion)


for a in  pbx_solid_entries:
    print(a, a.conc)
    print('hkjo')


### Calculate DFT reference energy for ions (See Persson et al, PRB (2012))

pbx_ion_entries = []

for id in ion_dict:
    comp = Ion.from_formula(id['Name']) # Ion name-> Ion comp name (ex. Fe[3+] -> Ion: Fe1 +3)
    energy = id['Energy'] #+ ion_correction * factor
    print(id['Name'], comp, energy)
    pbx_entry_ion = PourbaixEntry(IonEntry(comp, energy))
    ### AP pbx_entry_ion.name = id['Name']
    pbx_entry_ion.conc = 1.0e-6
    if pbx_entry_ion.name not in ['jfdksl']: # ['H2RuO2[2+]']: #["RuO4(aq)"]:
    	pbx_ion_entries.append(pbx_entry_ion)

all_entries = pbx_solid_entries + pbx_ion_entries


replot(all_entries)

#__|
