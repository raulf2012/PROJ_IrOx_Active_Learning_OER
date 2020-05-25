# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_irox] *
#     language: python
#     name: conda-env-PROJ_irox-py
# ---

# + Collapsed="false" jupyter={}
import sys
import os
print(os.getcwd())

sys.path.insert(
    0,
    os.path.join(
        os.environ["PROJ_irox"],
        "data"))
from proj_data_irox import proj_dir_name, irox_bulk_color_map

# #############################################################################
import copy

import pickle

import numpy as np

import chart_studio.plotly as py
import plotly.graph_objs as go

import plotly.offline as py_off
from plotly.offline import (
    init_notebook_mode,
    iplot,
    )

from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixPlotter

# #############################################################################
from methods import (
    create_pourbaix_plot,
    create_outside_borders,
    create_pourb_entry_outline,
    create_oer_equil_line,

    get_base_spec,
    get_spec_entries,
    )

init_notebook_mode(connected=True)


# #############################################################################
from plotting.my_plotly import (
    add_minor_ticks,
    my_plotly_plot,
    add_duplicate_axes)

# + Collapsed="false"
# # %%capture

if True:
    !python sc_create_all_entries.py

# + Collapsed="false"
# #############################################################################
path_i = os.path.join(
    "out_data",
    "all_entries.pickle")
with open(path_i, "rb") as fle:
    all_entries = pickle.load(fle)
# #############################################################################
# -

ir_entry = get_base_spec("Ir", all_entries)
iro2_entry = get_base_spec("IrO2", all_entries)
ir_ion_entry = get_base_spec("IrO4-", all_entries)


# # Find transition method

def find_pour_trans(
    PourbaixDiagram=None,
    range=[0, 3],
    num=500,
    ):
    """
    """
    PD = PourbaixDiagram

    transition_V = None

    entries_sweep = []
    for i_cnt, V_i in enumerate(np.linspace(range[0], range[1], num=num)):
        stable_entry = PD.get_stable_entry(0, V_i)
        entry_name = stable_entry.name

        if i_cnt != 0:
            prev_entry = entries_sweep[i_cnt - 1]

            if prev_entry != entry_name:
                # print("V_i:", V_i)
                transition_V = V_i

                break

        entries_sweep.append(entry_name)

    return(transition_V)


# # Ir --> IrO2 Transition

# +
PD = PourbaixDiagram([
    ir_entry,
    iro2_entry,
    ])

# find_pour_trans(PourbaixDiagram=PD)
ir_iro2_trans = find_pour_trans(PourbaixDiagram=PD, range=[0.5, 0.9], num=1000)
print("ir_iro2_trans:", ir_iro2_trans)
# -

# # IrO2 --> IrO3 transitions
# ---

# # IrO2 --> a-IrO3 Transition

# +
# # %%capture

out_dict = get_spec_entries(
    ["IrO3_a-AlF3"],
    all_entries)
a_iro3_entry = out_dict["IrO3_a-AlF3"]

PD = PourbaixDiagram([
    iro2_entry,
    a_iro3_entry,
    ])

iro2_a_iro3_trans = find_pour_trans(PourbaixDiagram=PD, range=[1, 1.5], num=1000)
print("iro2_a_iro3_trans:", iro2_a_iro3_trans)

# +
# assert False
# -

# # IrO2 --> rutile-IrO3 Transition

# +
# # %%capture

out_dict = get_spec_entries(
    ["IrO3_rutile-like"],
    all_entries)
r_iro3_entry = out_dict["IrO3_rutile-like"]

PD = PourbaixDiagram([
    iro2_entry,
    r_iro3_entry,
    ])

iro2_r_iro3_trans = find_pour_trans(PourbaixDiagram=PD, range=[1, 1.5], num=1000)
print("iro2_r_iro3_trans:", iro2_r_iro3_trans)
# -

# # IrO2 --> b-IrO3 Transition

# +
# # %%capture

out_dict = get_spec_entries(
    ["IrO3_battery"],
    all_entries)
b_iro3_entry = out_dict["IrO3_battery"]

PD = PourbaixDiagram([
    iro2_entry,
    b_iro3_entry,
    ])

iro2_b_iro3_trans = find_pour_trans(PourbaixDiagram=PD, range=[1, 1.5], num=1000)
print("b_iro2_a_iro3_trans:", iro2_b_iro3_trans)
# -

# # IrO3 --> Ir Ion transitions
# ---

# # a-IrO3 --> Ir[4+] Transition

# +
PD = PourbaixDiagram([
    a_iro3_entry,
    ir_ion_entry,
    ])

a_iro3_ir_ion_trans = find_pour_trans(PourbaixDiagram=PD, range=[1.5, 2], num=1000)
print("a_iro3_ir_ion_trans:", a_iro3_ir_ion_trans)
# -

# # r-IrO3 --> Ir[4+] Transition

# +
PD = PourbaixDiagram([
    r_iro3_entry,
    ir_ion_entry,
    ])

r_iro3_ir_ion_trans = find_pour_trans(PourbaixDiagram=PD, range=[1.5, 2], num=1000)
print("r_iro3_ir_ion_trans:", r_iro3_ir_ion_trans)
# -

# # b-IrO3 --> Ir[4+] Transition

# +
PD = PourbaixDiagram([
    b_iro3_entry,
    ir_ion_entry,
    ])

b_iro3_ir_ion_trans = find_pour_trans(PourbaixDiagram=PD, range=[1.5, 2], num=1000)
print("b_iro3_ir_ion_trans:", b_iro3_ir_ion_trans)

# + active=""
#
#
#
#
#

# +
out_dict = dict(
    ir_iro2_trans=ir_iro2_trans,

    iro2_a_iro3_trans=iro2_a_iro3_trans,
    iro2_r_iro3_trans=iro2_r_iro3_trans,
    iro2_b_iro3_trans=iro2_b_iro3_trans,

    a_iro3_ir_ion_trans=a_iro3_ir_ion_trans,
    r_iro3_ir_ion_trans=r_iro3_ir_ion_trans,
    b_iro3_ir_ion_trans=b_iro3_ir_ion_trans,
    )



# Pickling data ######################################################
import os; import pickle
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)
with open(os.path.join(directory, "bulk_pourb_transitions.pickle"), "wb") as fle:
    pickle.dump(out_dict, fle)
# #####################################################################

# +
# # IrO2(s)
# # Ir(s)

# entries_sweep = []
# for i_cnt, V_i in enumerate(np.linspace(0, 1, num=1000)):
#     stable_entry = PD.get_stable_entry(0, V_i)
#     entry_name = stable_entry.name

#     if i_cnt != 0:
#         prev_entry = entries_sweep[i_cnt - 1]

#         if prev_entry != entry_name:
#             print("V_i:", V_i)

#     entries_sweep.append(entry_name)

# # PD.get_stable_entry(0, i).name
# # 0.7410741074107411
