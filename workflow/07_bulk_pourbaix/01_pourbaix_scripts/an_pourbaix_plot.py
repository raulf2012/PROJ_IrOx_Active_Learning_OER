# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_irox]
#     language: python
#     name: conda-env-PROJ_irox-py
# ---

# + [markdown] Collapsed="false" toc-hr-collapsed=true
# # Import Modules

# + [markdown] Collapsed="false"
# ## Notebook Settings

# + Collapsed="false" jupyter={}
# %%capture
# %load_ext autoreload
# %autoreload 2

# + [markdown] Collapsed="false"
# ## Python Packages

# + Collapsed="false"
import os
print(os.getcwd())
import sys

sys.path.insert(
    0,
    os.path.join(
        os.environ["PROJ_irox"],
        "data"))
from proj_data_irox import proj_dir_name, irox_bulk_color_map

# #############################################################################
import copy

import pickle

import chart_studio.plotly as py
import plotly.graph_objs as go

import plotly.offline as py_off
from plotly.offline import (
    init_notebook_mode,
    iplot,
    )

# #############################################################################
from methods import (
    create_pourbaix_plot,
    create_outside_borders,
    create_pourb_entry_outline,
    create_oer_equil_line,
    )

init_notebook_mode(connected=True)


# #############################################################################
from plotting.my_plotly import (
    add_minor_ticks,
    my_plotly_plot,
    add_duplicate_axes)

# + [markdown] Collapsed="false"
# # Script Inputs

# + Collapsed="false"
run_create_entries_script = True

axis_ranges = {
        "x_axis": [0., 14],
        "y_axis": [0., 2.5],

        # "x_axis": [-0.0001, 0.0001],
        # "y_axis": [1.4920084, 1.4923185],
        }

entries_to_remove_list = [
    # "IrO3_TEMP00",
    # "IrO3_TEMP01",
    # "IrO3_TEMP02",
    # "IrO3_TEMP03",
    # "IrO3_TEMP04",
    # "IrO3_TEMP05",
    # "IrO3_TEMP06",
    # "IrO3_TEMP07",
    # "IrO3_TEMP08",
    # "IrO3_TEMP09",
    # "IrO3_TEMP10",
    # "IrO3_TEMP11",
    # "IrO3_TEMP12",
    # "IrO3_TEMP13",
    # "IrO3_TEMP14",
    # "IrO3_TEMP15",
    # "IrO3_TEMP16",
    # "IrO3_TEMP17",
    # "IrO3_TEMP18",
    # "IrO3_TEMP19",
    # "IrO3_TEMP20",
    # "IrO3_TEMP21",


    "IrO3_rutile-like",
    "IrO3_battery",

    # "IrO3_a-AlF3",
    ]

# + [markdown] Collapsed="false"
# # Load Entry Data

# + Collapsed="false"
# # %%capture

if run_create_entries_script:
    tmp = 42
#     !python sc_create_all_entries.py

# +
print("123 SDKISJDIFJSIDPJF")
print(sys.argv)

stoich_i = "AB2"
if sys.argv[-1] == "AB2" or sys.argv[-1] == "AB3":
    stoich_i = sys.argv[-1]
# -

# ## All Pourbaix entries

# + Collapsed="false"
# #############################################################################
path_i = os.path.join(
    "out_data",
    "all_entries_" + stoich_i + ".pickle")
with open(path_i, "rb") as fle:
    all_entries = pickle.load(fle)
# #############################################################################
# -

all_entries

# +
# assert False
# -

# ## Bulk Pourbaix transitions

# +
# path_i = os.path.join(
#     os.environ["PROJ_irox"],
#     "workflow/07_bulk_pourbaix/01_pourbaix_scripts",
#     "out_data/bulk_pourb_transitions.pickle")
# with open(path_i, "rb") as fle:
#     pourb_trans = pickle.load(fle)

# pourb_trans

# + [markdown] Collapsed="false"
# # Process Data for Plotting
# -

all_entries

entries_to_remove_list

# + Collapsed="false"
# %%capture

data_entry_borders = create_pourb_entry_outline(
    entries_to_remove_list=entries_to_remove_list,
    all_entries=all_entries)

data = create_pourbaix_plot(entries=all_entries, axis_ranges=axis_ranges)

oer_equil_line = create_oer_equil_line(axis_ranges=axis_ranges)

data += data_entry_borders
data += [oer_equil_line]

# + [markdown] Collapsed="false"
# # Plotting

# +
from layout import layout

fig = go.Figure(
    data=data,
    layout=layout,
    )


x_axis_range_tmp = [
    axis_ranges["x_axis"][0] + 0.03,
    axis_ranges["x_axis"][1] - 0.03,
    ]


y_axis_range_tmp = [
    axis_ranges["y_axis"][0] + 0.03,
    axis_ranges["y_axis"][1] - 0.03,
    ]

fig.layout.yaxis.range = y_axis_range_tmp
fig.layout.xaxis.range = x_axis_range_tmp

# fig.layout.yaxis.range = axis_ranges["y_axis"]
# fig.layout.xaxis.range = axis_ranges["x_axis"]

# +
# len(layout.annotations)

annotations_new = []
for annot_i in layout.annotations:
    # print(annot_i)
    
    if "<sub>3<" in annot_i.text:
        # print(annot_i)
        
        annotations_new.append(annot_i)

    print("")

# +
# %%capture

fig.layout.update(annotations=annotations_new)
# -

# # Adding 2nd x-axis for minor ticks

# +
shared_dupl_axis = dict(
    tickcolor="black",
    ticklen=3)

ticks_props_new_x = go.layout.XAxis(
    tick0=0,
    **shared_dupl_axis)

ticks_props_new_y = go.layout.YAxis(
    dtick=0.25,
    **shared_dupl_axis)


add_duplicate_axes(
    fig, axis_type='x',
    axis_data=ticks_props_new_x.to_plotly_json())

add_duplicate_axes(
    fig, axis_type='y',
    axis_data=ticks_props_new_y.to_plotly_json())

# +
my_plotly_plot(
    figure=fig,
    # plot_name="bulk_pourbaix_wo_alpha",
    plot_name="bulk_pourbaix_" + stoich_i,
    write_html=True,
    write_png=False,
    png_scale=6.0,
    write_pdf=True,
    write_svg=False,
    try_orca_write=True,
    )

# fig.show()
# -

# # Saving figure to pickle

# + Collapsed="false"
# Pickling data ######################################################
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)
with open(os.path.join(directory, "pourb_fig_" + stoich_i + ".pickle"), "wb") as fle:
    pickle.dump(fig, fle)
# #####################################################################
# -

print(20 * "# # ")
print("All done!")
assert False

fig.show()
