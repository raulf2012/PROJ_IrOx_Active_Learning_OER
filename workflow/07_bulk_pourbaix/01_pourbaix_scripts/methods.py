#!/usr/bin/env python

"""Methods to create plotly plot of bulk Pourbaix diagram.

Author: Raul A. Flores
"""

#| - Import Modules
import numpy as np

import plotly.graph_objs as go
    # import plotly.graph_objs as go

from pymatgen.analysis.pourbaix_diagram import (
    PourbaixDiagram,
    PourbaixPlotter,
    generate_entry_label,
    )

from proj_data_irox import (
    pymatgen_to_my_naming_convention,
    irox_bulk_color_map,
    )
#__|


plot_title = None
tick_lab_size = 8 * (4. / 3.)
axes_lab_size = 9 * (4. / 3.)
# legend_size = 18

# font_family="Computer Modern"  # "Courier New, monospace"
font_family = "Arial"  # "Courier New, monospace"


def create_pourbaix_plot(
    entries=None,
    axis_ranges=None,
    ):
    """
    """
    #| - create_pourbaix_plot
    all_entries = entries

    pourbaix = PourbaixDiagram(all_entries)  # comp_dict={'Ru':0.5,'Y':0.5})
    # plotter = PourbaixPlotter(pourbaix)

    stable_domains, stable_domain_vertices = pourbaix.get_pourbaix_domains(
        all_entries,
        limits=[
            axis_ranges["x_axis"],
            axis_ranges["y_axis"],
            ],
        )

    # print("9s9s9s9s9s09dfsd)(*())")
    # print(type(stable_domain_vertices))
    # print("9s9s9s9s9s09dfsd)(*())")
    # [[-2, 16], [-4, 4]]

    # for entry, vertices in plotter._pd._stable_domain_vertices.items():
    data = []
    for entry, vertices in stable_domain_vertices.items():
        trace_i, trace_center = process_sys(entry, vertices)
        data.append(trace_i)
        data.append(trace_center)

    layout = pourbaix_plot_layout(axis_ranges)

    return(data, layout)
    #__|

def process_sys(entry, vertices):
    """
    """
    #| - process_sys
    x, y = np.transpose(np.vstack([vertices, vertices[0]]))
    center = np.average(vertices, axis=0)


    #| - TMP

    if hasattr(entry.entry, 'attribute'):
        att_dict = entry.entry.attribute
        name_i = att_dict.get("full_name", "N/A")
    else:
        # name_i = "TMPTMP"

        name_i = generate_entry_label(entry)
        b = "_${}^"
        for char in b:
            name_i = name_i.replace(char, "")

    # print(name_i)

    color_i = irox_bulk_color_map.get(
        pymatgen_to_my_naming_convention.get(name_i, name_i),
        "black",
        )

    # print(name_i)
    # print(color_i)
    #__|

    #| - Scatter Instances
    trace_i = go.Scatter(
        x=x,
        y=y,
        mode='lines',
        name=name_i,
        showlegend=True,
        fill='toself',
        fillcolor=color_i,
        hoverinfo="none",
        line={
            "color": "black",
            }
        )

    trace_center = go.Scatter(
        x=[center[0]],
        y=[center[1]],
        hoverinfo="name",
        name=name_i,
        mode='lines+text',
        text=[name_i],
        textposition='bottom center',
        textfont={
            "family": font_family,
            "color": "black",
            "size": axes_lab_size,
            },
        )
    #__|

    return(trace_i, trace_center)
    #__|

def pourbaix_plot_layout(
    axis_ranges,
    ):
    """
    """
    #| - pourbaix_plot_layout

    plot_title = None
    tick_lab_size = 8 * (4. / 3.)
    axes_lab_size = 9 * (4. / 3.)
    # legend_size = 18

    # font_family="Computer Modern"  # "Courier New, monospace"
    font_family = "Arial"  # "Courier New, monospace"


    layout = go.Layout(

        font={
            "family": font_family,
            "color": "black",
            },

        title=plot_title,
        titlefont=None,

        xaxis={
            "title": "pH",
            "range": axis_ranges["x_axis"],
            "zeroline": False,
            "showline": True,
            "mirror": 'ticks',
            "linecolor": 'black',
            "showgrid": False,

            "ticks": 'inside',
            "tick0": 0,
            "tickcolor": 'black',

            "dtick": 1,

            "ticklen": 2,
            "tickwidth": 1,

            "titlefont": dict(size=axes_lab_size),
            "tickfont": dict(
                size=tick_lab_size,
                ),
            # "titlefont": dict(
            #     # family='Courier New, monospace',
            #     size=18,
            #     color="black",
            #     )
            },

        yaxis={
            "title": "E(V)",

            "range": axis_ranges["y_axis"],
            "zeroline": False,
            "showline": True,
            "mirror": 'ticks',
            "linecolor": 'black',
            "showgrid": False,
            "ticks": 'inside',
            "tick0": 0,
            "tickcolor": 'black',
            # "dtick": 0.25,
            "ticklen": 2,
            "tickwidth": 1,

            "titlefont": dict(size=axes_lab_size),
            "tickfont": dict(
                size=tick_lab_size,
                ),
            },

        margin={
            "b": 40.,
            "l": 40.,
            "r": 40.,
            "t": 40.,
            },

        # legend=None,
        showlegend=False,
        width=1 * 9. * 37.795275591,
        height=1 * 6. * 37.795275591,
        )

    return(layout)
    #__|

# #############################################################################

def create_outside_borders(axis_ranges=None):
    """
    """
    #| - create_outside_borders
    x_border = [
        axis_ranges["x_axis"][0],
        axis_ranges["x_axis"][1],
        axis_ranges["x_axis"][1],
        axis_ranges["x_axis"][0],

        axis_ranges["x_axis"][0],
        ]

    y_border = [
        axis_ranges["y_axis"][0],
        axis_ranges["y_axis"][0],
        axis_ranges["y_axis"][1],
        axis_ranges["y_axis"][1],

        axis_ranges["y_axis"][0],
        ]

    outside_border = go.Scatter(
        x=x_border,
        y=y_border,
        mode='lines',
        marker=dict(
            color="black",
            line=dict(
                width=4,
                color="black",
                )
            )
        )

    return(outside_border)
    #__|

def create_pourb_entry_outline(
    entries_to_remove_list=None,
    all_entries=None,
    ):
    """
    """
    #| - create_pourb_entry_outline
    # Common species that all IrO3 species will be plotted with
    ir_entry = get_base_spec("Ir", all_entries)
    iro2_entry = get_base_spec("IrO2", all_entries)
    ir_ion_entry = get_base_spec("IrO4-", all_entries)

    out_dict = get_spec_entries(entries_to_remove_list, all_entries)

    # Create entries list, ready for Pourbaix plotting
    pourb_dict = {}
    for key, value in out_dict.items():
        pourb_dict[key] = [value, ir_entry, iro2_entry, ir_ion_entry]

    # #############################################################################
    # #############################################################################

    data_tmp = []
    for sys_i in entries_to_remove_list:
        color_i = irox_bulk_color_map[sys_i]

        tmp = pourb_dict.get(sys_i, None)
        if tmp is None:
            continue

        pourbaix_i = PourbaixDiagram(tmp)

        # pourbaix_i = PourbaixDiagram(pourb_dict[sys_i])
        pourb_vert_i, pourb_domains_i = pourbaix_i.get_pourbaix_domains(
            pourb_dict[sys_i],
            limits=None,
            )

        x_list = []
        y_list = []
        for side_i in pourb_domains_i[pourb_dict[sys_i][0]]:
            x_list.append(side_i[0])
            y_list.append(side_i[1])

        x_list.append(x_list[0])
        y_list.append(y_list[0])

        data_i = go.Scatter(
            x=x_list,
            y=y_list,
            mode='lines',

            line=dict(
                width=3,
                color=color_i,
                ),


            # marker=dict(
            #     # size = 10,
            #     color=color_i,
            #     line=dict(
            #         # width=4,
            #         width=10,
            #         # color=color_i,
            #         color="black",
            #         )
            #     )

            )

        data_tmp.append(data_i)

    return(data_tmp)
    #__|

def get_base_spec(base_species, all_entries):
    """
    """
    #| - get_base_spec
    entry_out = None
    for entry_i in all_entries:

        if hasattr(entry_i.entry, 'attribute'):
            att_dict = entry_i.entry.attribute
            name_i = att_dict.get("full_name", "N/A")
        else:
            # name_i = "TMPTMP"

            name_i = generate_entry_label(entry_i)
            b = "_${}^"
            for char in b:
                name_i = name_i.replace(char, "")

        if name_i == base_species:
            entry_out = entry_i

#         if hasattr(i.entry, 'attribute'):
#             name_i = i.entry.attribute["full_name"]

#             if name_i == base_species:
#                 entry_i = i

    return(entry_out)
    #__|

def get_spec_entries(entries_list, all_entries):
    """
    """
    #| - get_spec_entries
    out_dict = {}
    for j in entries_list:
        for i in all_entries:
            if hasattr(i.entry, 'attribute'):
                name_i = i.entry.attribute["full_name"]
                if name_i == j:
                    out_dict[name_i] = i
    return(out_dict)
    #__|


# #############################################################################

def create_oer_equil_line(axis_ranges=None):
    """Create plotly data trace for OER equilibrium line on Pourbaix plot.


    """
    #| - create_oer_equil_line
    # import plotly.graph_objs as go
    PREFAC = 0.0591

    oer_equil_line = go.Scatter(
        x=[axis_ranges["x_axis"][0], axis_ranges["x_axis"][1]],
        y=[
            -PREFAC * axis_ranges["x_axis"][0] + 1.23,
            -PREFAC * axis_ranges["x_axis"][1] + 1.23,
            ],
        mode='lines',
        name='lines',

        # color="black",
        # width=1,
        # dash="dot",

        line=dict(
            color=('black'),
            width=0.75,
            dash='dot'),
        )

    return(oer_equil_line)
    #__|
