#!/usr/bin/env python

"""Methods for surface energy notebook.

Author(s):Raul A. Flores
"""

# | - Import Modules
import numpy as np

import plotly.graph_objs as go

from oxr_reaction.oxr_rxn import ORR_Free_E_Plot

from surface_energy.surface_energy import surf_e_4
#__|



def make_bulk_stability_shading_dict(
    subplot_num=1,
    x_range=[0., 1.],
    fill_color="black",
    opacity=0.3,
    ):
    """
    """
    # | - make_bulk_stability_shading_dict
    out_dict = {
        'type': 'rect',
        'x0': x_range[0],
        'x1': x_range[1],

        # #####################################################################
        'y0': 0.,
        'y1': 0.05,
        'xref': "x" + str(subplot_num),
        'yref': "paper",

        # 'y0': -0.25,
        # 'y1': 0.5,
        # 'xref': "x" + str(subplot_num),
        # 'yref': "y" + str(subplot_num),
        # #####################################################################

        'line': {
            'color': "black",
            'width': 0.,
            },

        'fillcolor': fill_color,
        # 'opacity': 0.15,  # 0.5
        'opacity': opacity,  # 0.5

        'layer': 'below',
        }

    return(out_dict)
    #__|

def make_color_subplot_list(
    subplot_num=1,
    plot_range=[0, 3],
    bulk_pourb_trans_dict=None,
    irox_bulk_color_map=None,
    opacity=0.3,
    ):
    """
    """
    # | - make_color_subplot_list
    out_list = []

    # | - Bulk Pourbaix Entries
    out_list.append(
        make_bulk_stability_shading_dict(
            subplot_num=subplot_num,
            x_range=[
                plot_range[0],
                bulk_pourb_trans_dict["ir_iro2"]],
            fill_color=irox_bulk_color_map["Ir"],
            opacity=opacity,
            ))

    out_list.append(
        make_bulk_stability_shading_dict(
            subplot_num=subplot_num,
            x_range=[
                bulk_pourb_trans_dict["ir_iro2"],
                bulk_pourb_trans_dict["iro2_iro3"]],
            fill_color=irox_bulk_color_map["IrO2"],
            opacity=opacity,
            ))

    out_list.append(
        make_bulk_stability_shading_dict(
            subplot_num=subplot_num,
            x_range=[
                bulk_pourb_trans_dict["iro2_iro3"],
                bulk_pourb_trans_dict["iro3_iro4-"]],
            fill_color=irox_bulk_color_map["IrO3"],
            opacity=opacity,
            ))

    out_list.append(
        make_bulk_stability_shading_dict(
            subplot_num=subplot_num,
            x_range=[
                bulk_pourb_trans_dict["iro3_iro4-"],
                plot_range[1]],
            fill_color=irox_bulk_color_map["IrO4-"],
            opacity=opacity,
            ))
    #__|

    # | - IrO3 Metastable Polymorph Lines
    height_i = 0.05
    # width_i = 2.
    width_i = 0.02

    # | - Battery IrO3

    key_i = "IrO3_battery"
    out_list.append(
        {
            'type': 'rect',
            'x0': bulk_pourb_trans_dict[key_i + "_low"] - width_i,
            'x1': bulk_pourb_trans_dict[key_i + "_low"] + width_i,
            'y0': 0.,
            'y1': 0.05,
            'xref': "x",
            'yref': "paper",

            'fillcolor': irox_bulk_color_map[key_i],
            'layer': 'below',
            'line': {
                'color': irox_bulk_color_map[key_i],
                # 'width': width_i,
                'width': 0.05,
                },
            }
        )

    out_list.append(
        {
            'type': 'rect',
            'x0': bulk_pourb_trans_dict[key_i + "_high"] - width_i,
            'x1': bulk_pourb_trans_dict[key_i + "_high"] + width_i,
            'y0': 0.,
            'y1': 0.05,
            'xref': "x",
            'yref': "paper",

            'fillcolor': irox_bulk_color_map[key_i],
            'layer': 'below',
            'line': {
                'color': irox_bulk_color_map[key_i],
                # 'width': width_i,
                'width': 0.05,
                },
            }
        )

    # out_list.append(
    #     {
    #         'type': 'line',
    #         'x0': bulk_pourb_trans_dict["iro3_battery_high"],
    #         'x1': bulk_pourb_trans_dict["iro3_battery_high"],
    #         'y0': 0.,
    #         'y1': height_i,
    #         'xref': "x",
    #         'yref': "paper",
    #         'layer': 'below',
    #         'line': {
    #             'color': irox_bulk_color_map["IrO3_battery"],
    #             'width': width_i,
    #             },
    #         }
    #     )
    #__|

    # | - Rutile IrO3
    key_i = "IrO3_rutile-like"
    out_list.append(
        {
            'type': 'rect',
            'x0': bulk_pourb_trans_dict[key_i + "_low"] - width_i,
            'x1': bulk_pourb_trans_dict[key_i + "_low"] + width_i,
            'y0': 0.,
            'y1': 0.05,
            'xref': "x",
            'yref': "paper",

            'fillcolor': irox_bulk_color_map[key_i],
            'layer': 'below',
            'line': {
                'color': irox_bulk_color_map[key_i],
                # 'width': width_i,
                'width': 0.05,
                },
            }
        )

    out_list.append(
        {
            'type': 'rect',
            'x0': bulk_pourb_trans_dict[key_i + "_high"] - width_i,
            'x1': bulk_pourb_trans_dict[key_i + "_high"] + width_i,
            'y0': 0.,
            'y1': 0.05,
            'xref': "x",
            'yref': "paper",

            'fillcolor': irox_bulk_color_map[key_i],
            'layer': 'below',
            'line': {
                'color': irox_bulk_color_map[key_i],
                # 'width': width_i,
                'width': 0.05,
                },
            }
        )
    #__|

    #__|

    return(out_list)
    #__|

def add_convex_hull(
    df_i,
    O_mu_range=None,
    bulk_e_per_atom_dict=None,
    h2_ref=None,
    h2o_ref=None,
    smart_format_dict=None,
    irox_surface_e_color_map=None,
    num_mesh_points=50,
    ):
    """
    """
    # | - add_convex_hull
    list_tmp = []
    for i_cnt, row_i in df_i.iterrows():
        list_tmp.append(
            process_row(
                row_i,
                mesh_eval=True,
                xy_axis=("x", "y"),
                O_mu_range=O_mu_range,
                num_mesh_points=num_mesh_points,
                bulk_e_per_atom_dict=bulk_e_per_atom_dict,
                h2_ref=h2_ref,
                h2o_ref=h2o_ref,
                smart_format_dict=smart_format_dict,
                irox_surface_e_color_map=irox_surface_e_color_map,
                )
            )

    min_e_list = []
    for V_i in zip(*list_tmp):
        min_e = min(V_i)
        min_e_list.append(min_e)
    min_e_list = np.array(min_e_list)
    min_e_list -= 0.003

    trace_ch = go.Scatter(
        x=np.linspace(
            O_mu_range[0], O_mu_range[1],
            num=num_mesh_points, endpoint=True,
            retstep=False, dtype=None,
            ),
        y=min_e_list,
        hoverinfo='skip',
        mode='lines',
        name="Convex Hull",
        line=dict(
            color="black",
            # color="#A92F41",
            width=1.,
            ),
        )

    return(trace_ch)
    #__|

def process_row(
    row_i,
    mesh_eval=False,
    num_mesh_points=50,
    xy_axis=("x", "y"),
    O_mu_range=None,
    bulk_e_per_atom_dict=None,
    h2_ref=None,
    h2o_ref=None,
    smart_format_dict=None,
    irox_surface_e_color_map=None,
    ):
    """

    Args:
        row_i:
        mesh_eval:
            Whether to evaluate the surface energy on a discritized mesh for V
    """
    # | - process_row

    # | - Selecting Bulk Reference
    if row_i["bulk_system"] == "IrO2":
        bulk_e_per_atom = bulk_e_per_atom_dict["IrO2"]
        xy_axis = ("x1", "y1")
    elif row_i["bulk_system"] == "IrO3":
        bulk_e_per_atom = bulk_e_per_atom_dict["IrO3"]
        xy_axis = ("x2", "y2")
    elif row_i["bulk_system"] == "IrO3_rutile-like":
        bulk_e_per_atom = bulk_e_per_atom_dict["IrO3_rutile-like"]
        xy_axis = ("x3", "y3")
    elif row_i["bulk_system"] == "IrO3_battery":
        bulk_e_per_atom = bulk_e_per_atom_dict["IrO3_battery"]
        xy_axis = ("x4", "y4")
    else:
        print("Bulk system not found")
    #__|

    v_list = np.linspace(O_mu_range[0], O_mu_range[1])
    energy_list = []
    for v_i in v_list:
        energy_i = surf_e_4(
            row_i,
            G_H2=h2_ref,
            G_H2O=h2o_ref,
            bias=v_i,
            pH=0.,
            bulk_e_per_atom=bulk_e_per_atom,
            )
        energy_list.append(energy_i)

    if mesh_eval:
        # | - mesh_eval = True
        surf_e_i_vs_V = []
        x_mesh = np.linspace(O_mu_range[0], O_mu_range[1],
            num=num_mesh_points, endpoint=True, retstep=False, dtype=None)
        for V_i in x_mesh:

            surf_e_i = surf_e_4(
                row_i,
                G_H2=h2_ref,
                G_H2O=h2o_ref,
                bias=V_i,
                pH=0.,
                bulk_e_per_atom=bulk_e_per_atom,
                )

            surf_e_i_vs_V.append(surf_e_i)

        return(surf_e_i_vs_V)
        #__|

    # | - Constructing Series Properties
    series_format_dict_i = ORR_Free_E_Plot().__create_smart_format_dict__(
        row_i.to_dict(),
        smart_format_dict,
        )

    dash_i = series_format_dict_i.get("dash", "solid")
#     color_i = irox_bulk_color_map[row_i["bulk_system"]]

    system_i = row_i["bulk_system"] + "_" + row_i["coverage_type"]
    color_i = irox_surface_e_color_map[system_i]

    name_i = row_i["name_i_3"]
    name_i = name_i.replace("IrO2_", "")
    name_i = name_i.replace("IrO3_rutile-like_", "")
    name_i = name_i.replace("IrO3_battery_", "")
    name_i = name_i.replace("IrO3_", "")
    #__|

    # | - Plotly Scatter Trace
    trace = go.Scatter(
        x=v_list,
        y=energy_list,

        # x=[O_mu_range[0], O_mu_range[1]],
        # y=[surf_e_LHS, surf_e_RHS],

        xaxis=xy_axis[0],
        yaxis=xy_axis[1],
        mode='lines',
        hoverinfo='text',
        name=name_i,
        text=name_i,
        line=dict(
            color=color_i,
            dash=dash_i,
            width=1.,
            ),
        )
    #__|

    return(trace)




    # | - __old__
    # if row_i["bulk_system"] == "IrO2":
    #     traces_IrO2.append(trace)
    # elif row_i["bulk_system"] == "IrO3":
    #     traces_IrO3.append(trace)
    # elif row_i["bulk_system"] == "IrO3_rutile-like":
    #     traces_IrO3_rutile_like.append(trace)
    # elif row_i["bulk_system"] == "IrO3_battery":
    #     traces_IrO3_battery.append(trace)
    #__|

    #__|






# | - __old__

# def surf_e_4(
#     row_i,
#     G_H2=0.,
#     G_H2O=0.,
#     bias=0.,
#     pH=0.,
#     bulk_e_per_atom=None,
#     get_e0_from_row=False,
#     norm_mode="area",  # 'area' or 'atoms'
#     num_atoms=None,
#     ):
#     """
#     Calculate surface energy assuming a water reference state
#     and using the computational hydrogen electrode.
#     """
#     # | - surf_e_4
#
#     # | - Read info from row_i
#     # COMBAK This shouldn't be hard coded in
#     metal = "Ir"
#
#     atoms_i = row_i.get("init_atoms", default=None)
#     elec_energy = row_i.get("elec_energy", default=0.)
#     nonstoich_Os = row_i.get("nonstoich_Os", default=0)
#     elems_dict = row_i.get("elem_num_dict", default={})
#
#     if bulk_e_per_atom is None:
#         bulk_e_per_atom = row_i.get("bulk_e_per_atom_DFT", default=0.)
#
#     # print("bulk_e_per_atom: ", bulk_e_per_atom)
#
#     N_stoich_in_slab = elems_dict[metal] + elems_dict.get("O", 0) - nonstoich_Os
#     nonstoich_Hs = elems_dict.get("H", 0)
#     #__|
#
#     # | - Calculate Standard State Surface Energy
#     if "surf_e_0" in row_i:
#         surf_e_0 = row_i.get("surf_e_0", default=0.)
#     else:
#         surf_e_0 = 0. + \
#             +elec_energy + \
#             -N_stoich_in_slab * bulk_e_per_atom + \
#             -nonstoich_Os * (G_H2O - G_H2) + \
#             -nonstoich_Hs * (G_H2) + \
#             +0.
#
#         if norm_mode == "area":
#             surf_e_0 = surf_e_0 / (2 * row_i["slab_area"])
#         elif norm_mode == "atoms":
#             if num_atoms is not None:
#                 surf_e_0 = surf_e_0 / (2 * num_atoms)
#     #__|
#
#     # | - Calculate V, pH Dependant Surface Energy
#     slope = 2 * nonstoich_Os - nonstoich_Hs
#
#     surf_e = 0. + \
#         +surf_e_0 + \
#         +(slope * 0.0591 * pH) / (2 * row_i["slab_area"]) + \
#         -(slope * bias) / (2 * row_i["slab_area"]) + \
#         +0.
#
# #     surf_e = surf_e / (2 * row_i["slab_area"])
#     #__|
#
#     return(surf_e)
#     #__|


#__|
