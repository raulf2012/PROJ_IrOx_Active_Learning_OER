#!/usr/bin/env python

"""My plot settings and methods.

Author: Raul A. Flores
"""

# | - Import Modules
import os

import copy

# Plotly imports
import plotly

import chart_studio.plotly as py
import plotly.graph_objs as go

from plotly import io as pyio
#__|

from plotting.my_plotly import get_xy_axis_info

def add_duplicate_axes(
    fig,
    axis_type="x",  # 'x' or 'y'
    axis_data=dict(),
    axis_num_list=None,
    ):
    """
    """
    # | - add_duplicate_axes

    # This is necessary to make sure that the original traces are still visible after adding the new traces
    fig.update_layout(
        # paper_bgcolor="white",
        plot_bgcolor="rgba(255,255,255,0.)",
        )

    axis_info_dict = get_xy_axis_info(fig)[axis_type]
    num_of_axis = axis_info_dict["num_of_axis"]

    # #########################################################################
    if axis_num_list is None:
        # axis_list = axis_info_dict["axis_list"]
        axis_num_list = axis_info_dict["axis_num_list"]


    # axis_num_list_new = [i + len(axis_num_list) for i in axis_num_list]
    axis_num_list_new = [i + num_of_axis for i, j in enumerate(axis_num_list)]

    # [(i, j) for i, j in enumerate(mylist)]

    iterator = enumerate(zip(axis_num_list, axis_num_list_new))
    for i_cnt, (old_index, new_index) in iterator:
        old_Axis = fig.layout[axis_type + "axis" + str(old_index)]

        new_axis = copy.deepcopy(old_Axis)
        new_axis = new_axis.update(
            # dtick=0.1,
            showticklabels=False,
            title=dict(
                font=None,
                standoff=None,
                text="",
                ),
            )

        new_axis = new_axis.update(**axis_data)

        axis_key = axis_type + "axis" + str(new_index)
        new_layout = go.Layout({
            axis_key: new_axis,
            })

        fig.update_layout(new_layout)

        fig.add_scatter(
            **go.Scatter({
                axis_type + "axis": axis_type + str(new_index)
                }).to_plotly_json())

    #__|
