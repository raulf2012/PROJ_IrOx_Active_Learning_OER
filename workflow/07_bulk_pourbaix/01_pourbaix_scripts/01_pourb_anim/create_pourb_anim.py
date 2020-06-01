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

# # Import Modules

# +
import os
import pickle

import numpy as np
# -

duration_long = 1000 * 3
duration_short = 800 * 3

# +
# #############################################################################
path_i = os.path.join(
    "out_data",
    "pourb_fig_0.pickle")
with open(path_i, "rb") as fle:
    fig_0 = pickle.load(fle)
# #############################################################################


# #############################################################################
path_i = os.path.join(
    "out_data",
    "pourb_fig_1.pickle")
with open(path_i, "rb") as fle:
    fig_1 = pickle.load(fle)
# #############################################################################

# +
tmp_x = fig_1.data[2]["x"]

first_x = tmp_x[0]
tmp_x = tmp_x[1:]
x_new = np.append(tmp_x, first_x)

fig_1.data[2]["x"] = x_new



tmp_y = fig_1.data[2]["y"]

first_y = tmp_y[0]
tmp_y = tmp_y[1:]
y_new = np.append(tmp_y, first_y)

fig_1.data[2]["y"] = y_new

# +
x_new = [1.40000000e+01,  1.22371036e-02, -8.88178420e-16,  8.88178420e-16, 1.40000000e+01,  1.40000000e+01]
y_new = [0.38893608,      1.4911718 ,      1.49261823,      2.5       ,     2.5       ,      0.38893608]

# x_new = [0, 0, 10, 10, 0]
# y_new = [4, 6, 6,  4,  4]

fig_0.data[6]["x"] = x_new
fig_0.data[6]["y"] = y_new

# +
# x_new = [1.4000e+01,      1.26107660e+01,  1.77635684e-15,  0.0000e+00,     1.40000000e+01,  1.40000000e+01]
# y_new = [0.38893608,      0.498407720000,  1.989000260000,  2.5       ,     2.5           ,  0.38893608]

x_new = [1.40000000e+01,  1.22371036e-02, -8.88178420e-16,  8.88178420e-16, 1.40000000e+01,  1.40000000e+01]
y_new = [0.38893608,      1.4911718 ,      1.49261823,      2.5       ,     2.5       ,      0.38893608]

# x_new = [0, 0, 10, 10, 0]
# y_new = [4, 6, 6,  4,  4]

fig_1.data[6]["x"] = x_new
fig_1.data[6]["y"] = y_new

# +
# on_top_traces = []
# for trace_i in fig_1.data:
#     name = trace_i.name

#     if name is not None:
#         if "AlF3" in name:
#             tmp = 42
#             on_top_traces.append(trace_i)
            

# new_data =fig_1.data + tuple(on_top_traces)


# fig_1.data = new_data

# +
fig_data_list = list(fig_1.data)

fig_data_list.append(fig_data_list.pop(4)) 

fig_1.data = tuple(fig_data_list)


fig_data_list = list(fig_0.data)

fig_data_list.append(fig_data_list.pop(4)) 

fig_0.data = tuple(fig_data_list)

# + jupyter={}
import plotly.graph_objs as go


def get_layout(
    duration_long=1000 * 2,
    duration_short=800 * 2,
    ):
    """
    """
    # | - get_layout

    # | - updatemenus
    updatemenus = [
        {
            'buttons': [
                {
                    'args': [
                        None,
                        {
                            'frame': {
                                'duration': duration_long,  # TEMP_DURATION
                                'redraw': False,
                                },
                            'fromcurrent': True,
                            'transition': {
                                'duration': duration_short,  # TEMP_DURATION
                                'easing': 'quadratic-in-out',
                                }
                            }
                        ],
                    'label': 'Play',
                    'method': 'animate'
                    },
                {
                    'args': [
                        [None],
                        {
                            'frame': {
                                'duration': 0,
                                'redraw': False,
                                },
                            'mode': 'immediate',
                            'transition': {'duration': 0},
                            }
                        ],
                    'label': 'Pause',
                    'method': 'animate'
                    }
                ],
            'direction': 'left',
            'pad': {'r': 10, 't': 87},
            'showactive': False,
            'type': 'buttons',
            'x': 0.1,
            'xanchor': 'right',
            'y': 0,
            'yanchor': 'top'
            }
        ]

    #__|

    layout = go.Layout(
        # title='Material Discovery Training',
        showlegend=False,
        font=dict(
            # family='Courier New, monospace',

            family='Arial',
            size=20,
            color='black'
            ),

        xaxis={
            'title': 'Candidate Structures',

            # 'range': [0 - 5, len(models_list[0]) + 5],
            # 'autorange': False,
            'autorange': True,

            'showgrid': False,
            'zeroline': False,
            'showline': True,
            'ticks': '',
            'showticklabels': True,
            'mirror': True,
            'linecolor': 'black',

            },

        yaxis={
            'title': 'Formation Energy (eV)',

            # 'range': [global_y_min, global_y_max],
            # 'range': [-1.5, 2.4],
            # 'autorange': False,
            'autorange': True,
            'fixedrange': False,

            'showgrid': False,
            'zeroline': True,
            'showline': True,
            'ticks': '',
            'showticklabels': True,
            'mirror': True,
            'linecolor': 'black',

            },

        updatemenus=updatemenus,
        )

    return(layout)
    #__|

    
def get_sliders_init_dict(duration_short):
    """
    """
    # | - get_sliders_init_dict
    sliders_dict = {
        # 'active': 0,
        'active': 0,
        'yanchor': 'top',
        'xanchor': 'left',
        'currentvalue': {
            'font': {'size': 20},
            'prefix': 'Loop #:',
            'visible': True,
            'xanchor': 'right'
            },
        'transition': {
            'duration': duration_short,  # TEMP_DURATION
            'easing': 'cubic-in-out'},
        'pad': {'b': 10, 't': 50},
        'len': 0.8,
        'x': 0.1,
        'y': 0,
        'steps': [],
        }

    return(sliders_dict)
    #__|

def get_slider_step_i(i_cnt, duration_short):
    """
    """
    # | - get_slider_step_i
    slider_step_i = {
        'args': [
            [str(i_cnt)],
            {
                'frame': {
                    'duration': duration_short,  # TEMP_DURATION
                    'redraw': False},
                'mode': 'immediate',
                'transition': {
                    'duration': duration_short,  # TEMP_DURATION
                    },
                },
            ],
        'label': str(i_cnt),
        'method': 'animate',
        }

    return(slider_step_i)
    #__|

# TEMP


# + jupyter={}
layout_anim = get_layout(
    duration_long=duration_long,
    duration_short=duration_short)
sliders_dict = get_sliders_init_dict(duration_short)

frames = []; data = []
# for i_cnt, df_i in enumerate(df_list):

# + jupyter={}
fig_0.layout.template = None
layout_anim = layout_anim.update(fig_0.layout)

# + jupyter={}
i_cnt = 0
traces_i = fig_0.data
# traces_i = get_traces(df_i)

# #################################################################
if i_cnt == 0: data.extend(traces_i)
data_i = []; data_i.extend(traces_i)
frame_i = go.Frame(data=data_i, name=str(i_cnt))
frames.append(frame_i)
slider_step_i = get_slider_step_i(i_cnt, duration_short)
sliders_dict['steps'].append(slider_step_i)

# + jupyter={}
i_cnt = 1
traces_i = fig_1.data
# traces_i = get_traces(df_i)

# #################################################################
if i_cnt == 0: data.extend(traces_i)
data_i = []; data_i.extend(traces_i)
frame_i = go.Frame(data=data_i, name=str(i_cnt))
frames.append(frame_i)
slider_step_i = get_slider_step_i(i_cnt, duration_short)
sliders_dict['steps'].append(slider_step_i)

# + jupyter={}
layout_anim["showlegend"] = False

fig = go.Figure(
    data=data,
    layout=layout_anim,
    frames=frames)
fig['layout']['sliders'] = [sliders_dict]

fig.show()
# -

assert False

fig_0.data[6]

fig_1.data[6]

# +
# tmp = fig_1.data[2]["x"]

# # tmp = 

# # tmp = tmp[1:].append(tmp[0])
# # # + [tmp[0]]
