
#| - Import Modules
import os
import pickle

import plotly.graph_objs as go
#__|

#| - Script Inputs
# duration_long = 500
# duration_short = 300
#__|

# #############################################################################
# #############################################################################

#| - Read Data
# path_i = os.path.join("..", "models_list.pickle")
path_i = os.path.join(
    "..",
    "data",
    "models_list_iro2.pickle")
with open(path_i, "rb") as fle:
    models_list = pickle.load(fle)

# models_list = models_list[-5:]
#__|

#| - Process Data

#| - Calculate y-axis range
max_ys = []
min_ys = []
for model_i in models_list:
    max_y_i = (model_i["prediction"] + model_i["uncertainty"]).max()
    min_y_i = (model_i["prediction"] - model_i["uncertainty"]).min()

    max_ys.append(max_y_i)
    min_ys.append(min_y_i)

global_y_max = max(max_ys)
global_y_min = min(min_ys)
#__|

reference_model = models_list[-1]
reference_model = reference_model.sort_values("prediction")

main_index_order = reference_model.index.tolist()

reference_model["color_order_rank"] = [i / len(reference_model) for i in range(len(reference_model))]

color_mapping_dict = dict(zip(
    reference_model.index,
    reference_model["color_order_rank"]))

def method(row_i):
    color_rank_i = color_mapping_dict[row_i.name]
    return(color_rank_i)


models_list_new = []
for model_i in models_list:
    # Sort model by energy temporarily
    model_i = model_i.sort_values("prediction", inplace=False)

    model_i["x_axis_index"] = [i for i in range(len(model_i))]

    model_i["color_order_rank"] = model_i.apply(method, axis=1)

    model_i = model_i.reindex(main_index_order)
    models_list_new.append(model_i)

models_list = models_list_new
#__|

# #############################################################################
# #############################################################################

#| - Layout

#| - updatemenus
updatemenus = [
    {
        'buttons': [
            {
                'args': [
                    None,
                    {
                        'frame': {
                            'duration': 500, 'redraw': False,
                            },
                        'fromcurrent': True,
                        'transition': {
                            'duration': 300, 'easing': 'quadratic-in-out'
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

    xaxis={
        'title': 'Candidate Structures',

#         'range': [0, 2],
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
#__|

# #############################################################################
# #############################################################################

#| - sliders
sliders_dict = {
#     'active': 0,
    'active': 0,
    'yanchor': 'top',
    'xanchor': 'left',
    'currentvalue': {
        'font': {'size': 20},
        'prefix': 'Loop #:',
        'visible': True,
        'xanchor': 'right'
    },
    'transition': {'duration': 300, 'easing': 'cubic-in-out'},
    'pad': {'b': 10, 't': 50},
    'len': 0.8,
    'x': 0.1,
    'y': 0,
    'steps': [],
    }

def get_slider_step_i(i_cnt):
    """
    """
    #| - get_slider_step_i
    slider_step_i = {
        'args': [
            [str(i_cnt)],
            {
                'frame': {'duration': 300, 'redraw': False},
                'mode': 'immediate',
                'transition': {'duration': 300},
                },
            ],
        'label': str(i_cnt),
        'method': 'animate',
        }

    return(slider_step_i)
    #__|
#__|

# #############################################################################
# #############################################################################

def get_trace_i(
    grid,
    column_names,
    index=None,
    ):
    """
    # 'prediction_4',
    # 'uncertainty_4',
    # 'computed_4',
    # 'x_axis_index_4',
    # 'color_order_rank_4',
    """
    #| - get_trace_i
    ref_cols = ["prediction", "x_axis_index"]

    data_i = {}
    for col_j in column_names:
        if "prediction_" in col_j:
            data_i["prediction"] = grid.get_column_reference(col_j)
        elif "uncertainty_" in col_j:
            data_i["uncertainty"] = grid.get_column(col_j).data
        elif "computed_" in col_j:
            data_i["computed"] = grid.get_column(col_j).data
        elif "x_axis_index_" in col_j:
            data_i["x_axis_index"] = grid.get_column_reference(col_j)
        elif "color_order_rank_" in col_j:
            data_i["color_order_rank"] = grid.get_column(col_j).data
        else:
            tmp = 42


    marker_size_list = []
    for i in data_i["computed"].tolist():
        if i is True: marker_size_list.append(12)
        elif i is False: marker_size_list.append(5)


    trace_i = go.Scatter(
        xsrc=data_i["x_axis_index"],
        ysrc=data_i["prediction"],

        error_y=dict(
            type='data',
            array=data_i["uncertainty"],
            visible=True,
            thickness=1.,
            width=3.,
            # color="rgba(230,230,230,1.0)",
            color="rgba(120,120,120,1.0)",
            ),
        name=index,
        mode="markers",

        hoverinfo="text",
        text=data_i["uncertainty"].index,
        marker={
            "size": marker_size_list,
            "color": data_i["color_order_rank"],
            "colorscale": "Viridis",
            "line": {
                "width": 0.05,
                "color": "rgb(0, 0, 0)",
                },
            },
        )

    return(trace_i)
    #__|










#| - __old__
#
# updatemenus = [
#     {
#         'buttons': [
#             {
#                 'args': [
#                     None,
#                     {
#                         'frame': {'duration': 500, 'redraw': False},
#                         'fromcurrent': True,
#                         'transition': {
#                             'duration': 300,
#                             'easing': 'quadratic-in-out',
#                             }
#                         },
#                     ],
#                 'label': 'Play',
#                 'method': 'animate'
#                 },
#
#             {
#                 'args': [
#                     [None],
#                     {
#                         'frame': {'duration': 0, 'redraw': False},
#                         'mode': 'immediate',
#                         'transition': {
#                             'duration': 0,
#                             },
#                         }
#                     ],
#                 'label': 'Pause',
#                 'method': 'animate'
#                 }
#             ],
#
#
#
#         'direction': 'left',
#         'pad': {'r': 10, 't': 87},
#         'showactive': False,
#         'type': 'buttons',
#         'x': 0.1,
#         'xanchor': 'right',
#         'y': 0,
#         'yanchor': 'top'
#         }
#     ]

# [
#         {
#             'buttons': [
#                 {
#                     'args': [None],
#                     'label': 'Play',
#                     'method': 'animate',
#                     }
#                 ],
#             'pad': {'r': 10, 't': 87},
#             'showactive': False,
#             'type': 'buttons',
#             }
#         ]


# import pandas as pd
# from plotly.offline import init_notebook_mode, iplot
# from IPython.display import display, HTML
# init_notebook_mode(connected=True)


# figure = go.Figure()

#
# #| - Layout
#
#
# #| - Buttons
# button_play = go.layout.updatemenu.Button(
#     args=[
#         None,
#         {
#             'frame': {'duration': duration_long, 'redraw': False},
#             'fromcurrent': True,
#             'transition': {'duration': duration_short, 'easing': 'quadratic-in-out'}
#             }
#         ],
#     label="Play",
#     method="animate",
#     )
# # #############################################################################
# # #############################################################################
# button_pause = go.layout.updatemenu.Button(
#     args=[
#         [None],
#         {
#             'frame': {'duration': 0, 'redraw': False},
#             'mode': 'immediate',
#             'transition': {'duration': 0}
#             }
#         ],
#     label="Pause",
#     method="animate",
#     )
# #__|
#
# # #############################################################################
# # #############################################################################
#
# layout = go.Layout(
#     xaxis={'range': [0, 800], 'title': 'Life Expectancy'},
#     yaxis={'title': 'GDP per Capita'},
#     hovermode="closest",
#     updatemenus= [go.layout.Updatemenu(
#         buttons=[button_play, button_pause],
#         direction="left",
#         pad={'r': 10, 't': 87},
#         showactive=False,
#         type="buttons",
#         x=0.1,
#         xanchor="right",
#         y=0,
#         yanchor="top",
#         )],
#     )
#
# figure.layout = layout
#
# sliders_dict = {
#     'active': 0,
#     'yanchor': 'top',
#     'xanchor': 'left',
#     'currentvalue': {
#         'font': {'size': 20},
#         'prefix': 'Year:',
#         'visible': True,
#         'xanchor': 'right'
#     },
#     'transition': {'duration': duration_short, 'easing': 'cubic-in-out'},
#     'pad': {'b': 10, 't': 50},
#     'len': 0.9,
#     'x': 0.1,
#     'y': 0,
#     'steps': [],
#     }
# #__|

#__|
