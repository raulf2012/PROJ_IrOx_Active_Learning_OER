
#| - Import Modules
import os
import pickle
import pandas as pd

import plotly.plotly as py
import plotly.graph_objs as go
#__|

#| - Script Inputs
# duration_long = 500
# duration_short = 300

# duration_long = 800
# duration_short = 600

duration_long = 1000 * 2
duration_short = 800 * 2
#__|

# #############################################################################
# #############################################################################


def get_data(stoich_i=None):
    """
    """

    #| - Read Data
    # stoich_i = "iro2"
    # stoich_i = "iro2"

    path_i = os.path.join(
        "..",
        "data",
        "models_list_" + stoich_i + ".pickle")
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

    if stoich_i == "iro2":
        global_y_min = -1.5
        global_y_max = 2.4

        # Shortening iterations to stop on "nicer" iter
        print("TEMP TEMP TEMP | dsjfisdjf77sfsd")
        print(len(models_list))
        # models_list = models_list[:-3]
        models_list = models_list[:-8]
        print(len(models_list))

    elif stoich_i == "iro3":
        global_y_min = -1.0
        global_y_max = 1.3

        # Shortening iterations to stop on "nicer" iter
        models_list = models_list[:-5]


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

    # len(models_list[0]) + 5

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

            'range': [0 - 5, len(models_list[0]) + 5],
            'autorange': False,
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

            'range': [global_y_min, global_y_max],
            # 'range': [-1.5, 2.4],
            'autorange': False,
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
        'transition': {
            'duration': duration_short,  # TEMP_DURATION
            'easing': 'cubic-in-out'},
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

    #__|

    data = {
        "models_list": models_list,
        "layout": layout,
        # "get_trace_i": get_trace_i,
        "sliders_dict": sliders_dict,
        "get_slider_step_i": get_slider_step_i,
        }

    return(data)

# #############################################################################
# #############################################################################

def get_trace_i(
    grid,
    column_names,
    index=None,
    num_structures=None,
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
            data_i["prediction_not_ref"] = grid.get_column(col_j)
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

    #| - Creating Horizontal Line Trace
    # x_array = [0 - 50, len(models_list[0]) + 50]
    x_array = [0 - 50, num_structures + 50]

    # #########################################################################
    df_tmp = pd.merge(
        data_i["computed"],
        data_i["prediction_not_ref"].data,
        left_index=True, right_index=True)
    min_computed_y = df_tmp[df_tmp["computed"] == True]["prediction"].min()

    y_array = 2 * [min_computed_y]
    trace_tmp_0 = go.Scatter(
        x=x_array,
        y=y_array,
        mode="lines",
        line=dict(
            color=("black"),
            width=2,
            dash="solid",
            )
        )

    # #########################################################################
    min_y = data_i["prediction_not_ref"].data.min()

    y_array = 2 * [min_y]
    trace_tmp_1 = go.Scatter(
        x=x_array,
        y=y_array,
        mode="lines",
        line=dict(
            color=("grey"),
            width=2,
            dash="dash",
            )
        )

    #__|

    #| - Setting Individual Marker Settings
    marker_size_list = []
    marker_line_size_list = []
    marker_line_color_list = []
    for i in data_i["computed"].tolist():
        if i is True:
            marker_size_list.append(10)
            marker_line_color_list.append("red")
            marker_line_size_list.append(1.5)
        elif i is False:
            marker_size_list.append(5)
            marker_line_color_list.append("rgb(0, 0, 0)")
            marker_line_size_list.append(0.05)
    #__|

    #| - go.Scatter instance
    trace_i = go.Scatter(
        xsrc=data_i["x_axis_index"],
        ysrc=data_i["prediction"],

        error_y=dict(
            type='data',
            array=data_i["uncertainty"],
            visible=True,
            # thickness=0.5,
            thickness=0.2,
            width=1.5,
            # color="rgba(230,230,230,1.0)",
            color="rgba(120,120,120,1.0)",
            ),
        name=index,
        mode="markers",

        hoverinfo="text",
        text=data_i["uncertainty"].index,
        marker={
            "opacity": 0.9,
            "size": marker_size_list,
            "color": data_i["color_order_rank"],
            "colorscale": "Viridis",
            "line": {
                # "width": 0.05,
                # "width": 0.5,
                "width": marker_line_size_list,
                # "color": "rgb(0, 0, 0)",
                "color": marker_line_color_list,
                },
            },
        )
    #__|

    return([
        trace_i,
        trace_tmp_0,
        trace_tmp_1,
        ])
    #__|
