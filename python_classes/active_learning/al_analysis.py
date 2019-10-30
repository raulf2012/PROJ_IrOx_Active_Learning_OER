#!/usr/bin/env python

"""Module to TEMP TEMP.

Author: Raul A. Flores
"""

#| - IMPORT MODULES
import os
import sys

import time

import numpy as np
import pandas as pd

from plotly import io as pyio
import plotly.graph_objs as go
#__|


class ALAnalysis:
    """
    """

    #| - ALAnalysis ***********************************************************
    _TEMP = "TEMP"


    def __init__(self,
        ALBulkOpt=None,
        ):
        """
        """
        #| - __init__

        #| - Setting Argument Instance Attributes
        self.ALBulkOpt = ALBulkOpt
        #__|

        #| - Initializing Internal Instance Attributes

        #__|

        #__|

    #__| **********************************************************************


class ALAnimation:
    """
    """

    #| - ALAnimation **********************************************************
    _TEMP = "TEMP"


    def __init__(self,
        ALBulkOpt=None,
        duration_long=1000 * 6,
        duration_short=800 * 6,
        verbose=True,
        ):
        """
        """
        #| - __init__

        #| - Setting Argument Instance Attributes
        self.ALBulkOpt = ALBulkOpt
        self.duration_long = duration_long
        self.duration_short = duration_short
        self.verbose = verbose
        #__|

        #| - Initializing Internal Instance Attributes
        self.traces = dict()
        #__|

        verbose = self.verbose
        #__|


    def create_animation(self,
        duration_long=1000 * 6,
        duration_short=800 * 6,
        serial_parallel="parallel",  # 'serial' or 'parallel'

        marker_color_dict=None,
        ):
        """
        """
        #| - create_animation

        #| - Attributes #######################################################
        ALBulkOpt = self.ALBulkOpt
        verbose = self.verbose

        get_trace_j = self.get_trace_j
        get_layout = self.get_layout
        get_sliders_init_dict = self.get_sliders_init_dict
        get_slider_step_i = self.get_slider_step_i
        __save_figure_to_file__ = self.__save_figure_to_file__
        # __| #################################################################

        if verbose:
            print("\n", "Creating animation...")

        # #####################################################################
        self.__create_traces__(
            marker_color_dict=marker_color_dict,
            serial_parallel=serial_parallel)

        self.__create_figure__(
            duration_long=duration_long,
            duration_short=duration_short)

        # Save figure (HTML) to file
        __save_figure_to_file__()

        if verbose:
            print("DONE!")
        # __|

    def __create_traces__(self,
        marker_color_dict=None,
        serial_parallel=None
        ):
        """
        """
        #| - __create_traces__

        #| - Attributes #######################################################
        ALBulkOpt = self.ALBulkOpt
        verbose = self.verbose

        get_trace_j = self.get_trace_j
        get_layout = self.get_layout
        get_sliders_init_dict = self.get_sliders_init_dict
        get_slider_step_i = self.get_slider_step_i
        __save_figure_to_file__ = self.__save_figure_to_file__

        from multiprocessing import Pool
        from functools import partial
        # __| #################################################################


        # Shared kwargs for 'get_trace_j' method
        get_trace_j_kwargs = dict(
            prediction_key="y",
            uncertainty_key="err",
            plot_dft_instead_of_pred=True,
            trace_all_dft=True,
            trace_horiz_lines=True,
            marker_color_dict=marker_color_dict,
            # marker_size=8,
            )

        ALBulkOpt.al_gen_dict.items()

        if serial_parallel == "parallel":
            # | - Parallel execution
            t0 = time.time()

            # models = [i.model for i in ALBulkOpt.al_gen_dict.values()]
            AL_i_list = [i for i in ALBulkOpt.al_gen_dict.values()]
            traces_all = Pool().map(partial(
                get_trace_j,  # METHOD

                # KWARGS
                **get_trace_j_kwargs,
                ), AL_i_list)
                # ), models)

            run_time = time.time() - t0
            run_time_per = run_time / len(ALBulkOpt.al_gen_dict)

            self._create_traces__run_time = run_time
            self._create_traces__run_time_per = run_time_per

            iterator = enumerate(ALBulkOpt.al_gen_dict.items())
            for i_cnt, (gen_i, AL_gen_i) in iterator:
                traces_i = traces_all[i_cnt]
                traces_dict_i = {gen_i: traces_i}
                self.traces.update(traces_dict_i)
            # __|

        elif serial_parallel == "serial":
            #| - Serial execution
            t0 = time.time()

            iterator = enumerate(ALBulkOpt.al_gen_dict.items())
            for i_cnt, (gen_i, AL_gen_i) in iterator:
                if verbose:
                    print(4 * "", "Processing generation #", gen_i)

                # model_i = AL_gen_i.model

                # #############################################################
                traces_i = self.traces.get(gen_i, None)
                if traces_i is None:
                    traces_i = get_trace_j(
                        AL_gen_i,
                        # model_i,

                        # KWARGS
                        **get_trace_j_kwargs,
                        )

                    traces_dict_i = {gen_i: traces_i}
                    self.traces.update(traces_dict_i)

            run_time = time.time() - t0
            run_time_per = run_time / len(ALBulkOpt.al_gen_dict)

            self._create_traces__run_time = run_time
            self._create_traces__run_time_per = run_time_per
            # __|

        # __|

    def __create_figure__(self,
        duration_long=None,
        duration_short=None,
        ):
        """
        """
        #| - __create_figure__
        # #####################################################################
        ALBulkOpt = self.ALBulkOpt
        verbose = self.verbose

        traces_dict = self.traces
        get_trace_j = self.get_trace_j
        get_layout = self.get_layout
        get_sliders_init_dict = self.get_sliders_init_dict
        get_slider_step_i = self.get_slider_step_i
        __save_figure_to_file__ = self.__save_figure_to_file__
        # #####################################################################

        layout_anim = get_layout(
            duration_long=duration_long,
            duration_short=duration_short)
        sliders_dict = get_sliders_init_dict(duration_short)

        frames = []; data = []
        iterator = enumerate(ALBulkOpt.al_gen_dict.items())
        for i_cnt, (gen_i, AL_gen_i) in iterator:
            traces_i = traces_dict[gen_i]

            # #################################################################
            if i_cnt == 0: data.extend(traces_i)
            data_i = []; data_i.extend(traces_i)
            frame_i = go.Frame(data=data_i, name=str(i_cnt))
            frames.append(frame_i)
            slider_step_i = get_slider_step_i(i_cnt, duration_short)
            sliders_dict['steps'].append(slider_step_i)

        layout_anim["showlegend"] = True

        fig = go.Figure(
            data=data,
            layout=layout_anim,
            frames=frames)
        fig['layout']['sliders'] = [sliders_dict]

        self.fig = fig
        # __|

    def __save_figure_to_file__(self):
        """Save figure to file."""
        #| - __save_figure_to_file__
        # #####################################################################
        # duration_long = self.duration_long
        # duration_short = self.duration_short
        # get_trace_j = self.get_trace_j
        # get_layout = self.get_layout
        # get_sliders_init_dict = self.get_sliders_init_dict
        # get_slider_step_i = self.get_slider_step_i

        ALBulkOpt = self.ALBulkOpt
        fig = self.fig
        # #####################################################################


        file_path_i = os.path.join(
            "out_plot", "al_anim_" + ALBulkOpt.name + ".html")

        # print(file_path_i)
        pyio.write_html(
            fig,
            file_path_i)
        # __|

    def get_trace_j(self,
        AL_i,
        # model,
        prediction_key="prediction",
        uncertainty_key="uncertainty",
        plot_dft_instead_of_pred=True,
        trace_all_dft=True,
        trace_horiz_lines=True,
        marker_color_dict=None,
        ):
        """

        Args:
            plot_dft_instead_of_pred:
                Plot the actual DFT energy instead of the predicted value
        """
        #| - get_trace_j
        # #####################################################################
        ALBulkOpt = self.ALBulkOpt
        verbose = self.verbose

        traces_dict = self.traces
        get_trace_j = self.get_trace_j
        get_layout = self.get_layout
        get_sliders_init_dict = self.get_sliders_init_dict
        get_slider_step_i = self.get_slider_step_i
        __save_figure_to_file__ = self.__save_figure_to_file__
        # #####################################################################

        duplicate_swap_lists_all = ALBulkOpt.duplicate_swap_dict

        # id_color_dict

        # duplicate_swap_dict = AL.duplicate_swap_dict

        al_gen = AL_i.al_gen
        # al_gen = 1

        duplicate_swap_lists_prev = []
        for i in range(al_gen):
            # print(i)
            duplicate_swap_lists_prev += duplicate_swap_lists_all[i]

        model = AL_i.model


        # Must have consistent order of rows in DF for animation to work
        model = model.sort_index()

        # print("AL_i.al_gen:", AL_i.al_gen)
        # print("9yz2mt8hbh: ", model.index.get_loc("9yz2mt8hbh"))
        # print("6r716sxr9t: ", model.index.get_loc("6r716sxr9t"))


        #| - Rearrange indices to account for duplicates
        # #####################################################################
        # AL_i = AL.al_gen_dict[2]
        # model_i = AL_i.model
        model_i = model

        # model_i = model_i.sort_index()
        model_index = model_i.index.tolist()
        i = model_index

        duplicate_swap_lists = AL_i.duplicate_swap_lists

        # duplicate_swap_lists = duplicate_swap_lists + duplicate_swap_lists_all
        duplicate_swap_lists = duplicate_swap_lists + duplicate_swap_lists_prev
        if duplicate_swap_lists is None:
            duplicate_swap_lists = []


        for dupl_swap_i in duplicate_swap_lists:
            # a = model_index.index('9yz2mt8hbh')
            # b = model_index.index('6r716sxr9t')

            a = model_index.index(dupl_swap_i[0])
            b = model_index.index(dupl_swap_i[1])
            # print(a)
            # print(b)
            # print("---")

            model_index[b], model_index[a] = model_index[a], model_index[b]

            # a = model_index.index('9yz2mt8hbh')
            # b = model_index.index('6r716sxr9t')
            # print(a)
            # print(b)
            # print("---")


            # TEMP
            # print("ifjsjfisdifjkds")
            # print(dupl_swap_i)
            # print(marker_color_dict.keys())
            # print("ifjsjfisdifjkds")

            if dupl_swap_i[0] in marker_color_dict.keys():
                # print(40 * "TEMP")
                # print(marker_color_dict)
                print(marker_color_dict)
                marker_color_dict[dupl_swap_i[1]] = \
                    marker_color_dict.pop(dupl_swap_i[0])
                print(marker_color_dict)
                print(40 * "TEMP")

        model_i = model_i.reindex(
            labels=model_index,
            # index=
            )
        model = model_i

            # print(a)
            # print(b)
        # #####################################################################
        # __|


        # print("1", 4 * "*****----- ")
        # print("9yz2mt8hbh: ", model.index.get_loc("9yz2mt8hbh"))
        # print("6r716sxr9t: ", model.index.get_loc("6r716sxr9t"))




        # TEMP
        verbose = False
        ti = time.time()

        data = []

        # #####################################################################
        #| - If computed use DFT energy
        def method(row_i):
            #| - method
            computed = row_i["acquired"]
            y_real = row_i["y_real"]

            actually_computed = False
            if not np.isnan(y_real):
                actually_computed = True
            dft_energy = y_real

            predicted_energy = row_i[prediction_key]
            predicted_uncert = row_i[uncertainty_key]

            # #####################################################################
            new_column_values_dict = {
                "out_energy": None,
                "out_uncert": None}

            # #####################################################################
            if computed and actually_computed:
                new_column_values_dict["Y_main"] = dft_energy
                new_column_values_dict["Y_uncer"] = 0.
            else:
                new_column_values_dict["Y_main"] = predicted_energy
                new_column_values_dict["Y_uncer"] = predicted_uncert

            # #####################################################################
            for key, value in new_column_values_dict.items():
                row_i[key] = value
            return(row_i)
            #__|

        # if plot_dft_instead_of_pred:
        #     model = model.apply(method, axis=1)
        # else:
        #     model["Y_main"] = model[prediction_key]
        #     model["Y_uncer"] = model[uncertainty_key]

        model = ALBulkOpt.add_main_Y_to_model(
            model,
            plot_dft_instead_of_pred=True,
                prediction_key=prediction_key,
            uncertainty_key=uncertainty_key,

            # prediction_key="y",
            # uncertainty_key="err",
            )

        #__|

        # print("2", 4 * "*****----- ")
        # print("9yz2mt8hbh: ", model.index.get_loc("9yz2mt8hbh"))
        # print("6r716sxr9t: ", model.index.get_loc("6r716sxr9t"))


        # #####################################################################
        #| - Applying formating to df
        def method(row_i,
            # marker_size,
            marker_color_dict,
            ):
            #| - method
            new_column_values_dict = {}

            id_i = row_i.name

            computed_bool = row_i["acquired"]
            if computed_bool:
                # new_column_values_dict["marker_size"] = 10
                new_column_values_dict["marker_size"] = 6
                new_column_values_dict["marker_line_color"] = "black"
                new_column_values_dict["marker_color"] = "red"
                new_column_values_dict["marker_line_size"] = 0.1
            else:
                new_column_values_dict["marker_size"] = 6
                new_column_values_dict["marker_line_color"] = "black"
                new_column_values_dict["marker_color"] = "grey"
                new_column_values_dict["marker_line_size"] = 0.1

            if marker_color_dict is not None:
                if id_i in marker_color_dict.keys():
                    new_column_values_dict["marker_size"] = 10
                    new_column_values_dict["marker_line_color"] = \
                        marker_color_dict.get(id_i, "orange")
                    new_column_values_dict["marker_line_size"] = 4.5

            # #########################################################################
            for key, value in new_column_values_dict.items():
                row_i[key] = value
            return(row_i)
            #__|

        model = model.apply(method, axis=1,
            marker_color_dict=marker_color_dict,
            )
        #__|

        # print("3", 4 * "*****----- ")
        # print("9yz2mt8hbh: ", model.index.get_loc("9yz2mt8hbh"))
        # print("6r716sxr9t: ", model.index.get_loc("6r716sxr9t"))


        # #####################################################################
        #| - Removing duplicates from main trace (set them aside)

        # ████████ ███████ ███    ███ ██████
        #    ██    ██      ████  ████ ██   ██
        #    ██    █████   ██ ████ ██ ██████
        #    ██    ██      ██  ██  ██ ██
        #    ██    ███████ ██      ██ ██


        model_i_tmp = model[model["duplicate"] == False]
        # model_i_tmp["x_axis_ind_new"] = [i for i in range(model_i_tmp.shape[0])]
        model_i_tmp = model_i_tmp.sort_values("Y_main")
        model_i_tmp["x_axis_ind"] = [i for i in range(model_i_tmp.shape[0])]

        model = pd.concat([
            # model.drop("x_axis_ind", axis=1),
            model,
            model_i_tmp["x_axis_ind"]
            ],
            axis=1,
            sort=False,
            )
        model["x_axis_ind"] = model["x_axis_ind"].fillna("-60")


        # ████████ ███████ ███    ███ ██████
        #    ██    ██      ████  ████ ██   ██
        #    ██    █████   ██ ████ ██ ██████
        #    ██    ██      ██  ██  ██ ██
        #    ██    ███████ ██      ██ ██

        # __|


        # print("4", 4 * "*****----- ")
        # print("9yz2mt8hbh: ", model.index.get_loc("9yz2mt8hbh"))
        # print("6r716sxr9t: ", model.index.get_loc("6r716sxr9t"))
        # print("")

        # #####################################################################
        #| - Main data trace
        trace_i = go.Scatter(
            x=model["x_axis_ind"],
            # y=model[prediction_key],
            y=model["Y_main"],
            error_y=dict(
                type='data',
                array=model["Y_uncer"],
                visible=True,
                thickness=0.3,
                width=1.5,
                # color="rgba(120,120,120,1.0)",
                color="rgba(80,80,60,1.0)",
                ),
            # name=model["id"],
            mode="markers",

            text=model.index.tolist(),
            hoverinfo="text",

            marker={
                "opacity": 0.9,
                "size": model["marker_size"],
                "color": model["marker_color"],
                # "color": model[prediction_key],
                # "color": model["energy_pa"],
                # "colorscale": "Viridis",
                "line": {
                    "width": model["marker_line_size"],
                    "color": model["marker_line_color"],
                    },
                },
            )
        data.append(trace_i)
        #__|


        # #####################################################################
        #| - Horizontal lines at E minimum
        min_e = model[prediction_key].min()
        min_e_w_uncert = (model[prediction_key] - model[uncertainty_key]).min()

        # #########################################################################
        trace_i = go.Scatter(
            x=[0, len(model["x_axis_ind"].tolist())],
            y=[min_e, min_e],
            mode="lines",
            line=dict(
                color="firebrick",
                width=2,
                dash="dot",
                ),
            )

        if trace_horiz_lines:
            data.append(trace_i)

        # #########################################################################
        trace_i = go.Scatter(
            x=[0, len(model["x_axis_ind"].tolist())],
            y=[min_e_w_uncert, min_e_w_uncert],
            mode="lines",
            line=dict(
                color="black",
                width=1.5,
                dash="dash",
                ),
            )
        if trace_horiz_lines:
            data.append(trace_i)
        #__|


        # #####################################################################
        #| - Validation DFT Trace
        df_dft_val = model[~model["y_real"].isna()]
        df_dft_val = df_dft_val[["y_real", "acquired", "x_axis_ind"]]

        def method(row_i):
            #| - method
            acquired = row_i["acquired"]
            y_real = row_i["y_real"]

            actually_computed = False
            if not np.isnan(y_real): actually_computed = True

            # dft_energy = y_real

            # Computed and Acquired (HIDE)
            if acquired and actually_computed:
                row_i["marker_opacity"] = 0.
                row_i["marker_color"] = "black"
                row_i["marker_line_color"] = "pink"
                row_i["marker_size"] = 8.
                row_i["marker_symbol"] = "diamond"

                # | - __old__
                # row_i["marker_symbol"] = "cross-open"
                # row_i["marker_symbol"] = "asterisk-open"
                # row_i["marker_symbol"] = "diamond-open"
                # row_i["marker_symbol"] = "circle"
                # row_i["energy_pa"] = -4.6
                # if redundant_global:
                #     row_i["marker_opacity"] = 1.
                #     row_i["marker_symbol"] = "star"
                #     row_i["marker_color"] = "green"
                # __|

            # Validation data not available (HIDE)
            elif actually_computed == False:
                row_i["marker_opacity"] = 1.
                row_i["marker_color"] = "grey"
                row_i["marker_size"] = 0.
                row_i["marker_symbol"] = "diamond"

            # Main points (Actually show up)
            elif not acquired and actually_computed:
                row_i["marker_opacity"] = 1.
                row_i["marker_color"] = "grey"
                row_i["marker_line_color"] = "black"
                row_i["marker_size"] = 8.
                # row_i["marker_size"] = 4.
                row_i["marker_symbol"] = "diamond-open"

                # if redundant_global:
                #     row_i["marker_color"] = "green"
                #     row_i["marker_symbol"] = "star-open"

            else:
                print("NOT CAUGHT BEHAVIOR")

            return(row_i)
            #__|
        df_dft_val = df_dft_val.apply(method, axis=1)

        if verbose:
            tf = time.time()
            print("    ", "Validation DFT Trace | apply method:", tf - ti); ti = tf

        trace_validation = go.Scatter(
            y=df_dft_val["y_real"],
            x=df_dft_val["x_axis_ind"],
            mode="markers",
            marker=dict(
                symbol=df_dft_val["marker_symbol"],
                color=df_dft_val["marker_color"],
                opacity=df_dft_val["marker_opacity"],
                size=df_dft_val["marker_size"],
                line=dict(
                    color=df_dft_val["marker_line_color"],
                    width=1.,
                    ),
                ),
            )

        # data.append(trace_validation)

        # Add trace to beginning of list (Show below other traces)
        data.insert(0, trace_validation)


        if verbose:
            tf = time.time()
            print("    ", "Validation DFT Trace | create trace:", tf - ti); ti = tf
        # __|


        return(data)


        #| - __old__
        # if verbose:
        #     tf = time.time()
        #     print("If computed use DFT energy:", tf - ti); ti = tf

        # if verbose:
        #     tf = time.time()
        #     print("Applying formating to df:", tf - ti); ti = tf

        # if verbose:
        #     tf = time.time()
        #     print("Removing duplicates from main trace:", tf - ti); ti = tf

        # if verbose:
        #     tf = time.time()
        #     print("Main data trace:", tf - ti); ti = tf

        # if verbose:
        #     tf = time.time()
        #     print("Horizontal lines:", tf - ti); ti = tf

        # if verbose:
        #     tf = time.time()
        #     print("Validation DFT Trace:", tf - ti); ti = tf
        # __|

        #__|

    def get_layout(self,
        duration_long=1000 * 2,
        duration_short=800 * 2,
        ):
        """
        """
        #| - get_layout

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

    def get_sliders_init_dict(self, duration_short):
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

    def get_slider_step_i(self, i_cnt, duration_short):
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


    #__| **********************************************************************
