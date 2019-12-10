#!/usr/bin/env python

"""Module to TEMP TEMP.

Author: Raul A. Flores
"""

#| - IMPORT MODULES
import os
import sys

import copy
import time

from pathlib import Path

import numpy as np
import pandas as pd

from multiprocessing import Pool
from functools import partial

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
        # duration_long=1000 * 6,
        # duration_short=800 * 6,
        marker_color_dict=None,
        verbose=True,
        ):
        """
        """
        #| - __init__

        #| - Setting Argument Instance Attributes
        self.ALBulkOpt = ALBulkOpt
        # self.duration_long = duration_long
        # self.duration_short = duration_short
        self.marker_color_dict = marker_color_dict
        self.verbose = verbose
        #__|

        #| - Initializing Internal Instance Attributes
        self.traces = dict()

        self.swap_histories = None
        #__|


        # self.duplicate_system_history_analysis()

        verbose = self.verbose
        #__|


    def create_animation(self,
        duration_long=1000 * 6,
        duration_short=800 * 6,
        serial_parallel="parallel",  # 'serial' or 'parallel'
        filename=None,
        # marker_color_dict=None,
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
            # marker_color_dict=marker_color_dict,
            serial_parallel=serial_parallel,
            read_traces_from_file=True,
            )

        self.__create_figure__(
            duration_long=duration_long,
            duration_short=duration_short)

        # Save figure (HTML) to file
        __save_figure_to_file__(filename=filename)


        if verbose:
            print("DONE!")
        # __|

    def __create_traces__(self,
        # marker_color_dict=None,
        serial_parallel=None,
        read_traces_from_file=True,
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

        # __| #################################################################


        traces_read_succ = False

        # TEMP
        read_traces_from_file = False
        if read_traces_from_file:
            print("TEMP | 0 | ijs8hf8sdf")
            # #################################################################
            import pickle; import os
            path_i = os.path.join(
                # os.environ[""],
                "out_data",
                "ALAnim__traces.pickle")

            my_file = Path(path_i)
            if my_file.exists():
                print("TEMP | 1 | ijs8hf8sdf")

                with open(path_i, "rb") as fle:
                    traces = pickle.load(fle)
                    self.traces = traces
                    traces_read_succ = True
            # #################################################################

        if not traces_read_succ:
            #| - Create traces
            # Shared kwargs for 'get_trace_j' method
            get_trace_j_kwargs = dict(
                prediction_key="y",
                uncertainty_key="err",
                plot_dft_instead_of_pred=True,
                # trace_all_dft=False,
                trace_horiz_lines=False,
                plot_validation_dft=False,
                # marker_color_dict=marker_color_dict,
                # marker_size=8,
                )

            # ALBulkOpt.al_gen_dict.items()

            if serial_parallel == "parallel":
                # | - Parallel execution
                t0 = time.time()

                # models = [i.model for i in ALBulkOpt.al_gen_dict.values()]
                AL_i_list = [i for i in ALBulkOpt.al_gen_dict.values()]

                # TEMP
                # AL_i_list = AL_i_list[0:5]

                # TEMP_PRINT
                # print("IJSIFJISDJFSDF*SDF&DS")
                # print(get_trace_j_kwargs)

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
                    # TEMP
                    try:
                        traces_i = traces_all[i_cnt]
                        traces_dict_i = {gen_i: traces_i}
                        self.traces.update(traces_dict_i)
                    except:
                        print(i_cnt)
                        pass
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


            # Pickling data #######################################################
            import os; import pickle
            directory = "out_data"
            if not os.path.exists(directory): os.makedirs(directory)
            file_path = os.path.join(directory, "ALAnim__traces.pickle")
            with open(file_path, "wb") as fle:
                pickle.dump(self.traces, fle)
            # #####################################################################

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

        al_gen_dict = ALBulkOpt.al_gen_dict


        layout_anim = get_layout(
            duration_long=duration_long,
            duration_short=duration_short)
        sliders_dict = get_sliders_init_dict(duration_short)

        frames = []; data = []
        # iterator = enumerate(ALBulkOpt.al_gen_dict.items())
        # for i_cnt, (gen_i, AL_gen_i) in iterator:
        for i_cnt, (gen_i, traces_i) in enumerate(traces_dict.items()):
            # traces_i = traces_dict[gen_i]

            #| - Adding custom layout to frames (Work in progress)
            # AL_i = al_gen_dict[gen_i]
            # model = AL_i.model
            # num_dft_calcs = model[model["acquired"] == True].shape[0]
            # print(num_dft_calcs)
            # layout_i = go.Layout(
            #
            #     annotations=[
            #         go.layout.Annotation(
            #             x=0,
            #             y=1,
            #             xref="paper",
            #             yref="paper",
            #             text="DFT: " + str(num_dft_calcs).zfill(3),
            #             showarrow=False,
            #             # arrowhead=7,
            #             # ax=0,
            #             # ay=-40
            #             )
            #         ],
            #
            #     )
            # __|

            # #################################################################
            if i_cnt == 0: data.extend(traces_i)
            data_i = []; data_i.extend(traces_i)

            frame_i = go.Frame(data=data_i, name=str(i_cnt),
                # layout=layout_i,
                )

            frames.append(frame_i)
            slider_step_i = get_slider_step_i(i_cnt, duration_short)
            sliders_dict['steps'].append(slider_step_i)


        # layout_anim["showlegend"] = True

        fig = go.Figure(
            data=data,
            layout=layout_anim,
            frames=frames)
        fig['layout']['sliders'] = [sliders_dict]

        self.fig = fig
        # __|

    def __save_figure_to_file__(self, filename=None):
        """Save figure to file."""
        #| - __save_figure_to_file__
        # #####################################################################
        ALBulkOpt = self.ALBulkOpt
        fig = self.fig
        # #####################################################################


        # fig.layout["showlegend"] = False

        if filename is None:
            file_path_i = os.path.join(
                "out_plot", "al_anim_" + ALBulkOpt.name + ".html")
        else:

            if ".html" not in filename:
                filename = filename + ".html"

            file_path_i = os.path.join(
                "out_plot",
                filename,
                )


        print(file_path_i)
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
        # trace_all_dft=True,
        trace_horiz_lines=True,
        plot_validation_dft=True,
        internally_order_df=False,
        dft_calc_al_gen_text_overlay=True,
        # marker_color_dict=None,
        ):
        """

        Args:
            plot_dft_instead_of_pred:
                Plot the actual DFT energy instead of the predicted value
            internally_order_df:
                If True, order the rows of the model dataframe by whether
                they've been acquired and then by the target value
                NOTE: This messes up the animation feature and is only intended
                for generating single generation figures
        """
        #| - get_trace_j
        # #####################################################################
        ALBulkOpt = self.ALBulkOpt
        verbose = self.verbose
        marker_color_dict = self.marker_color_dict
        traces_dict = self.traces
        get_trace_j = self.get_trace_j
        get_layout = self.get_layout
        get_sliders_init_dict = self.get_sliders_init_dict
        get_slider_step_i = self.get_slider_step_i
        __save_figure_to_file__ = self.__save_figure_to_file__
        # #####################################################################

        # print("plot_validation_dft:", plot_validation_dft)

        # TEMP
        verbose = False
        ti = time.time()

        duplicate_swap_lists_all = ALBulkOpt.duplicate_swap_dict

        al_gen = AL_i.al_gen
        print("al_gen:", al_gen)

        if marker_color_dict is not None:
            marker_color_dict = self.__get_color_dict__(
                gen_i=al_gen,
                )

        duplicate_swap_lists_prev = []
        for i in range(al_gen):
            duplicate_swap_lists_prev += duplicate_swap_lists_all[i]

        model = AL_i.model

        # Must have consistent order of rows in DF for animation to work
        model = model.sort_index()


        #| - Rearrange indices to account for duplicates
        model_i = model

        model_index = model_i.index.tolist()
        i = model_index

        duplicate_swap_lists = AL_i.duplicate_swap_lists

        duplicate_swap_lists = duplicate_swap_lists_prev + duplicate_swap_lists
        if duplicate_swap_lists is None:
            duplicate_swap_lists = []


        for dupl_swap_i in duplicate_swap_lists:
            a = model_index.index(dupl_swap_i[0])
            b = model_index.index(dupl_swap_i[1])

            model_index[b], model_index[a] = model_index[a], model_index[b]

        # self.marker_color_dict = marker_color_dict

        model_i = model_i.reindex(labels=model_index)
        model = model_i
        # __|


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


        model = ALBulkOpt.add_main_Y_to_model(
            model,
            plot_dft_instead_of_pred=True,
            prediction_key=prediction_key,
            uncertainty_key=uncertainty_key,
            )

        #__|


        # #####################################################################
        #| - Applying formating to df
        def method(row_i,
            marker_color_dict,
            ):
            #| - method
            new_column_values_dict = {}

            id_i = row_i.name

            computed_bool = row_i["acquired"]
            if computed_bool:
                # new_column_values_dict["marker_size"] = 10

                new_column_values_dict["marker_symbol"] = "circle"
                new_column_values_dict["marker_size"] = 6
                new_column_values_dict["marker_line_color"] = "black"
                new_column_values_dict["marker_color"] = "rgba(255,0,0,0.6)"
                new_column_values_dict["marker_line_size"] = 0.1
            else:
                new_column_values_dict["marker_symbol"] = "circle-open"
                # new_column_values_dict["marker_symbol"] = "circle"
                new_column_values_dict["marker_size"] = 6
                new_column_values_dict["marker_line_color"] = "black"
                new_column_values_dict["marker_color"] = "grey"
                new_column_values_dict["marker_line_size"] = 0.2


            # if marker_color_dict is not None:
            if False:
                if id_i in marker_color_dict.keys():
                    # new_column_values_dict["marker_symbol"] = "circle"
                    # new_column_values_dict["marker_symbol"] = "circle-cross"
                    new_column_values_dict["marker_symbol"] = "x"

                    # new_column_values_dict["marker_size"] = 6
                    new_column_values_dict["marker_size"] = 8
                    # new_column_values_dict["marker_line_color"] = \
                    #     marker_color_dict.get(id_i, "orange")

                    new_column_values_dict["marker_line_color"] = "green"
                    new_column_values_dict["marker_line_size"] = 2.0

                    if computed_bool:
                        new_column_values_dict["marker_color"] = \
                            "rgba(255,0,0,1.0)"
                    else:
                        new_column_values_dict["marker_color"] = \
                            "white"


            # #########################################################################
            for key, value in new_column_values_dict.items():
                row_i[key] = value
            return(row_i)
            #__|

        model = model.apply(method, axis=1,
            marker_color_dict=marker_color_dict,
            )
        #__|


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
        # TEMP Why is this a string of -60?
        # model["x_axis_ind"] = model["x_axis_ind"].fillna("-60")
        model["x_axis_ind"] = model["x_axis_ind"].fillna(-60)


        # ████████ ███████ ███    ███ ██████
        #    ██    ██      ████  ████ ██   ██
        #    ██    █████   ██ ████ ██ ██████
        #    ██    ██      ██  ██  ██ ██
        #    ██    ███████ ██      ██ ██

        # __|


        if internally_order_df:
            model_acq_t = model[model["acquired"] == True]
            model_acq_f = model[model["acquired"] == False]

            key_tmp = "Y_main"
            model = pd.concat([
                model_acq_f.sort_values(key_tmp),
                model_acq_t.sort_values(key_tmp)])

            # if id_i in marker_color_dict.keys():

            # Bringing 'special' systems up front
            # special_ids = list(marker_color_dict.keys())
            # special_ids = \
            #     [i for i in model.index.tolist() if i in special_ids]
            #
            # model_0 = model.loc[special_ids]
            # # model_1 = model.loc[~special_ids]
            #
            # model_1 = model[~model.index.isin(special_ids)]
            #
            # model = pd.concat([
            #     model_1.sort_values(key_tmp),
            #     model_0.sort_values(key_tmp),
            #     ])



        special_ids = list(marker_color_dict.keys())
        special_ids = \
            [i for i in model.index.tolist() if i in special_ids]

        model_0 = model.loc[special_ids]
        # print(3 * "\n")
        # print(model_0)

        # TEMP
        # Pickling data ######################################################
        import os; import pickle
        # directory = "out_data"
        # if not os.path.exists(directory): os.makedirs(directory)
        # with open(os.path.join(os.environ["HOME"], "__temp__", "TEMP.pickle"), "wb") as fle:
        #     pickle.dump(model_0, fle)
        # #####################################################################


        # Adding vertical traces to track top 10
        for i_ind, row_i in model_0.iterrows():
            tmp1 = 42

            Y_main = row_i["Y_main"]
            x_ind = row_i["x_axis_ind"]
            acquired = row_i["acquired"]

            if acquired:
                color = "red"
                y = [4.7, 6]
                width = 0.8
            else:
                color = "grey"
                y = [5., 6]
                width = 0.5

            trace_i = go.Scatter(
                mode="lines",
                x=[x_ind, x_ind],
                y=y,
                line=dict(
                    width=width,
                    color=color,
                    ),
                )

            data.append(trace_i)


        # trace_i = go.Scatter(
        #     mode="lines",
        #     x=[100, 100],
        #     y=[0, 4],
        #     )
        #
        #
        # data.append(trace_i)





        # #####################################################################
        #| - Main data trace


        #| - Error filled area trace
        shared_scatter_props = dict(
            mode="lines",

            # marker=dict(
            #     symbol="circle",
            #     size=6,
            #     opacity=0.5,
            #     line=dict(
            #         color='black',
            #         width=1,
            #         )
            #     ),

            line=dict(
                width=0.,
                # color="grey",
                color="rgba(150,150,150,.8)",
                # dash="dash",
                ),

            # error_y=dict(
            #     thickness=.5,
            #     width=1.0,
            #     ),
            )


        # TEMP
        # print(model["x_axis_ind"].tolist())
        model_tmp = model.sort_values("x_axis_ind")

        trace = go.Scatter(
            x=model_tmp["x_axis_ind"],
            y=model_tmp["Y_main"] + model_tmp["Y_uncer"],
            **shared_scatter_props,
            # fill="tonexty",
            )
        data.append(trace)

        trace = go.Scatter(
            x=model_tmp["x_axis_ind"],
            y=model_tmp["Y_main"] - model_tmp["Y_uncer"],
            fill="tonexty",
            **shared_scatter_props,
            )
        data.append(trace)
        # __|


        trace_i = go.Scatter(
            x=model["x_axis_ind"],
            # y=model[prediction_key],
            y=model["Y_main"],

            # error_y=dict(
            #     type='data',
            #     array=model["Y_uncer"],
            #     visible=True,
            #     thickness=0.3,
            #     width=1.5,
            #     # color="rgba(120,120,120,1.0)",
            #     color="rgba(80,80,60,1.0)",
            #     ),
            # name=model["id"],
            mode="markers",

            text=model.index.tolist(),
            hoverinfo="text",

            marker={
                "opacity": 1.,
                "symbol": model["marker_symbol"],
                "size": model["marker_size"],
                "color": model["marker_color"],
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
        if plot_validation_dft:
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


        # #####################################################################
        #| - DFT Calc and AL Gen Num Text Overlay
        # dft_calc_al_gen_text_overlay = True
        if dft_calc_al_gen_text_overlay:
            num_dft_calcs = model[model["acquired"] == True].shape[0]

            trace_i = go.Scatter(
                x=[0],
                y=[model["y_real"].max()],
                mode="text",
                # name="Lines, Markers and Text",
                text=["DFT Calcs: " + str(num_dft_calcs).zfill(3)],

                textposition="middle right",
                textfont=dict(
                    family="Arial",
                    size=20,
                    color="black"
                    )
                )
            data.append(trace_i)



            trace_i = go.Scatter(
                x=[0],
                y=[model["y_real"].max() - 1],
                mode="text",
                text=["AL gen: " + str(al_gen).zfill(3)],

                textposition="middle right",
                textfont=dict(
                    family="Arial",
                    size=20,
                    color="black"
                    )
                )
            data.append(trace_i)
        # __|


        return(data)
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
                family='Arial',
                size=20,
                color='black'
                ),

            xaxis={
                'title': 'Candidate Space',

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
                'title': 'ΔH<sub>f</sub> (eV)',

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

            paper_bgcolor="white",
            plot_bgcolor="white",

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

    def duplicate_system_history_analysis(self):
        """
        """
        #| - duplicate_system_history_analysis
        # #####################################################################
        ALBulkOpt = self.ALBulkOpt
        verbose = self.verbose
        marker_color_dict = self.marker_color_dict
        traces_dict = self.traces
        get_trace_j = self.get_trace_j
        get_layout = self.get_layout
        get_sliders_init_dict = self.get_sliders_init_dict
        get_slider_step_i = self.get_slider_step_i
        __save_figure_to_file__ = self.__save_figure_to_file__
        # #####################################################################

        self.__create_swap_histories__()

        self.__color_dict_progression__()
        # __|

    def __get_color_dict__(self,
        # id_color_dict=None,
        gen_i=None,
        # color_dict_progression=None,
        ):
        """
        """
        #| - get_color_dict
        # #####################################################################
        ALBulkOpt = self.ALBulkOpt
        verbose = self.verbose
        marker_color_dict = self.marker_color_dict
        traces_dict = self.traces
        get_trace_j = self.get_trace_j
        get_layout = self.get_layout
        get_sliders_init_dict = self.get_sliders_init_dict
        get_slider_step_i = self.get_slider_step_i

        swap_histories = self.swap_histories
        # color_dict_progression = self.color_dict_progression
        # self.swap_histories = swap_histories

        __save_figure_to_file__ = self.__save_figure_to_file__
        # #####################################################################

        color_dict_progression = ALBulkOpt.color_dict_progression

        id_color_dict = marker_color_dict

        id_color_dict_gen_i = dict()
        for id_i, color_i in id_color_dict.items():
            if id_i in color_dict_progression.keys():
                tmp1 = color_dict_progression.get(id_i).get(gen_i, None)
                id_color_dict_gen_i[tmp1] = color_i

            else:
                id_color_dict_gen_i[id_i] = color_i

        return(id_color_dict_gen_i)
        #__|

    #__| **********************************************************************



class ALPerformance:
    """
    """

    #| - ALPerformance ********************************************************
    _TEMP = "TEMP"


    def __init__(self,
        ALBulkOpt=None,
        verbose=False,
        ):
        """
        """
        #| - __init__

        #| - Setting Argument Instance Attributes
        self.ALBulkOpt = ALBulkOpt
        self.verbose = verbose
        #__|


        #| - Initializing Internal Instance Attributes

        #__|

        #__|

    def num_sys_discovered(self,
        mode="perc",  # 'perc', 'num', or 'user_specified'
        perc_of_structs=10,
        num_structs=None,
        ids_to_track=None,
        account_duplicates=True,
        ):
        """
        Args:
            mode: 'perc' or 'num'
                If 'perc', track the best Nth percent structures
                If 'num' track the best N structures
                If 'user_specified', then user given list is used as ids to
              track
            perc_of_structs:
                Percent of total number of structures to track
            num_structs;
                Best N structures to track
            ids_to_track:
                User defined list of top ids to track

        """
        #| - num_sys_discovered

        # #####################################################################
        ALBulkOpt = self.ALBulkOpt
        verbose = self.verbose
        # #####################################################################
        al_gen_dict = ALBulkOpt.al_gen_dict
        # #####################################################################

        swap_histories = ALBulkOpt.swap_histories
        last_gen = list(al_gen_dict.keys())[-1]
        AL_last = al_gen_dict[last_gen]
        model = AL_last.model

        # #####################################################################
        # #####################################################################
        # #####################################################################

        #| - Gettings ids to track
        if mode != "user_specified":
            if mode == "perc":
                num_candidates_init = model.shape[0]
                num_track_structs = \
                    round(num_candidates_init * (perc_of_structs * 0.01))
            elif mode == "num":
                if num_structs is not None:
                    num_track_structs = num_structs
                else:
                    assert False, "Must define 'num_structs' when mode='num'"

            model_tmp = model[model["duplicate"] == False]
            model_tmp = model_tmp.sort_values("y_real")
            model_tmp = model_tmp.iloc[0:num_track_structs]
            top_ids = model_tmp.index.tolist()

        if mode == "user_specified":
            top_ids = ids_to_track

        top_ids_static = copy.deepcopy(top_ids)
        top_ids_working = copy.deepcopy(top_ids)
        # __|

        # print("len(top_ids_working):", len(top_ids_working))

        # #####################################################################
        # #####################################################################
        # #####################################################################
        duplicates_of_top_ids = []
        for id_i in top_ids_working:
            if id_i in swap_histories.keys():
                swap_lists = swap_histories.get(id_i, "TEMP")

                swap_ids_i = []
                for gen_j, swap_list_j in swap_lists.items():
                    swap_ids_i.extend(swap_list_j)

                duplicates_of_top_ids.extend(swap_ids_i)

        duplicates_of_top_ids = list(set(duplicates_of_top_ids))

        # #####################################################################
        # #####################################################################
        # #####################################################################
        top_ids_w_dupl = list(set(duplicates_of_top_ids + top_ids_working))

        # #####################################################################
        # #####################################################################
        # #####################################################################
        new_swap_dict = dict()
        for id_i, swap_history_i in swap_histories.items():
            for gen_j, swap_list_j in swap_history_i.items():
                for swap_id in swap_list_j:
                    # #########################################################
                    if swap_id in new_swap_dict.keys():
                        if new_swap_dict[swap_id] != id_i:
                            print("This id corresponds to more than 1 final id")

                    new_swap_dict[swap_id] = id_i


        # #####################################################################
        # #####################################################################
        # #####################################################################

        data_list_master = []
        for gen_i, AL_i in al_gen_dict.items():
            data_dict_i = dict()

            model_i = AL_i.model
            model_tmp = model_i[
                (model_i["acquired"] == True) & \
                (model_i["duplicate"] == False)
                ]

            # Number of DFT experiments
            num_dft_calcs = model_i[model_i["acquired"] == True].shape[0]
            data_dict_i["num_dft"] = num_dft_calcs

            for id_i in model_tmp.index:

                if id_i in top_ids_working:
                    top_ids_working.remove(id_i)

                if account_duplicates:
                    final_swap_id = new_swap_dict.get(id_i, None)
                    if final_swap_id is not None and final_swap_id in top_ids_working:
                        top_ids_working.remove(final_swap_id)

            num_ids_disc = len(top_ids_static) - len(top_ids_working)
            data_dict_i["num_ids_discovered"] = num_ids_disc

            data_list_master.append(data_dict_i)


        df = pd.DataFrame(data_list_master)

        self.num_sys_discovered_df = df
        # __|

    #__| **********************************************************************
