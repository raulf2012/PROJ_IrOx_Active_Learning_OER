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

# +
import os
print(os.getcwd())
import sys

import pickle

# #########################################################
# Local Imports
sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import tmp

from al_data import main_AB2_run, main_AB3_run
# -

# stoich_i = "AB2"
stoich_i = "AB3"

# +
if stoich_i == "AB3":
    path_i = main_AB3_run
elif stoich_i == "AB2":
    path_i = main_AB2_run

# #############################################################################
with open(path_i, "rb") as fle:
    AL = pickle.load(fle)

# + active=""
#
#
#
#
#
#
# -

al_gen_dict = AL.al_gen_dict

# +
last_gen_key = list(al_gen_dict.keys())[-1]

# if gens_to_plot[-1] == "last":
#     gen_4 = last_gen_key
#     gens_to_plot[-1] = gen_4

lowest_N_sys_to_track = 10

AL_last = al_gen_dict[last_gen_key]

model = AL_last.model

model_i = model[
    (model["duplicate"] == False) & \
    (model["acquired"] == True)
    ].sort_values("y_real")
top_ids_to_track = model_i.iloc[0:lowest_N_sys_to_track].index.tolist()


color_list = [

    # "rgba(12,0,127,1.0)",
    # "rgba(0,0,172,1.0)",
    # "rgba(0,1,215,1.0)",
    # "rgba(0,51,233,1.0)",
    # "rgba(0,83,255,1.0)",
    # "rgba(0,115,255,1.0)",
    # "rgba(0,141,243,1.0)",
    # "rgba(0,181,246,1.0)",
    # "rgba(0,220,245,1.0)",
    # "rgba(0,255,243,1.0)",

    "#fde725",
    "#b8de29",
    "#74d055",
    "#3cbc75",
    "#20a386",
    "#238a8d",
    "#2d708e",
    "#39558c",
    "#453781",
    "#481568",

    # "rgb(202,88,66)",
    # "rgb(71,189,198)",
    # "rgb(210,70,147)",
    # "rgb(120,181,66)",
    # "rgb(157,99,201)",
    # "rgb(81,163,108)",
    # "rgb(189,104,138)",
    # "rgb(131,128,57)",
    # "rgb(101,130,203)",
    # "rgb(209,154,68)",
    ]


marker_color_dict = dict(zip(
    top_ids_to_track,
    color_list,
    ))

# + active=""
#
#
#
#
#
#
# -

from active_learning.al_analysis import ALAnimation

ALAnim = ALAnimation(
    ALBulkOpt=AL,
    marker_color_dict=marker_color_dict,
    verbose=True,
    color_custom_points=True,
    )

# +
duration_long=1000 * 6
duration_short=800 * 6
serial_parallel="parallel"  # 'serial' or 'parallel'
filename=None


# | - Attributes #######################################################
ALBulkOpt = ALAnim.ALBulkOpt
verbose = ALAnim.verbose

get_trace_j = ALAnim.get_trace_j
get_layout = ALAnim.get_layout
get_sliders_init_dict = ALAnim.get_sliders_init_dict
get_slider_step_i = ALAnim.get_slider_step_i
__save_figure_to_file__ = ALAnim.__save_figure_to_file__
#__| #################################################################

if verbose:
    print("\n", "Creating animation...")

# #####################################################################
get_trace_j_kwargs = dict(
    prediction_key="y",
    uncertainty_key="err",
    plot_dft_instead_of_pred=True,
    # trace_all_dft=False,
    trace_horiz_lines=False,
    plot_validation_dft=False,
    # marker_color_dict=marker_color_dict,
    # marker_size=8,
    add_vertical_track_lines=True,
    just_traces=True,
    )


ALAnim.__create_traces__(
    # marker_color_dict=marker_color_dict,
    serial_parallel=serial_parallel,
    read_traces_from_file=True,
    get_trace_j_kwargs=get_trace_j_kwargs,
    )

ALAnim.__create_figure__(
    duration_long=duration_long,
    duration_short=duration_short)

# Save figure (HTML) to file
__save_figure_to_file__(filename=filename)


if verbose:
    print("DONE!")
# -

assert False

# +
# ALAnim.create_animation?

# + active=""
#
#
#
#

# + jupyter={}
# ALAnim = ALAnimation(
#     ALBulkOpt=AL,
#     # marker_color_dict=id_color_dict,
#     verbose=True)

# ALAnim.create_animation(
#     duration_long=1000 * 4,
#     duration_short=800 * 4,
#     serial_parallel="parallel",  # 'serial' or 'parallel'
# #     marker_color_dict=id_color_dict,
#     )

# + jupyter={}
# sys.path.insert(0,
#     os.path.join(
#         os.environ["PROJ_irox"],
#         "python_classes"))

# sys.path.insert(0,
#     os.path.join(
#         os.environ["PROJ_irox"],
#         "python_classes/active_learning"))

# sys.path.insert(0,
#     os.path.join(
#         os.environ["PROJ_irox"],
#         "python_classes/active_learning/"))

# from active_learning.al_algeneration import ALGeneration
# from active_learning.al_bulkopt import ALBulkOpt
# from active_learning import al_bulkopt
