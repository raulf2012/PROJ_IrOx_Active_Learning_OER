# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
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

# #########################################################
from active_learning.al_analysis import ALAnimation
# -

stoich_i = "AB2"
# stoich_i = "AB3"

# +
if stoich_i == "AB3":
    path_i = main_AB3_run
elif stoich_i == "AB2":
    path_i = main_AB2_run

# #########################################################
with open(path_i, "rb") as fle:
    AL = pickle.load(fle)

# + active=""
#
#
#
#
#
#

# +
al_gen_dict = AL.al_gen_dict
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

ALAnim = ALAnimation(
    ALBulkOpt=AL,
    marker_color_dict=marker_color_dict,
    verbose=True,
    color_custom_points=True,
    # gens_to_plot = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 22],
    )

# +
# duration_long=1000 * 6
# duration_short=800 * 6

duration_long=1000 * 2
duration_short=800 * 2

serial_parallel="parallel"  # 'serial' or 'parallel'
# serial_parallel="serial"  # 'serial' or 'parallel'
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
    dft_calc_al_gen_text_overlay=False,
    # marker_color_dict=marker_color_dict,
    # marker_size=8,
    add_vertical_track_lines=True,
    just_traces=True,
    error_type="normal",  # filled or normal
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

ALAnim.get_global_props()

# +
# ALAnim.y_max_global

AL = self.ALBulkOpt

AL_i = AL.al_gen_dict[0]

num_xaxis = AL_i.model.shape[0]
# -

import numpy as np

# +
y_max_list = []
for gen_i, AL_i in AL.al_gen_dict.items():
    model = AL_i.model

    model_tmp = model[model.duplicate == False]

    y_max = model_tmp.y.max()
    y_max_list.append(y_max)
    
y_max_global = np.max(y_max_list)

# +
AL_i = AL.al_gen_dict[3]

model = AL_i.model

# + active=""
#

# +
# model_tmp = model[model.duplicate == False]

# y_min = model_tmp.y_real.min()
# y_max = model_tmp.y_real.max()
