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

# + [markdown] Collapsed="false"
# # Import Modules

# +
import os
print(os.getcwd())
import sys

import pickle

import numpy as np
import pandas as pd

# #########################################################
sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import tmp
from al_data import al_data_files_dict

# + [markdown] Collapsed="false"
# # Script Inputs

# +
stoich_i = "AB2"
# stoich_i = "AB3"

# gp_or_random = "gp"
gp_or_random = "random"


if gp_or_random == "gp":
    files_list = al_data_files_dict[stoich_i]["files_list_gp_ucb"]
elif gp_or_random == "random":
    files_list = al_data_files_dict[stoich_i]["files_list_random"]
# -

files_list[0:4]

# + Collapsed="false"
out_data_dict = dict()

data_dict = dict()
for file_i in files_list:
    # #########################################################################
    num = file_i.split("_")[-1].split(".")[0]
    with open(file_i, "rb") as fle:
        AL_i = pickle.load(fle)

    data_dict[num] = AL_i
out_data_dict["AL_dict"] = data_dict


# -

def get_gen_acquired(
    AL_i,
    id_of_most_stable
    # init_id_of_most_stable,
    ):
    out_dict = dict()
    out_dict["first_acq_id"] = None

    color_progression_i = AL_i.color_dict_progression.get(id_of_most_stable, None)

    gen_acquired = None

    if color_progression_i is not None:
        init_id_of_most_stable = color_progression_i[0]

        for gen_i, AL_gen_i in AL_i.al_gen_dict.items():
            # print("init_id_of_most_stable:", init_id_of_most_stable)

            acquired_i = AL_gen_i.model.loc[init_id_of_most_stable].acquired


            if acquired_i:
                print("init_id_of_most_stable:", init_id_of_most_stable)
                out_dict["first_acq_id"] = init_id_of_most_stable
                gen_acquired = gen_i
                break

    elif color_progression_i is None:
        # gen_acquired = 69

        AL_gen_i = AL_i.al_gen_dict[
            list(AL_i.al_gen_dict.keys())[-1]
            ]
        model = AL_gen_i.model
        
        out_dict["first_acq_id"] = id_of_most_stable

        gen_acquired = model.loc[id_of_most_stable].gen_acquired

    print("Acquired:", gen_acquired)

    out_dict["gen_acq"] = gen_acquired

    return(out_dict)

# +
# %%capture

AL_dict = out_data_dict["AL_dict"]

runs_list = list(AL_dict.keys())

data_list = []
gen_acquired_list = []
for run_i in runs_list:
    print("run_i:", run_i)

    data_row_i = dict()

    AL_i = AL_dict[run_i]

    last_gen = list(AL_i.al_gen_dict.keys())[-1]

    tmp = AL_i.al_gen_dict[last_gen]
    duplicates = tmp.indices_that_are_duplicates

    not_dupl_list = []
    for index in tmp.model.index:
        if index not in duplicates:
            not_dupl_list.append(index)

    id_of_most_stable = tmp.model.loc[not_dupl_list].sort_values("y").iloc[0].name
    print("id_of_most_stable:", id_of_most_stable)

    # gen_acq = get_gen_acquired(
    out_dict = get_gen_acquired(
        AL_i,
        id_of_most_stable,
        )
    gen_acq = out_dict["gen_acq"]
    gen_acquired_list.append(gen_acq)

    data_list.append(out_dict)

    print("")
# -

print("Average generatios to acquire structure", "\n", np.mean(gen_acquired_list))

# +
# #####################################
# IrO2 ################################
# 3.74 | GP
# 6.33 | Random

# #####################################
# IrO3 ################################
# 4.33 | GP
# 4.82 | Random

# +
df = pd.DataFrame(data_list)

# df.first_acq_id.value_counts()
# df

for i in df.first_acq_id.unique().tolist():
    print(i)
    
    from IPython.display import display
    df_i = df[df.first_acq_id == i]
    # display(df_i)
    
    mean_gen_acq = df_i.gen_acq.mean()
    print("mean_gen_acq:", mean_gen_acq)
    print("")

# + active=""
#
#
#

# + Collapsed="false" jupyter={}
# import chart_studio.plotly as py
# import plotly.graph_objs as go

# from ccf_similarity.ccf import CCF

# from active_learning.al_analysis import ALPerformance

# from plotting.my_plotly import my_plotly_plot

# from layout import layout
