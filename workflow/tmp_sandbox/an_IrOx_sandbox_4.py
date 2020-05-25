# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernel_info:
#     name: python3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Testing new ORR_Series code to handle more than 1 energy per state

# +
# %%capture
# %load_ext autoreload
# %autoreload 2

import sys
import os

sys.path.insert(
    0,
    os.path.join(
        os.environ["PROJ_irox"],
        "scripts/01_Michal_OER_Plot_Script",
        ),
    )

sys.path.insert(
    0,
    os.path.join(
        os.environ["PROJ_irox"],
        "workflow/an_analysis_dir",
        ),
    )

sys.path.insert(
    0,
    os.path.join(
        os.environ["PROJ_irox"],
        "workflow/data",
        ),
    )


import plotly as py
import pandas as pd
pd.set_option("display.max_columns", None)
pd.set_option('display.max_rows', None)

import copy

# My Modules
from orr_reaction.orr_fed_plot import ORR_Free_E_Plot
from orr_reaction.orr_series import ORR_Free_E_Series

from an_data_processing import load_df

df_pourbaix, df_ads, df_surf = load_df(
    from_file=False,
    root_dir=os.path.join(os.environ["PROJ_irox"], "workflow/data"),
    data_dir=os.path.join(os.environ["PROJ_irox"], "workflow/data"),
    file_name="df_master.pickle",
    process_df=True,
    )
df_m = df_ads

smart_format = [
    
    [
        {"bulk_system": "IrO3"},
        {"color": "red"},
        ],

    
    [
        {"bulk_system": "IrO2"},
        {"color": "blue"},
        ],
    
    [
        {"coverage_type": "o_covered"},
        {"dash": "dot"},
        ],

    [
        {"coverage_type": "h_covered"},
        {"dash": None},
        ],    
    ]

ORR_PLT = ORR_Free_E_Plot(
    free_energy_df=None,
    state_title="adsorbate",
    free_e_title="ads_e",  
    smart_format=smart_format,
    bias=0.,
#     bias=1.23 + 0.330,   
    color_list=None,   
    show_H_e_pairs_annotations=True,
    show_legend=True,
    rxn_type="OER",
    )

ORR_PLT.ideal_ORR_series()

# +
df_i = df_m
df_i = df_i[df_i["bulk_system"] == "IrO3"]
df_i = df_i[df_i["facet"] == "110"]
df_i = df_i[df_i["job_type"] == "ORR_adsorption"]
df_i = df_i[df_i["coverage_type"] == "o_covered"]

df_2 = copy.deepcopy(df_i)
df_2 = df_2.drop([234])
df_2 = df_2.drop([229])
# -

ORR_PLT.add_series(
    df_i,
    plot_mode="all",
    opt_name=df_i["name_i"].tolist()[0],
    overpotential_type="OER",
    )

ORR_PLT.series_list[-1].energy_lst_new

# +
# ORR_Free_E_Series(
#     free_energy_df=df_2,
#     # system_properties=None,
#     state_title="adsorbate",
#     free_e_title="ads_e",

#     bias=0.,
#     rxn_x_coord_array=None,
#     opt_name=None,
#     properties=None,
#     color_list=None,
#     color=None,
#     i_cnt=0,
#     hover_text_col=None,
#     plot_mode="all",
#     smart_format=None,
#     # overpotential_type="ORR",
#     rxn_type="ORR",
#     )

# +
# dir(ORR_PLT.series_list[-1])
# orr_i = ORR_PLT.series_list[-1]
# orr_i.energy_states_dict

# df = df_i
# rxn_mech_states = ["bulk", "ooh", "o", "oh", "bulk"]
# # rxn_mech_states --> self.rxn_mech_states
# for state in rxn_mech_states:
#     df_state = df.loc[df["adsorbate"] == state]

# + active=""
#  
#   
#    
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#     
#    
#        
#         
#           
#           
# -

# # Testing with a correctly formated DF

# +
df_2 = df_2.drop([234])
df_2 = df_2.drop([229])

ORR_PLT.add_series(
    df_2,
    plot_mode="all",
    opt_name=df_i["name_i"].tolist()[0],
    overpotential_type="OER",
    )
