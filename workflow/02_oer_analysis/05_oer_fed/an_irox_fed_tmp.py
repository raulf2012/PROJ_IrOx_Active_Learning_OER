# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Free Energy Diagrams for OER on IrOx
#
# ***

# # Import Modules

# +
# %%capture
# %load_ext autoreload
# %autoreload 2

from IPython.display import display

import sys
import os

sys.path.insert(
    0,
    os.path.join(
        os.environ["PROJ_irox"],
        "workflow",
        ),
    )

sys.path.insert(
    0,
    os.path.join(
        os.environ["PROJ_irox"],
        "data",
        ),
    )


from an_data_processing import load_df

###########################################################

import plotly as py
import copy

# My Modules
from orr_reaction.orr_fed_plot import ORR_Free_E_Plot

# Project Data
from proj_data_irox import (
    smart_format_dict_FED,
    proj_dir_name,
    
    data_dir,
    )

smart_format_dict = smart_format_dict_FED
# -

# # Load Data

df_pourbaix, df_ads, df_surf = load_df(
    from_file=False,
    root_dir=data_dir,
    data_dir=data_dir,
    file_name="df_master.pickle",
    process_df=True,
    )
df_m = df_ads

# # Script Inputs

# +
save_plot = False

display_df = False
# -

# # Instantiate ORR_Free_E_Plot Class

# +
ORR_PLT = ORR_Free_E_Plot(
    free_energy_df=None,
    state_title="adsorbate",
    free_e_title="ads_e",  
    smart_format=smart_format_dict,
    bias=1.23,
    color_list=None,   
    show_H_e_pairs_annotations=True,
    show_legend=True,
    rxn_type="OER",
    )

ORR_PLT.ideal_ORR_series()

# + [markdown] toc-hr-collapsed=true
# # IrO3
# -

# ## 110

# ### IrO3 | 110 | O-covered

# +
df_i = df_m
df_i = df_i[df_i["bulk_system"] == "IrO3"]
df_i = df_i[df_i["facet"] == "110"]
df_i = df_i[df_i["job_type"] == "ORR_adsorption"]
df_i = df_i[df_i["coverage_type"] == "o_covered"]

if len(df_i) > 4:
    index_up = df_i[(df_i["adsorbate"] == "ooh") & (df_i["ooh_direction"] == "up")].index.tolist()[0]
    df_1 = df_i.drop([index_up])

    index_down = df_i[(df_i["adsorbate"] == "ooh") & (df_i["ooh_direction"] == "down")].index.tolist()[0]
    df_2 = df_i.drop([index_down])

    ORR_PLT.add_series(
        df_1,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0] + " | up",
        overpotential_type="OER",
        )

    ORR_PLT.add_series(
        df_2,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0] + " | down",
        overpotential_type="OER",
        )

elif len(df_i) == 4:

    ORR_PLT.add_series(
        df_i,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0],
        overpotential_type="OER",
        )

if display_df:
    display(df_i)
# -

# ### IrO3 | 110 | H-covered

# +
df_i = df_m
df_i = df_i[df_i["bulk_system"] == "IrO3"]
df_i = df_i[df_i["facet"] == "110"]
df_i = df_i[df_i["job_type"] == "ORR_adsorption"]
df_i = df_i[df_i["coverage_type"] == "h_covered"]


if len(df_i) > 4:
    index_up = df_i[(df_i["adsorbate"] == "ooh") & (df_i["ooh_direction"] == "up")].index.tolist()[0]
    df_1 = df_i.drop([index_up])

    index_down = df_i[(df_i["adsorbate"] == "ooh") & (df_i["ooh_direction"] == "down")].index.tolist()[0]
    df_2 = df_i.drop([index_down])

    ORR_PLT.add_series(
        df_1,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0] + " | up",
        overpotential_type="OER",
        )

    ORR_PLT.add_series(
        df_2,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0] + " | down",
        overpotential_type="OER",
        )

elif len(df_i) == 4:

    ORR_PLT.add_series(
        df_i,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0],
        overpotential_type="OER",
        )
if display_df:
    display(df_i)
# -

# ## 221

# ### IrO3 | 211 | O-covered

# +
df_i = df_m
df_i = df_i[df_i["bulk_system"] == "IrO3"]
df_i = df_i[df_i["facet"] == "211"]
df_i = df_i[df_i["job_type"] == "ORR_adsorption"]
df_i = df_i[df_i["coverage_type"] == "o_covered"]


if len(df_i) > 4:
    index_up = df_i[(df_i["adsorbate"] == "ooh") & (df_i["ooh_direction"] == "up")].index.tolist()[0]
    df_1 = df_i.drop([index_up])

    index_down = df_i[(df_i["adsorbate"] == "ooh") & (df_i["ooh_direction"] == "down")].index.tolist()[0]
    df_2 = df_i.drop([index_down])

    ORR_PLT.add_series(
        df_1,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0] + " | up",
        overpotential_type="OER",
        )

    ORR_PLT.add_series(
        df_2,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0] + " | down",
        overpotential_type="OER",
        )

elif len(df_i) == 4:

    ORR_PLT.add_series(
        df_i,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0],
        overpotential_type="OER",
        )

if display_df:
    display(df_i)
# -

# ### IrO3 | 211 | H-covered

# +
df_i = df_m
df_i = df_i[df_i["bulk_system"] == "IrO3"]
df_i = df_i[df_i["facet"] == "211"]
df_i = df_i[df_i["job_type"] == "ORR_adsorption"]
df_i = df_i[df_i["coverage_type"] == "h_covered"]

if len(df_i) > 4:
    index_up = df_i[(df_i["adsorbate"] == "ooh") & (df_i["ooh_direction"] == "up")].index.tolist()[0]
    df_1 = df_i.drop([index_up])

    index_down = df_i[(df_i["adsorbate"] == "ooh") & (df_i["ooh_direction"] == "down")].index.tolist()[0]
    df_2 = df_i.drop([index_down])

    ORR_PLT.add_series(
        df_1,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0] + " | up",
        overpotential_type="OER",
        )

    ORR_PLT.add_series(
        df_2,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0] + " | down",
        overpotential_type="OER",
        )

elif len(df_i) == 4:

    ORR_PLT.add_series(
        df_i,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0],
        overpotential_type="OER",
        )

if display_df:
    display(df_i)

# + [markdown] toc-hr-collapsed=true
# # IrO2
# -

# ## 100

# ### IrO2 | 100 | O-covered

# +
df_i = df_m
df_i = df_i[df_i["bulk_system"] == "IrO2"]
df_i = df_i[df_i["facet"] == "100"]
df_i = df_i[df_i["job_type"] == "ORR_adsorption"]
df_i = df_i[df_i["coverage_type"] == "o_covered"]

if len(df_i) > 4:
    index_up = df_i[(df_i["adsorbate"] == "ooh") & (df_i["ooh_direction"] == "up")].index.tolist()[0]
    df_1 = df_i.drop([index_up])

    index_down = df_i[(df_i["adsorbate"] == "ooh") & (df_i["ooh_direction"] == "down")].index.tolist()[0]
    df_2 = df_i.drop([index_down])

    ORR_PLT.add_series(
        df_1,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0] + " | up",
        overpotential_type="OER",
        )

    ORR_PLT.add_series(
        df_2,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0] + " | down",
        overpotential_type="OER",
        )

elif len(df_i) == 4:

    ORR_PLT.add_series(
        df_i,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0],
        overpotential_type="OER",
        )

if display_df:
    display(df_i)
# -

# ### IrO2 | 100 | H-covered

# +
df_i = df_m
df_i = df_i[df_i["bulk_system"] == "IrO2"]
df_i = df_i[df_i["facet"] == "100"]
df_i = df_i[df_i["job_type"] == "ORR_adsorption"]
df_i = df_i[df_i["coverage_type"] == "h_covered"]

if len(df_i) > 4:
    index_up = df_i[(df_i["adsorbate"] == "ooh") & (df_i["ooh_direction"] == "up")].index.tolist()[0]
    df_1 = df_i.drop([index_up])

    index_down = df_i[(df_i["adsorbate"] == "ooh") & (df_i["ooh_direction"] == "down")].index.tolist()[0]
    df_2 = df_i.drop([index_down])

    ORR_PLT.add_series(
        df_1,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0] + " | up",
        overpotential_type="OER",
        )

    ORR_PLT.add_series(
        df_2,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0] + " | down",
        overpotential_type="OER",
        )

elif len(df_i) == 4:

    ORR_PLT.add_series(
        df_i,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0],
        overpotential_type="OER",
        )

if display_df:
    display(df_i)
# -

# ## 110

# ### IrO2 | 110 | O-covered

# +
df_i = df_m
df_i = df_i[df_i["bulk_system"] == "IrO2"]
df_i = df_i[df_i["facet"] == "110"]
df_i = df_i[df_i["job_type"] == "ORR_adsorption"]
df_i = df_i[df_i["coverage_type"] == "o_covered"]

if len(df_i) > 4:
    index_up = df_i[(df_i["adsorbate"] == "ooh") & (df_i["ooh_direction"] == "up")].index.tolist()[0]
    df_1 = df_i.drop([index_up])

    index_down = df_i[(df_i["adsorbate"] == "ooh") & (df_i["ooh_direction"] == "down")].index.tolist()[0]
    df_2 = df_i.drop([index_down])

    ORR_PLT.add_series(
        df_1,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0] + " | up",
        overpotential_type="OER",
        )

    ORR_PLT.add_series(
        df_2,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0] + " | down",
        overpotential_type="OER",
        )

elif len(df_i) == 4:

    ORR_PLT.add_series(
        df_i,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0],
        overpotential_type="OER",
        )

if display_df:
    display(df_i)
# -

# ### IrO2 | 110 | H-covered

# +
df_i = df_m
df_i = df_i[df_i["bulk_system"] == "IrO2"]
df_i = df_i[df_i["facet"] == "110"]
df_i = df_i[df_i["job_type"] == "ORR_adsorption"]
df_i = df_i[df_i["coverage_type"] == "h_covered"]

if len(df_i) > 4:
    index_up = df_i[(df_i["adsorbate"] == "ooh") & (df_i["ooh_direction"] == "up")].index.tolist()[0]
    df_1 = df_i.drop([index_up])

    index_down = df_i[(df_i["adsorbate"] == "ooh") & (df_i["ooh_direction"] == "down")].index.tolist()[0]
    df_2 = df_i.drop([index_down])

    ORR_PLT.add_series(
        df_1,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0] + " | up",
        overpotential_type="OER",
        )

    ORR_PLT.add_series(
        df_2,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0] + " | down",
        overpotential_type="OER",
        )

elif len(df_i) == 4:

    ORR_PLT.add_series(
        df_i,
        plot_mode="all",
        opt_name=df_i["name_i"].tolist()[0],
        overpotential_type="OER",
        )

if display_df:
    display(df_i)
# -

# # Plotting

# +
if save_plot:
    save_dir = proj_dir_name
else:
    save_dir = "__temp__"

layout = ORR_PLT.plotly_fed_layout(
    plot_title="",
    plot_width=1.5 * 680,
    plot_height=2.5 * 510,
    )

py.plotly.iplot(
    {
        "data": ORR_PLT.plotly_data(),
        "layout": layout,
        },
        filename=os.path.join(save_dir, "pl_irox_fed_oer"),
    )

# +
# # VASP Gas-phase References
# h2_ref = -6.77014123
# h2o_ref = -14.21744725

# # h2_ref =  -6.759300
# # h2o_ref = -14.019771

# # Free Energy Corrections
# corrections_dict = {
#     "ooh": 0.34475,
#     "o": -0.0145,
#     "oh": 0.30225,
#     "bare": 0.,
#     }

# # corrections_dict = {
# #     "ooh": 0.,
# #     "o": 0.,
# #     "oh": 0.,
# #     "bare": 0.,
# #     }

# Jobs = DFT_Jobs_Analysis(
#     update_job_state=False,
#     job_type_class=None,
#     load_dataframe=True,
#     working_dir="./data",
#     dataframe_dir="./data",
#     )

# df_m = Jobs.filter_early_revisions(Jobs.data_frame)


# # Shorter path
# root_dir = "/global/cscratch1/sd/flores12/IrOx_Project"
# def calc_short_path(row, root_dir):
#     short_path = row["path"].replace(root_dir, "")[1:]
#     return(short_path)

# df_m["path_short"] = df_m.apply(calc_short_path, args=(root_dir,), axis=1)

# # System Specific Name (for legends)
# df_m["name_i"] = df_m["facet"] + " | " + df_m["coverage_type"] + " | " + df_m["bulk_system"]

# sys.path.insert(
#     0,
#     os.path.join(
#         os.environ["PROJ_irox"],
#         "scripts/01_Michal_OER_Plot_Script",
#         ),
#     )

# sys.path.insert(
#     0,
#     os.path.join(
#         os.environ["PROJ_irox"],
#         "workflow/an_analysis_dir",
#         ),
#     )

# sys.path.insert(
#     0,
#     os.path.join(
#         os.environ["PROJ_irox"],
#         "workflow/data",
#         ),
#     )

# import pandas as pd
# pd.set_option("display.max_columns", None)
# pd.set_option('display.max_rows', None)

# smart_format = [
    
#     [
#         {"bulk_system": "IrO3"},
#         {"color": "red"},
#         ],

    
#     [
#         {"bulk_system": "IrO2"},
#         {"color": "blue"},
#         ],
    
#     [
#         {"coverage_type": "o_covered"},
#         {"dash": "dot"},
#         ],

#     [
#         {"coverage_type": "h_covered"},
#         {"dash": None},
#         ],    
#     ]
