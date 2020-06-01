# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Surface Energy Pourbaix Analysis
#
# ***
#
# Remember to remove TEMP | Artifically add IrO3_battery row

# + [markdown] toc-hr-collapsed=true
# # Notebook Setup
# -

# ## Import Modules

# %%capture
# %load_ext autoreload
# %autoreload 2

# +
import sys
import os

sys.path.insert(
    0,
    os.path.join(
        os.environ["PROJ_irox"],
        "workflow"))

sys.path.insert(
    0,
    os.path.join(
        os.environ["PROJ_irox"],
        "data"))

sys.path.insert(
    0,
    ".")
from an_data_processing import load_df

###########################################################
###########################################################
import pandas as pd

import plotly as py
import chart_studio.plotly as py

import plotly.offline as py_off
import plotly.graph_objs as go
from plotly import tools


# My Modules

###########################################################
###########################################################

# Data Script Variable Import
from proj_data_irox import (
    bulk_e_per_atom_dict,
    )

from proj_data_irox import (
    h2_ref,
    h2o_ref,
    proj_dir_name,
    system_names_dict,
    data_dir,
    irox_bulk_color_map,
    irox_surface_e_color_map,
    bulk_pourb_trans_dict,
    )

from methods_surf_e import (
#     surf_e_4,
    make_color_subplot_list,
    process_row,
    add_convex_hull,
    )

# +
h2_ref

h2o_ref
# -

from surface_energy.surface_energy import surf_e_4

# ## Script Inputs

# +
O_mu_range = [-0., 2.4]
surf_e_range = [-0.1, 0.3]

opacity_rect = 0.3

save_plot = False

['001', '010', '100', '110', '111', '211']

smart_format_dict = [
    [{"facet": "001"}, {"dash": "solid"}],
    [{"facet": "010"}, {"dash": "32px,2px,32px,2px"}],
    [{"facet": "100"}, {"dash": "16px,2px,16px,2px"}],
    [{"facet": "110"}, {"dash": "8px,2px,8px,2px"}],
    [{"facet": "111"}, {"dash": "4px,2px,4px,2px"}],
    [{"facet": "211"}, {"dash": "2px,2px,2px,2px"}],
    ]

smart_format_dict = [
    [{"facet": "001"}, {"dash": "solid"}],
    [{"facet": "010"}, {"dash": "solid"}],
    [{"facet": "100"}, {"dash": "solid"}],
    [{"facet": "110"}, {"dash": "solid"}],
    [{"facet": "111"}, {"dash": "solid"}],
    [{"facet": "211"}, {"dash": "solid"}],
    ]

# + [markdown] toc-hr-collapsed=true
# # Load and Process Data
# -

# ## Load Data

# +
# %%capture

df_pourbaix, df_ads, df_surf = load_df(
    from_file=False, root_dir=data_dir,
    data_dir=data_dir, file_name="df_master.pickle",
    process_df=True)

df_m = df_surf
# -

# ## Process Data

# +
# Filter the jobs that were unsuccessful
df_m = df_m[[not i for i in pd.isna(df_m["elec_energy"].tolist())]]
df_m["name_i_3"] = df_m["name_i_2"] + "_" + df_m["layers"].apply(str)
df_m["surf_e_0"] = df_m.apply(
    surf_e_4,
    G_H2=h2_ref,
    G_H2O=h2o_ref,
    axis=1,
    )

df_m = df_m[df_m["job_type"] == "surface_coverage_energy"]
# -

# ## Removing Some Data to Simplify Plot

# +
rows_to_drop_indices = []

# indices_i = df_m[
#     (df_m["bulk_system"] == "IrO3") &
#     (df_m["coverage_type"] == "o_covered") &
#     (df_m["facet"] == "100") &
#     [True for i in range(len(df_m))]
#     ].index.tolist()
# assert len(indices_i) == 1, "More than 1 structure with these props"
# rows_to_drop_indices.append(indices_i[0])

# indices_i = df_m[
#     (df_m["bulk_system"] == "IrO3") &
#     (df_m["coverage_type"] == "o_covered") &
#     (df_m["facet"] == "211") &
#     [True for i in range(len(df_m))]
#     ].index.tolist()
# assert len(indices_i) == 1, "More than 1 structure with these props"
# rows_to_drop_indices.append(indices_i[0])

# #############################################################################
# #############################################################################

indices_i = df_m[
    (df_m["bulk_system"] == "IrO3") &
    (df_m["coverage_type"] == "bare") &
    (df_m["facet"] == "211") &
    [True for i in range(len(df_m))]
#     ]
    ].index.tolist()
assert len(indices_i) == 1, "More than 1 structure with these props"
rows_to_drop_indices.append(indices_i[0])

indices_i = df_m[
    (df_m["bulk_system"] == "IrO3") &
    (df_m["coverage_type"] == "bare") &
    (df_m["facet"] == "100") &
    [True for i in range(len(df_m))]
#     ]
    ].index.tolist()
assert len(indices_i) == 1, "More than 1 structure with these props"
rows_to_drop_indices.append(indices_i[0])


# #############################################################################
# #############################################################################

indices_i = df_m[
    (df_m["bulk_system"] == "IrO3_rutile-like") &
    (df_m["coverage_type"] == "bare") &
    (df_m["facet"] == "001") &
    [True for i in range(len(df_m))]
#     ]
    ].index.tolist()
assert len(indices_i) == 1, "More than 1 structure with these props"
rows_to_drop_indices.append(indices_i[0])

indices_i = df_m[
    (df_m["bulk_system"] == "IrO3_rutile-like") &
    (df_m["coverage_type"] == "o_covered") &
    (df_m["facet"] == "001") &
    [True for i in range(len(df_m))]
    ].index.tolist()
assert len(indices_i) == 1, "More than 1 structure with these props"
rows_to_drop_indices.append(indices_i[0])

indices_i = df_m[
    (df_m["bulk_system"] == "IrO3_rutile-like") &
    (df_m["coverage_type"] == "h_covered") &
    (df_m["facet"] == "001") &
    [True for i in range(len(df_m))]
    ].index.tolist()
assert len(indices_i) == 1, "More than 1 structure with these props"
rows_to_drop_indices.append(indices_i[0])

# #############################################################################
# #############################################################################


df_m = df_m.drop(labels=rows_to_drop_indices)
# -

# # Surface vs Oxygen Chemical Potential Plots

# +
traces_IrO2 = []
traces_IrO3 = []
traces_IrO3_rutile_like = []
traces_IrO3_battery = []

for i_cnt, row_i in df_m.iterrows():   
    trace = process_row(
        row_i,
        mesh_eval=False,
        xy_axis=("x", "y"),
        O_mu_range=O_mu_range,
        bulk_e_per_atom_dict=bulk_e_per_atom_dict,
        h2_ref=h2_ref,
        h2o_ref=h2o_ref,
        smart_format_dict=smart_format_dict,
        irox_surface_e_color_map=irox_surface_e_color_map,
        )

    if row_i["bulk_system"] == "IrO2":
        traces_IrO2.append(trace)
    elif row_i["bulk_system"] == "IrO3":
        traces_IrO3.append(trace)
    elif row_i["bulk_system"] == "IrO3_rutile-like":
        traces_IrO3_rutile_like.append(trace)
    elif row_i["bulk_system"] == "IrO3_battery":
        traces_IrO3_battery.append(trace)
# -

# # IrO2 Surface Energies vs V_RHE

# # Adding OER Equilibrium Line (1.23)

# +
trace_1_23 = go.Scatter(
    x= 2 * [1.23],
    y=surf_e_range,
    mode='lines',
    name="1.23",
    line = dict(
        color="black",
        width=1,
        dash="dot",
        ))

traces_IrO2.append(trace_1_23)
traces_IrO3.append(trace_1_23)
traces_IrO3_rutile_like.append(trace_1_23)
traces_IrO3_battery.append(trace_1_23)
# -

# # Adding Convex Hull

# +
num_mesh_points = 250

conv_hull_i = add_convex_hull(
    df_m[(df_m["bulk_system"] == "IrO2")],
    O_mu_range=O_mu_range,
    bulk_e_per_atom_dict=bulk_e_per_atom_dict,
    h2_ref=h2_ref,
    h2o_ref=h2o_ref,
    smart_format_dict=smart_format_dict,
    irox_surface_e_color_map=irox_surface_e_color_map,
    num_mesh_points=num_mesh_points,
    )
traces_IrO2.append(conv_hull_i)


conv_hull_i = add_convex_hull(df_m[
    (df_m["bulk_system"] == "IrO3")],
    O_mu_range=O_mu_range,
    bulk_e_per_atom_dict=bulk_e_per_atom_dict,
    h2_ref=h2_ref,
    h2o_ref=h2o_ref,
    smart_format_dict=smart_format_dict,
    irox_surface_e_color_map=irox_surface_e_color_map,
    num_mesh_points=num_mesh_points,
    )
traces_IrO3.append(conv_hull_i)


conv_hull_i = add_convex_hull(
    df_m[(df_m["bulk_system"] == "IrO3_rutile-like")],
    O_mu_range=O_mu_range,
    bulk_e_per_atom_dict=bulk_e_per_atom_dict,
    h2_ref=h2_ref,
    h2o_ref=h2o_ref,
    smart_format_dict=smart_format_dict,
    irox_surface_e_color_map=irox_surface_e_color_map,
    num_mesh_points=num_mesh_points,
    )
traces_IrO3_rutile_like.append(conv_hull_i)


conv_hull_i = add_convex_hull(
    df_m[(df_m["bulk_system"] == "IrO3_battery")],
    O_mu_range=O_mu_range,
    bulk_e_per_atom_dict=bulk_e_per_atom_dict,
    h2_ref=h2_ref,
    h2o_ref=h2o_ref,
    smart_format_dict=smart_format_dict,
    irox_surface_e_color_map=irox_surface_e_color_map,
    num_mesh_points=num_mesh_points,
    )
traces_IrO3_battery.append(conv_hull_i)

# + [markdown] toc-hr-collapsed=true
# # Plotting
# -

# ## Instantiate subplots

# +
fig = tools.make_subplots(
    rows=2,
    cols=2,
    vertical_spacing=0.1,
    subplot_titles=(
        system_names_dict["IrO2"],
        system_names_dict["IrO3"],
        system_names_dict["IrO3_rutile-like"],
        system_names_dict["IrO3_battery"],
        ),
    )

for trace_i in traces_IrO2:
    fig.append_trace(trace_i, 1, 1)    
for trace_i in traces_IrO3:
    fig.append_trace(trace_i, 1, 2)
for trace_i in traces_IrO3_rutile_like:
    fig.append_trace(trace_i, 2, 1)
for trace_i in traces_IrO3_battery:
    fig.append_trace(trace_i, 2, 2)

# fig = go.Figure(data=data, layout=layout)

# + [markdown] toc-hr-collapsed=false
# ## Plot Layout Settings
# -

# ### Common axis settings

# +
common_axis_dict = {
    "mirror": 'ticks',
    "zeroline": True,
    "showline": True,
    "linecolor": "black",
    "showgrid": False,
    
#     "autotick": False,
    "ticks": 'inside',
    "tick0": 0,
    "dtick": 0.5,
    "ticklen": 2,
    "tickwidth": 1,
    "tickcolor": 'black',
    
    "zerolinecolor": 'black',
    "zerolinewidth": 0.4,
    "linecolor": 'black',
    "linewidth": 1,
    }
        
common_xaxis_dict = {
    }

common_yaxis_dict = {
    "range": surf_e_range,
    }

# layout["title"] = "Surface Pourbaix Plots"
# layout["font"] = {"family": "Arial", "color": "black"}
# -

# ### Additional layout settings

# +
font_size_axis_title = 16 * (4/3)
font_size_subplot_title = 12. * (4/3)

# Changing the subplot title font size
for i in fig["layout"]["annotations"]:
    i["font"]["size"] = font_size_subplot_title 

fig["layout"].update({

    "font":dict(
        family='Arial',
#         size=18,
        color='black',
        ),
    
    "xaxis": common_axis_dict,
    "xaxis1": common_axis_dict,
    "xaxis2": common_axis_dict,
    "xaxis3": common_axis_dict,
    "xaxis4": common_axis_dict,
    
    "yaxis": dict(common_axis_dict, **common_yaxis_dict),
    "yaxis1": dict(common_axis_dict, **common_yaxis_dict),
    "yaxis2": dict(common_axis_dict, **common_yaxis_dict),
    "yaxis3": dict(common_axis_dict, **common_yaxis_dict),
    "yaxis4": dict(common_axis_dict, **common_yaxis_dict),

    
    "width": 18.7 * 37.795275591,
    "height": 18.7 * 37.795275591,

#     "width": 12.7 * 37.795275591,
#     "height": 12.7 * 37.795275591,

    "showlegend": False,
    })

fig["layout"]["annotations"] = fig["layout"]["annotations"] + \
    (
        dict(
            x=0.5,
            y=-0.11,
            showarrow=False,
            text='Voltage (V)',
            xref='paper',
            yref='paper',
            font=dict(
                color="black",
                size=font_size_axis_title,
                ),
            ),
        dict(
            x=-0.12,
            y=0.27,
            showarrow=False,
            text='Surface Free Energy (eV / A<sup>2</sup>)',
            textangle=-90,
            xref='paper',
            yref='paper',
            font=dict(
                color="black",
                size=font_size_axis_title,
                ),
            ),
        )
# -

# ### Adding in shade rectangle

# +
bulk_stability_shading = [] + \
    make_color_subplot_list(
        subplot_num=1,
        plot_range=O_mu_range,
        bulk_pourb_trans_dict=bulk_pourb_trans_dict,
        irox_bulk_color_map=irox_bulk_color_map, opacity=opacity_rect) + \
    make_color_subplot_list(
        subplot_num=2,
        plot_range=O_mu_range,
        bulk_pourb_trans_dict=bulk_pourb_trans_dict,
        irox_bulk_color_map=irox_bulk_color_map, opacity=opacity_rect) + \
    make_color_subplot_list(
        subplot_num=3,
        plot_range=O_mu_range,
        bulk_pourb_trans_dict=bulk_pourb_trans_dict,
        irox_bulk_color_map=irox_bulk_color_map, opacity=opacity_rect) + \
    make_color_subplot_list(
        subplot_num=4,
        plot_range=O_mu_range,
        bulk_pourb_trans_dict=bulk_pourb_trans_dict,
        irox_bulk_color_map=irox_bulk_color_map, opacity=opacity_rect)

fig["layout"]["shapes"] = bulk_stability_shading
# -

# ## Plot Out

# +
import datetime
now = datetime.datetime.now()
date_i = now.strftime("%Y-%m-%d")

if save_plot:
    save_dir = proj_dir_name
else:
    save_dir = "__temp__"

# py_off.plot(
#     fig,
#     filename=os.path.join("local_plots", "out_plot_00_" + date_i + ".html"))

# py.plotly.iplot(
#     fig,
#     filename=os.path.join(save_dir, "surface_pourbaix", "pl_surface_pourbaix_subplots"),
#     )
fig.show()

# +
fig["layout"]["width"] = 1.5 * 18.7 * 37.795275591
fig["layout"]["height"] = 1.5 * 18.7 * 37.795275591

# py.plotly.iplot(
#     fig,
#     filename=os.path.join(save_dir, "surface_pourbaix", "pl_surface_pourbaix_subplots__large"),
#     )

fig.show()
# -

# # Plotting Plots Separately

# +
common_yaxis_dict = {
    **common_yaxis_dict,
    **{
        "zeroline": False,
        "dtick": 0.1,
        }
    }

common_xaxis_dict = {
    **common_xaxis_dict,
    **{
        "zeroline": False,
        "showticklabels": False,
        }
    }

def plot_surf_pourb_indiv(traces_list, name="190320_TEMP"):
    data_m = traces_list

    for i in data_m:
        if "xaxis" in i.to_plotly_json().keys():
            i["xaxis"] = "x"
        if "yaxis" in i.to_plotly_json().keys():
            i["yaxis"] = "y"

            
    shapes_i = make_color_subplot_list(
        subplot_num=1,
        plot_range=O_mu_range,
        bulk_pourb_trans_dict=bulk_pourb_trans_dict,
        irox_bulk_color_map=irox_bulk_color_map, opacity=0.5)

    fig1 = dict(
        data=data_m,
        layout={
            "font":dict(
                family='Arial',
        #         size=18,
                color='black',
                ),
            "xaxis": dict(common_axis_dict, **common_xaxis_dict),
            "yaxis": dict(common_axis_dict, **common_yaxis_dict),

            "width": 7.5 * 37.795275591,
            "height": 5.6225 * 37.795275591,

            "margin": go.layout.Margin(
                l=30.0,
                r=10.0,
                b=50.0,
                t=50.0,
    #             pad=150.,
                ),
            "showlegend": False,

            
            "shapes": shapes_i,
            },
        )

    return((
        fig1,
        os.path.join(
            save_dir,
#             "__temp__",
            "surface_pourbaix", name)
        ))
# -

fig, filename = plot_surf_pourb_indiv(
    traces_IrO2,
    name="traces_IrO2")
py.plotly.iplot(fig, filename=filename)

fig, filename = plot_surf_pourb_indiv(
    traces_IrO3,
    name="traces_IrO3")
py.plotly.iplot(fig, filename=filename)

fig, filename = plot_surf_pourb_indiv(
    traces_IrO3_rutile_like,
    name="traces_IrO3_rutile_like")
# py.plotly.iplot(fig, filename=filename)

fig, filename = plot_surf_pourb_indiv(
    traces_IrO3_battery,
    name="traces_IrO3_battery")
# py.plotly.iplot(fig, filename=filename)

# + active=""
#
#
#
#
#
#
#
#

# +
# list(df_m)
# df_m["surf_e_0"]
# 100 rutile


df_m[
    (df_m["bulk_system"] == "IrO3") &
    (df_m["facet"] == "111") &
#     (df_m[""] == "") &
    [True for i in range(len(df_m))]
    ]

# df_m.loc[57]

# +
# shapes_i = make_color_subplot_list(
#     subplot_num=1,
#     plot_range=O_mu_range,
#     bulk_pourb_trans_dict=bulk_pourb_trans_dict,
#     irox_bulk_color_map=irox_bulk_color_map, opacity=0.8)


# fig["layout"]["shapes"] = bulk_stability_shading

# make_color_subplot_list(
#     subplot_num=1,
#     plot_range=O_mu_range,
#     bulk_pourb_trans_dict=bulk_pourb_trans_dict,
#     irox_bulk_color_map=irox_bulk_color_map, opacity=opacity_rect)
# -

surf_e_4
