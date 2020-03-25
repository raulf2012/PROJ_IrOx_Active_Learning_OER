
# | - Import Modules
# import os
# import sys

# import pickle

# import numpy as np

# from ase.visualize import view
# from ase_modules.ase_methods import view_in_vesta


# #############################################################################
# import chart_studio.plotly as py
import plotly.graph_objs as go

# sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
# from proj_data_irox import (
#     static_irox_structures_path,
#     bulk_dft_data_path,
#     )
# from IPython.display import display

# #############################################################################
#__|

def plot_dij_matrix_heatmap(
    df_dij,
    d_thresh,
    e_thresh,

    ):
    # df_d_comp = df_dij_dft
    df_d_comp = df_dij

    max_d = df_d_comp.max().max()
    min_d = df_d_comp.min().min()
    data_range_z = (max_d - min_d)
    # print("data_range_z: ", data_range_z)

    d_thresh_stan = d_thresh * (1 / data_range_z)
    # print(d_thresh_stan)

    tmp = df_d_comp.index.tolist()

    colorscale_i = [
        [0.0, "black"],
        [0.000000001, "red"],
        [d_thresh_stan, "pink"],
        [d_thresh_stan + 0.000001, "rgb(220,220,220)"],
        [1.0, "white"],
        # [1.2, "blue"],
        ]
    # print(colorscale_i)

    # #############################################################################
    # #############################################################################
    trace_i = go.Heatmap(
        z=df_d_comp.values,
        x=df_d_comp.columns.tolist(),
        y=df_d_comp.index.tolist(),

        colorscale=colorscale_i,

        # [
        #     [0.0, "black"],
        #     [0.000000001, "red"],
        #     [d_thresh_stan, "pink"],
        #     [d_thresh_stan + 0.000001, "rgb(220,220,220)"],
        #     [1.0, "white"],
        #     [1.2, "blue"],
        #     ],

        # #####################################################################
        transpose=True,
        # #####################################################################
        xgap=0.,
        # #####################################################################
        ygap=0.,


        # | - Additional heatmap arguments
        autocolorscale=None,
        coloraxis=None,
        colorbar=None,
        connectgaps=None,
        customdatasrc=None,
        dx=None,
        dy=None,
        hoverinfo=None,
        hoverinfosrc=None,
        hoverlabel=None,
        hovertemplate=None,
        hovertemplatesrc=None,
        hovertext=None,
        hovertextsrc=None,
        ids=None,
        idssrc=None,
        meta=None,
        metasrc=None,
        name=None,
        opacity=None,
        reversescale=None,
        showscale=None,
        stream=None,
        text=None,
        textsrc=None,
        uid=None,
        uirevision=None,
        visible=None,
        x0=None,
        xaxis=None,
        xcalendar=None,
        xsrc=None,
        xtype=None,
        y0=None,
        yaxis=None,
        ycalendar=None,
        ysrc=None,
        ytype=None,
        zauto=None,
        zhoverformat=None,
        zmax=None,
        zmid=None,
        zmin=None,
        zsmooth=None,
        #__|

        )

    data = [trace_i]


    return(data)

    # from plotting.my_plotly import my_plotly_plot
    # # layout = go.Layout(width=1000, height=1000)
    # layout = go.Layout(width=1100, height=1100)
    # fig = go.Figure(data=data, layout=layout)
    #
    # fig = my_plotly_plot(
    #     figure=fig,
    #     plot_name='irox_dij_heatmap',
    #     write_pdf_svg=True,
    #     write_html=True,
    #     write_png=False,
    #     write_pdf=True,
    #     write_svg=False,
    #     )

    return(fig)
