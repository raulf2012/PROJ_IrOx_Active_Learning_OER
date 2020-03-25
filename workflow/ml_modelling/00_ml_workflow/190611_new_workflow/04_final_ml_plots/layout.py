"""
"""

# | - Import Modules
import sys
import os

sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"],
    "data"))

from proj_data_irox import (
    axis_label_font_size,
    axis_tick_labels_font_size,
    font_family,
    base_font_color,
    irox_bulk_color_map,
    main_systems,
    system_names_dict,
    )

import plotly.graph_objs as go
#__|

# #############################################################################

def get_layout(model=None):
    """
    """
    # | - get_layout

    # | - Process model data
    if model is not None:
        e_max = model["prediction_unstandardized"].max()
        e_min = model["prediction_unstandardized"].min()

        extra_space_y = abs(0.05 * (e_max - e_min))

        len_x_axis = model["x_axis_ind"].shape[0]
    else:
        e_max = -4.5
        e_min = -7

        extra_space_y = 0.05
        len_x_axis=300

    #__|

    # | - Main layout object
    layout = go.Layout(
        angularaxis=None,
        annotations=None,
        annotationdefaults=None,
        autosize=None,
        bargap=None,
        bargroupgap=None,
        barmode=None,
        barnorm=None,
        boxgap=None,
        boxgroupgap=None,
        boxmode=None,
        calendar=None,
        clickmode=None,
        coloraxis=None,
        colorscale=None,
        colorway=None,
        datarevision=None,
        direction=None,
        dragmode=None,
        editrevision=None,
        extendfunnelareacolors=None,
        extendpiecolors=None,
        extendsunburstcolors=None,
        font=go.layout.Font(
            color=base_font_color,
            family=font_family,
            size=None,
            ),
        funnelareacolorway=None,
        funnelgap=None,
        funnelgroupgap=None,
        funnelmode=None,
        geo=None,
        grid=None,
        # #########################################################################
        height=6 * 37.795275591,
        hiddenlabels=None,
        hiddenlabelssrc=None,
        hidesources=None,
        hoverdistance=None,
        hoverlabel=None,
        hovermode=None,
        images=None,
        imagedefaults=None,
        legend=None,
        mapbox=None,
        # #########################################################################
        margin=go.layout.Margin(
            autoexpand=None,
            b=50,
            l=60,
            pad=None,
            r=5,
            t=10,
            ),
        meta=None,
        metasrc=None,
        modebar=None,
        orientation=None,
        # #########################################################################
        paper_bgcolor="rgba(240,240,240,0.95)",
        # paper_bgcolor="rgba(240,240,240,0.0)",
        piecolorway=None,
        # #########################################################################
        plot_bgcolor="rgba(250,250,250,0.98)",
        # plot_bgcolor="rgba(250,250,250,0.0)",
        polar=None,
        radialaxis=None,
        scene=None,
        selectdirection=None,
        selectionrevision=None,
        separators=None,
        shapes=None,
        shapedefaults=None,
        # #########################################################################
        showlegend=False,
        sliders=None,
        sliderdefaults=None,
        spikedistance=None,
        sunburstcolorway=None,
        template=None,
        ternary=None,
        title=None,
        titlefont=None,
        transition=None,
        uirevision=None,
        updatemenus=None,
        updatemenudefaults=None,
        violingap=None,
        violingroupgap=None,
        violinmode=None,
        waterfallgap=None,
        waterfallgroupgap=None,
        waterfallmode=None,
        # #########################################################################
        width=18 * 37.795275591,
        xaxis=None,
        yaxis=None,
        )
    #__|

    # | - Axis Layout  options

    # | - shared axis dict
    shared_axis_dict = dict(
        anchor=None,
        automargin=None,
        autorange=None,
        calendar=None,
        categoryarray=None,
        categoryarraysrc=None,
        categoryorder=None,
        color=None,
        constrain=None,
        constraintoward=None,
        dividercolor=None,
        dividerwidth=None,
        domain=None,
        dtick=None,
        exponentformat=None,
        fixedrange=None,
        gridcolor=None,
        gridwidth=None,
        hoverformat=None,
        layer=None,
        # #########################################################################
        linecolor=base_font_color,
        linewidth=None,
        matches=None,
        # #########################################################################
        mirror=True,
        nticks=None,
        overlaying=None,
        position=None,
        range=None,
        rangemode=None,
        scaleanchor=None,
        scaleratio=None,
        separatethousands=None,
        showdividers=None,
        showexponent=None,
        # #########################################################################
        showgrid=False,
        # #########################################################################
        showline=True,
        showspikes=None,
        showticklabels=None,
        showtickprefix=None,
        showticksuffix=None,
        side=None,
        spikecolor=None,
        spikedash=None,
        spikemode=None,
        spikesnap=None,
        spikethickness=None,
        tick0=None,
        tickangle=None,
        tickcolor=None,
        # #########################################################################
        tickfont=dict(
            color=base_font_color,
            family=None,
            size=axis_tick_labels_font_size),
        tickformat=None,
        tickformatstops=None,
        tickformatstopdefaults=None,
        ticklen=None,
        tickmode=None,
        tickprefix=None,
        ticks="outside",
        tickson=None,
        ticksuffix=None,
        ticktext=None,
        ticktextsrc=None,
        tickvals=None,
        tickvalssrc=None,
        tickwidth=None,
        title=None,
        titlefont=None,
        type=None,
        uirevision=None,
        visible=None,
        # #########################################################################
        zeroline=True,
        # #########################################################################
        zerolinecolor=base_font_color,
        # #########################################################################
        zerolinewidth=1.,
        )
    #__|

    xaxis_layout = go.layout.XAxis(shared_axis_dict)
    xaxis_layout.update(go.layout.XAxis(
        # range=[0, 700],
        range=[0 - 5, len_x_axis + 5],
        rangeselector=None,
        rangeslider=None,
        ))

    yaxis_layout = go.layout.YAxis(shared_axis_dict)
    yaxis_layout.update(go.layout.YAxis(
        # title="ISJFIDJSF",
        # range=[-1.5, 2.5],
        range=[e_min - extra_space_y, e_max + extra_space_y],
        ))


    layout.xaxis = xaxis_layout
    layout.yaxis = yaxis_layout
    #__|


    # | - Plot Annotations
    annotations = [

        # | - Axis Titles
        {
            # 'font': {'size': axis_label_font_size},
            'font': {'size': 12},
            'showarrow': False,
            'text': 'IrO2 Candidate Space',
            'x': 0.5,
            'xanchor': 'center',
            'xref': 'paper',
            'y': 0,
            'yanchor': 'top',
            'yref': 'paper',
            'yshift': -30,
            },

        # {
        #     # 'font': {'size': axis_label_font_size},
        #     'font': {'size': 12},
        #     'showarrow': False,
        #     'text': 'Surface Free Energy (eV / A<sup>2</sup>)',
        #     'textangle': -90,
        #     'x': 0,
        #     'xanchor': 'right',
        #     'xref': 'paper',
        #     'xshift': -40,
        #     'y': 0.5,
        #     'yanchor': 'middle',
        #     'yref': 'paper'
        #     },
        #__|

        ]

    layout.annotations = annotations
    #__|


    return(layout)

    #__|
