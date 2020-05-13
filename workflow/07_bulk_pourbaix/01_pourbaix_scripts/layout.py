"""
"""

# | - Import Modules
import os
import sys

import copy

import plotly.graph_objs as go

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import (
    axis_label_font_size,
    axis_tick_labels_font_size,
    font_family,
    base_font_color,
    irox_bulk_color_map,
    )

#__|

# #############################################################################


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
        color="black",
        family="Arial",
        size=None,
        ),
    funnelareacolorway=None,
    funnelgap=None,
    funnelgroupgap=None,
    funnelmode=None,
    geo=None,
    grid=None,
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
    margin=go.layout.Margin(
        autoexpand=None,
        b=0,
        l=0,
        pad=None,
        r=0,
        t=0,
        ),
    meta=None,
    metasrc=None,
    modebar=None,
    orientation=None,
    # #########################################################################
    paper_bgcolor="rgba(255,255,255,0.9)",
    plot_bgcolor="rgba(255,255,255,0.9)",

    # paper_bgcolor="rgba(230,230,230,0.5)",
    # plot_bgcolor="rgba(255,255,255,0.5)",

    piecolorway=None,
    polar=None,
    radialaxis=None,

    # #########################################################################
    scene=go.layout.Scene(
        annotations=None,
        annotationdefaults=None,
        aspectmode=None,
        aspectratio=None,
        bgcolor=None,
        camera=None,
        domain=None,
        dragmode=None,
        hovermode=None,
        uirevision=None,
        xaxis=None,
        yaxis=None,
        zaxis=None,
        ),

    selectdirection=None,
    selectionrevision=None,
    separators=None,
    shapes=None,
    shapedefaults=None,
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
    xaxis=None,
    yaxis=None,

    # #########################################################################
    # height=5.291667 * 37.795275591,
    # height=10 * 37.795275591,
    height=7 * 37.795275591,

    width=8.5 * 37.795275591,

    # height=None,
    # width=None,
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
    linecolor="black",
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
    tickcolor="black",
    tickfont=dict(
        color=base_font_color,
        family=None,
        size=axis_tick_labels_font_size,
        ),
    tickformat=None,
    tickformatstops=None,
    tickformatstopdefaults=None,
    ticklen=None,
    tickmode=None,
    tickprefix=None,
    ticks="outside",  # outside", "inside", ""
    tickson=None,
    ticksuffix=None,
    ticktext=None,
    ticktextsrc=None,
    tickvals=None,
    tickvalssrc=None,
    tickwidth=None,
    title=dict(
        font=dict(
            color="black",
            family=None,
            size=None,
            ),
        text="TEMP",
        ),

    # go.layout.Title(
    #     arg=None,
    #     font=None,
    #     pad=None,
    #     text=None,
    #     x=None,
    #     xanchor=None,
    #     xref=None,
    #     y=None,
    #     yanchor=None,
    #     yref=None,
    #     # font=None, text=None
    #     ),

    titlefont=None,
    type=None,
    uirevision=None,
    visible=None,

    # #########################################################################
    zeroline=False,
    zerolinecolor=None,
    zerolinewidth=None,
    )
#__|

#    tick0=None,

xaxis_layout = go.layout.XAxis(shared_axis_dict)
xaxis_layout.update(go.layout.XAxis(
    dtick=2.,
    tick0=1,
    title=dict(
        text="pH",
        ),
    ))

yaxis_layout = go.layout.YAxis(shared_axis_dict)
yaxis_layout.update(go.layout.YAxis(
    title=dict(
        text="U<sub>SHE</sub> (V)",
        ),
    ))

layout.xaxis = xaxis_layout
layout.yaxis = yaxis_layout
#__|


# | - Plot Annotations

shared_annot_props = go.layout.Annotation(
    font={"size": axis_label_font_size},
    showarrow=False,
    xanchor="left",
    xref="paper",
    yanchor="top",
    yref="paper",
    yshift=0.,
    )


dy = -0.13
dx0 = 0.11
annotations = [
    # #########################################################################
    # Ir aqueous ion ##########################################################
    copy.deepcopy(shared_annot_props).update(
        text="IrO<sub>4</sub><sup>-</sup> (aq)",
        font_color="white",
        x=0.55,
        y=0.75,
        overwrite=True),


    # #########################################################################
    # IrO3 ####################################################################

    # a-IrO3
    copy.deepcopy(shared_annot_props).update(
        text="α-AlF<sub>3</sub> IrO<sub>3</sub> (s)",
        font_color=irox_bulk_color_map["IrO3_a-AlF3"],
        x=0.02,
        y=0.98 + dy,
        overwrite=True),

    # rutile-IrO3
    copy.deepcopy(shared_annot_props).update(
        text="R-IrO<sub>3</sub> (s)",
        font_color=irox_bulk_color_map["IrO3_rutile-like"],
        x=0.02 + dx0,
        y=0.91 + dy - 0.005,
        overwrite=True),

    # beta-IrO3
    copy.deepcopy(shared_annot_props).update(
        text="β-IrO<sub>3</sub> (s)",
        font_color=irox_bulk_color_map["IrO3_battery"],
        x=0.02 + 2 * dx0,
        y=0.84 + dy - 0.01,
        overwrite=True),

    # #########################################################################
    # rutile-IrO2 #############################################################
    copy.deepcopy(shared_annot_props).update(
        # text="rutile IrO<sub>2</sub> (s)",
        text="R-IrO<sub>2</sub> (s)",
        x=0.5,
        y=0.23,
        overwrite=True),

    # #########################################################################
    # Ir metal ################################################################
    copy.deepcopy(shared_annot_props).update(
        text="Ir (s)",
        x=0.15,
        y=0.15,
        overwrite=True),


    ]


layout.annotations = annotations
#__|
