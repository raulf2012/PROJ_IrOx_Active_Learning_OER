"""
"""

# | - Import Modules
import os
import sys

import plotly.graph_objs as go

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import (
    axis_label_font_size,
    axis_tick_labels_font_size,
    font_family,
    base_font_color,
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
        b=10.,
        l=10.,
        r=5,
        t=5,
        ),
    meta=None,
    metasrc=None,
    modebar=None,
    orientation=None,
    # #########################################################################
    # paper_bgcolor="rgba(250,250,250,0.9)",
    paper_bgcolor="white",
    plot_bgcolor="rgba(0,0,0,0)",
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

    width=1.5 * 7.964 * 37.795275591,

    # height=1.5 * 5.6002 * 37.795275591,
    # height=8.4 * 37.795275591,
    # height=7.7 * 37.795275591,
    height=7.67 * 37.795275591,
    
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
        color="black",
        family=None,
        size=axis_tick_labels_font_size,
        ),
    tickformat=None,
    tickformatstops=None,
    tickformatstopdefaults=None,
    ticklen=None,
    tickmode=None,  # 'auto', 'linear', 'array'
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

            size=axis_label_font_size,
            ),
        # text="TEMP",
        ),

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

xaxis_layout = go.layout.XAxis(shared_axis_dict)
xaxis_layout.update(go.layout.XAxis(
    #  range=[1.0, 2.0],
    title=dict(
        # text="ΔG<sub>O</sub> - ΔG<sub>OH</sub> (eV)",
        text="ΔG<sub>OH</sub> (eV)",
        ),
    ))

yaxis_layout = go.layout.YAxis(shared_axis_dict)
yaxis_layout.update(go.layout.YAxis(
    #  range=[2.0, 1.4],
    title=dict(
        # text="U<sub>RHE</sub> (V)",
        text="ΔG<sub>OH</sub>, ΔG<sub>O</sub>, ΔG<sub>OOH</sub> (eV)",
        ),
    ))

layout.xaxis = xaxis_layout
layout.yaxis = yaxis_layout
#__|


# | - Plot Annotations


annotations = [


    dict(
        x=1.3 - 0.08,
        y=1.93,
        xref="x",
        yref="y",
        # text="G<sub>OOH</sub>=G<sub>OH</sub>+3.2 <br> G<sub>O</sub> = 2 G<sub>OH</sub>",
        text="G<sub>OOH</sub>=G<sub>OH</sub>+3.2",
        showarrow=False,
        textangle=-45,
        font=dict(
            color="black",
            size=8,
            ),

        ),


    dict(
        x=1.2 - 0.08,
        y=1.93,
        xref="x",
        yref="y",
        # text="G<sub>OOH</sub>=G<sub>OH</sub>+3.2 <br> G<sub>O</sub> = 2 G<sub>OH</sub>",
        text="G<sub>OOH</sub>=G<sub>OH</sub>+3.0",
        showarrow=False,
        textangle=-45,
        font=dict(
            color="gray",
            size=8,
            ),

        ),


    ]

# layout.annotations = annotations
#__|
