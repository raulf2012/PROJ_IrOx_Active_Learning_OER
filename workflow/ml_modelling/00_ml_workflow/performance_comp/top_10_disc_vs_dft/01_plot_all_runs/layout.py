"""
"""

#| - Import Modules
import os
import sys

import plotly.graph_objs as go

sys.path.insert(0, os.path.join(
    os.environ["PROJ_irox"], "data"))
from proj_data_irox import axis_label_font_size, axis_tick_labels_font_size
#__|

# #########################################################

#| - Main layout object
layout = go.Layout(
    angularaxis=None,
    # annotations=None,
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
	# b=None,
	# b=40,
	b=50,
	# l=10,
	l=60,
	pad=None,
	r=10,
	t=10,
	),
    meta=None,
    metasrc=None,
    modebar=None,
    orientation=None,
    # #########################################################################
    #  paper_bgcolor="rgba(255,255,255,0.5)",
    #  plot_bgcolor="rgba(255,255,255,0.5)",
    paper_bgcolor="white",
    plot_bgcolor="white",

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

    # #####################################################
    # height=5.291667 * 37.795275591,
    # height=10. * 37.795275591,
    #  height=9.5 * 37.795275591,
    #  height=9.7 * 37.795275591,
    height=9.8 * 37.795275591,

    # width=17.7 * 37.795275591,
    # width=16.5 * 37.795275591,
    width=16.3 * 37.795275591,

    )
#__|


#| - Axis Layout  options

#| - shared axis dict
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
    showticklabels=False,
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
        color=None,
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

xaxis_layout = go.layout.XAxis(shared_axis_dict)
xaxis_layout.update(go.layout.XAxis(
    title=dict(
        # text="DFT Calculations",
        ),
    # range=[0, 450],
    ))

yaxis_layout = go.layout.YAxis(shared_axis_dict)
yaxis_layout.update(go.layout.YAxis(
    title=dict(
        # text="N<sub>discovered</sub>",
        ),
    range=[-0.5, 11],
    ))

layout.xaxis = xaxis_layout
layout.yaxis = yaxis_layout
#__|


#| - Plot Annotations
annotations = [
    go.layout.Annotation({
	"font": {"size": 16},
	"showarrow": False,
	"text": "DFT Calculations",
	"x": 0.5,
	"xanchor": "center",
	"xref": "paper",
	"y": 0,
	"yanchor": "top",
	"yref": "paper",
	"yshift": -30
	}),

    go.layout.Annotation({
	"font": {"size": 16},
	"showarrow": False,
	"text": "N<sub>discovered</sub>",
	"textangle": -90,
	"x": 0,
	"xanchor": "right",
	"xref": "paper",
	"xshift": -40,
	"y": 0.5,
	"yanchor": "middle",
	"yref": "paper"
	}),

    ]

layout.annotations = annotations
#__|
