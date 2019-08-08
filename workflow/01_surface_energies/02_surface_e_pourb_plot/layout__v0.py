"""

TODO:
    zeroline color to black



8cmx15cm


mirror=True,
ticks='outside',
showline=True,

8.0 * 37.795275591
15.0 * 37.795275591

1cm/37.795275591px
"""

#| - Import Modules
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
# #############################################################################

#| - Layout *******************************************************************

#| - Main layout object
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
        # size,
        ),
    funnelareacolorway=None,
    funnelgap=None,
    funnelgroupgap=None,
    funnelmode=None,
    geo=None,
    grid=None,
    # #########################################################################
    # height=15.0 * 37.795275591,
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
    paper_bgcolor="rgba(240,240,240,0.95)",
    piecolorway=None,
    # #########################################################################
    plot_bgcolor="rgba(250,250,250,0.98)",
    polar=None,
    radialaxis=None,
    scene=None,
    selectdirection=None,
    selectionrevision=None,
    separators=None,
    shapes=None,
    shapedefaults=None,
    showlegend=None,
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
    # width=8.0 * 37.795275591,
    xaxis=None,
    yaxis=None,
    )
#__|

#| - Axis layout options

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
    # #########################################################################
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
        # color='crimson',
        size=axis_tick_labels_font_size),
    tickformat=None,
    tickformatstops=None,
    tickformatstopdefaults=None,
    ticklen=None,
    tickmode=None,
    tickprefix=None,
    ticks=None,
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
    zerolinecolor="black",
    # #########################################################################
    zerolinewidth=1.,
    )
#__|

xaxis_layout = go.layout.XAxis(shared_axis_dict)
xaxis_layout.update(go.layout.XAxis(
    range=[0, 2.4],
    rangeselector=None,
    rangeslider=None,
    ))

yaxis_layout = go.layout.YAxis(shared_axis_dict)
yaxis_layout.update(go.layout.YAxis(
    range=[-0.35, +0.3],
    ))

layout.xaxis = xaxis_layout
layout.yaxis = yaxis_layout

layout.xaxis2 = xaxis_layout
layout.yaxis2 = yaxis_layout

layout.xaxis3 = xaxis_layout
layout.yaxis3 = yaxis_layout

layout.xaxis4 = xaxis_layout
layout.yaxis4 = yaxis_layout
#__|

#| - Plot Annotations
annotations = [

    #| - Axis Titles
    {
        'font': {'size': axis_label_font_size},
        'showarrow': False,
        'text': 'Voltage (V vs RHE)',
        'x': 0.5,
        'xanchor': 'center',
        'xref': 'paper',
        'y': 0,
        'yanchor': 'top',
        'yref': 'paper',
        'yshift': -30,
        },

    {
        'font': {'size': axis_label_font_size},
        'showarrow': False,
        'text': 'Surface Free Energy (eV / A<sup>2</sup>)',
        'textangle': -90,
        'x': 0,
        'xanchor': 'right',
        'xref': 'paper',
        'xshift': -40,
        'y': 0.5,
        'yanchor': 'middle',
        'yref': 'paper'
        },
    #__|

    ]

subplot_label_dict = dict(
    text="TEMP",

    x=layout.xaxis.range[1],
    y=layout.yaxis.range[1],

    xref="x1",
    yref="y1",

    xanchor="right",
    yanchor="top",

    yshift=0.50,
    xshift=0.40,

    font=dict(
        family=font_family,
        size=axis_tick_labels_font_size,
        color="black",
        ),

    bgcolor=irox_bulk_color_map["IrO2"],
    opacity=1.0,
    showarrow=False,
    # align="center",
    )


for i_cnt, system_i in enumerate(main_systems):
    annotations.append(
        go.layout.Annotation(
            **subplot_label_dict,
            ).update(
                text=system_names_dict[system_i],
                xref="x" + str(i_cnt + 1),
                yref="y" + str(i_cnt + 1),
                bgcolor=irox_bulk_color_map[system_i],
                )
        )

layout.annotations = annotations
#__|

#| - Plot Shapes

subplot_label_rect_dict = dict(
    type="rect",
    # x0=1.5,
    # y0=0.2,
    x0=2,
    y0=0.2,
    # x1=2,
    # y1=0.1,
    x1=layout.xaxis.range[1],
    y1=layout.yaxis.range[1],
    xref="x2",
    yref="y2",
    line=dict(
        color="black",
        width=1,
        ),
    fillcolor="pink",
    )

shapes=[

    go.layout.Shape(
        **subplot_label_rect_dict,
        ).update(
            xref="x1",
            yref="y1",
            fillcolor=irox_bulk_color_map["IrO2"],
            ),

    go.layout.Shape(
        **subplot_label_rect_dict,
        ).update(
            xref="x2",
            yref="y2",
            fillcolor=irox_bulk_color_map["IrO3"],
            ),

    go.layout.Shape(
        **subplot_label_rect_dict,
        ).update(
            xref="x3",
            yref="y3",
            fillcolor=irox_bulk_color_map["IrO3_battery"],
            ),

    go.layout.Shape(
        **subplot_label_rect_dict,
        ).update(
            xref="x4",
            yref="y4",
            fillcolor=irox_bulk_color_map["IrO3_rutile-like"],
            ),

    ]

# layout.shapes = shapes
#__|

#__| **************************************************************************

# #############################################################################



#| - __old__
# annotations.append(
#     go.layout.Annotation(
#         **subplot_label_dict,
#         ).update(
#             text="TEMP",
#             xref="x1",
#             yref="y1",
#             bgcolor=irox_bulk_color_map["IrO2"],
#             )
#     )
#
# annotations.append(
#     go.layout.Annotation(
#         **subplot_label_dict,
#         ).update(
#             text="TEMP1",
#             xref="x2",
#             yref="y2",
#             bgcolor=irox_bulk_color_map["IrO3"],
#             )
#     )
#
# annotations.append(
#     go.layout.Annotation(
#         **subplot_label_dict,
#         ).update(
#             text="TEMP1",
#             xref="x3",
#             yref="y3",
#             bgcolor=irox_bulk_color_map["IrO3_battery"],
#             )
#     )
#
# annotations.append(
#     go.layout.Annotation(
#         **subplot_label_dict,
#         ).update(
#             text="TEMP1",
#             xref="x4",
#             yref="y4",
#             bgcolor=irox_bulk_color_map["IrO3_rutile-like"],
#             )
#     )
#__|

#| - __old__
# def get_layout(
#
#
#     ):
#     """
#     """
#     #| - get_layout
#
#
#     #__|
#__|
