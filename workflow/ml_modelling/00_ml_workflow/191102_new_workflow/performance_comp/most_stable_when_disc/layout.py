"""
"""

# | - Import Modules
import plotly.graph_objs as go
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
    # #########################################################################
    font=go.layout.Font(
        color="black",
        family="Arial",
        # size=4 * (4/3),
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
        b=10,
        l=10,
        pad=None,
        r=10,
        t=10,
        ),
    meta=None,
    metasrc=None,
    modebar=None,
    orientation=None,
    # #########################################################################
    # paper_bgcolor=None,
    paper_bgcolor="rgba(0,0,0,0.5)",
    piecolorway=None,
    # #########################################################################
    # plot_bgcolor=None,
    plot_bgcolor="rgba(0,0,0,0.5)",
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
    height=37.795275591 * 7.05,
    width=37.795275591 * 5,

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
    linecolor="black",
    linewidth=1.,

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

    # #########################################################################
    tickfont=dict(
        color="black",
        # family=None,
        # size=10 * (4/3),
        size=8 * (4/3),
        ),

    tickformat=None,
    tickformatstops=None,
    tickformatstopdefaults=None,
    ticklen=None,
    tickmode=None,
    tickprefix=None,
    # #########################################################################
    ticks="outside",

    tickson=None,
    ticksuffix=None,
    ticktext=None,
    ticktextsrc=None,
    tickvals=None,
    tickvalssrc=None,
    tickwidth=None,
    title=None,
    titlefont=dict(
        color="black",
        # family=None,
        # size=10 * (4/3),
        size=10 * (4/3),
        ),

    type=None,
    uirevision=None,
    visible=None,
    zeroline=False,
    zerolinecolor=None,
    zerolinewidth=None,
    )
#__|

xaxis_layout = go.layout.XAxis(shared_axis_dict)
xaxis_layout.update(go.layout.XAxis(
    title="DFT Calculations",
    rangeselector=None,
    rangeslider=None,
    ))

yaxis_layout = go.layout.YAxis(shared_axis_dict)
yaxis_layout.update(go.layout.YAxis(
    title="# of Top 10 Polymorphs Discovered",
    tick0=0,
    dtick=2,
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
        # 'font': {'size': axis_label_font_size},
        'font': {'size': 12},
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

# layout.annotations = annotations
#__|
