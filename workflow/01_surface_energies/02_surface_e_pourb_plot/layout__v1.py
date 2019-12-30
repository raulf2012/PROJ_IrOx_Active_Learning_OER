"""Plotly layout attributes for Surface Pourbaix Energy plot.

Author(s): Raul A. Flores
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
    voltage_name,
    )

import plotly.graph_objs as go
#__|

# #############################################################################
# #############################################################################
axis_num_list = [1, 3, 5, 7]
# axis_num_list = [1, 2, 3, 4]

#| - Main layout object
layout = go.Layout(

    #| - height/width
    width=8.0 * 37.795275591,

    # height=15.0 * 37.795275591,
    # height=14.8 * 37.795275591,

    # ########
    # height=14.85 * 37.795275591,  # short
    # height=14.92 * 37.795275591,  # long
    # height=14.89 * 37.795275591, # long

    # height=14.87 * 37.795275591, # long
    # height=14.86 * 37.795275591, # long
    height=14.85 * 37.795275591,  # long
    #__|

    font=go.layout.Font(
        color=base_font_color,
        family=font_family,
        # size,
        ),

    margin=go.layout.Margin(
        autoexpand=None,
        # b=40,
        # b=35,
        # b=34,
        # b=34.5,
        # b=34.5,  # <-------

        # b=38,
        # b=40,

        # TEMP
        # COMBAK

        # b=100,
        # b=60,
        # b=55,
        # b=50,
        # b=45,
        # b=40,
        # b=35,
        b=39,

        # l=60,
        # l=55,
        l=53,

        # r=8,
        r=6,

        # t=10,
        # t=15,
        # t=18,
        # t=15,
        # t=3,
        # t=9,
        # t=12,
        t=13,

        pad=None,
        ),

    # paper_bgcolor="rgba(240,240,240,0.95)",
    #  paper_bgcolor="rgba(250,250,250,0.4)",
    # paper_bgcolor="rgba(200,200,200,0.4)",
    paper_bgcolor="white",

    # plot_bgcolor="rgba(250,250,250,0.98)",

    showlegend=False,
    piecolorway=None,

    #| - __old__
    # angularaxis=None,
    # annotations=None,
    # annotationdefaults=None,
    # autosize=None,
    # bargap=None,
    # bargroupgap=None,
    # barmode=None,
    # barnorm=None,
    # boxgap=None,
    # boxgroupgap=None,
    # boxmode=None,
    # calendar=None,
    # clickmode=None,
    # coloraxis=None,
    # colorscale=None,
    # colorway=None,
    # datarevision=None,
    # direction=None,
    # dragmode=None,
    # editrevision=None,
    # extendfunnelareacolors=None,
    # extendpiecolors=None,
    # extendsunburstcolors=None,
    # funnelareacolorway=None,
    # funnelgap=None,
    # funnelgroupgap=None,
    # funnelmode=None,
    # geo=None,
    # grid=None,
    # hiddenlabels=None,
    # hiddenlabelssrc=None,
    # hidesources=None,
    # hoverdistance=None,
    # hoverlabel=None,
    # hovermode=None,
    # images=None,
    # imagedefaults=None,
    # legend=None,
    # mapbox=None,

    # meta=None,
    # metasrc=None,
    # modebar=None,
    # orientation=None,
    # polar=None,
    # radialaxis=None,
    # scene=None,
    # selectdirection=None,
    # selectionrevision=None,
    # separators=None,
    # shapes=None,
    # shapedefaults=None,
    # sliders=None,
    # sliderdefaults=None,
    # spikedistance=None,
    # sunburstcolorway=None,
    # template=None,
    # ternary=None,
    # title=None,
    # titlefont=None,
    # transition=None,
    # uirevision=None,
    # updatemenus=None,
    # updatemenudefaults=None,
    # violingap=None,
    # violingroupgap=None,
    # violinmode=None,
    # waterfallgap=None,
    # waterfallgroupgap=None,
    # waterfallmode=None,
    # xaxis=None,
    # yaxis=None,
    #__|

    )
#__|

#| - Axis layout options

#| - shared axis dict
shared_axis_dict = dict(

    #| - __old__
    # anchor=None,
    # automargin=None,
    # autorange=None,
    # calendar=None,
    # categoryarray=None,
    # categoryarraysrc=None,
    # categoryorder=None,
    # color=None,
    # constrain=None,
    # constraintoward=None,
    # dividercolor=None,
    # dividerwidth=None,
    # domain=None,
    # dtick=None,
    # exponentformat=None,
    # fixedrange=None,
    # gridcolor=None,
    # gridwidth=None,
    # hoverformat=None,
    # layer=None,
    #     linewidth=None,
    #     matches=None,
    # nticks=None,
    # overlaying=None,
    # position=None,
    # range=None,
    # rangemode=None,
    # scaleanchor=None,
    # scaleratio=None,
    # separatethousands=None,
    # showdividers=None,
    # showexponent=None,
    # showspikes=None,
    # showtickprefix=None,
    # showticksuffix=None,
    # side=None,
    # spikecolor=None,
    # spikedash=None,
    # spikemode=None,
    # spikesnap=None,
    # spikethickness=None,
    # tick0=None,
    # tickangle=None,
    # tickson=None,
    # ticksuffix=None,
    # ticktext=None,
    # ticktextsrc=None,
    # tickvals=None,
    # tickvalssrc=None,
    # tickwidth=None,
    # title=None,
    # titlefont=None,
    # type=None,
    # uirevision=None,
    # visible=None,
    # tickformat=None,
    # tickformatstops=None,
    # tickformatstopdefaults=None,
    # ticklen=None,
    # tickmode=None,
    # tickprefix=None,
    #__|

    linecolor="black",
    mirror=True,
    showgrid=False,
    showline=True,
    tickfont=dict(
        # color='crimson',
        size=axis_tick_labels_font_size),

    ticks="outside",
    tickcolor="black",

    zeroline=True,
    zerolinecolor="black",
    zerolinewidth=1.,
    )
#__|

xaxis_layout = go.layout.XAxis(shared_axis_dict)
xaxis_layout.update(go.layout.XAxis(
    # range=[0, 2.4],
    range=[0.5, 2.4],
    showticklabels=False,
    ))

yaxis_layout = go.layout.YAxis(shared_axis_dict)
yaxis_layout.update(go.layout.YAxis(
    # range=[-0.12, +0.22],
    # range=[-0.18, +0.22],
    # range=[-0.22, +0.22],
    # range=[-0.25, +0.22],
    range=[-0.28, +0.22],
    dtick=0.1,
    ))


#| - Setting x/y-axis to layout

tmp = dict(zip(
    ["xaxis" + str(i) for i in axis_num_list],
    4 * [xaxis_layout],
    ))
layout.update(tmp)

tmp = dict(zip(
    ["yaxis" + str(i) for i in axis_num_list],
    4 * [yaxis_layout],
    ))
layout.update(tmp)


#| - __old__
# layout.xaxis = xaxis_layout
# layout.yaxis = yaxis_layout
#
# layout.xaxis3 = xaxis_layout
# layout.yaxis3 = yaxis_layout
#
# layout.xaxis5 = xaxis_layout
# layout.yaxis5 = yaxis_layout
#
# layout.xaxis7 = xaxis_layout
# layout.yaxis7 = yaxis_layout


# layout.xaxis = xaxis_layout
# layout.yaxis = yaxis_layout
#
# layout.xaxis2 = xaxis_layout
# layout.yaxis2 = yaxis_layout
#
# layout.xaxis3 = xaxis_layout
# layout.yaxis3 = yaxis_layout
#
# layout.xaxis4 = xaxis_layout
# layout.yaxis4 = yaxis_layout
#__|

#__|

#__|

#| - Plot Annotations

#| - shared_annot_props
shared_annot_props = go.layout.Annotation(

    # ###########################################
    font=go.layout.annotation.Font(
        color="black",
        family=None,
        size=None,
        ),
    # ###########################################
    showarrow=False,

    # align=None,
    # arrowcolor=None,
    # arrowhead=None,
    # arrowside=None,
    # arrowsize=None,
    # arrowwidth=None,
    # ax=None,
    # axref=None,
    # ay=None,
    # ayref=None,
    # bgcolor=None,
    # bordercolor=None,
    # borderpad=None,
    # borderwidth=None,
    # captureevents=None,
    # clicktoshow=None,


    # height=None,
    # hoverlabel=None,
    # hovertext=None,
    # name=None,
    # opacity=None,


    # standoff=None,
    # startarrowhead=None,
    # startarrowsize=None,
    # startstandoff=None,
    # templateitemname=None,
    # text=None,
    # textangle=None,
    # valign=None,
    # visible=None,
    # width=None,
    # x=None,
    # xanchor=None,
    # xclick=None,
    # xref=None,
    # xshift=None,
    # y=None,
    # yanchor=None,
    # yclick=None,
    # yref=None,
    # yshift=None,
    )
#__|

annotations = [

    #| - Axis Titles

    go.layout.Annotation(
        # COMBAK TODO Is this actually a "free energy"?
        text="Surface Free Energy (eV / A<sup>2</sup>)",
        font=go.layout.annotation.Font(
            size=axis_label_font_size,),
        x=0.0,
        y=0.5,
        textangle=-90,
        xref="paper",
        yref="paper",
        # yanchor="top",
        xanchor="left",

        # COMBAK
        # xshift=-55,
        # xshift=-50,
        # xshift=-52,
        # xshift=-53,
        # xshift=-58,
        xshift=-56,

        ).update(**shared_annot_props.to_plotly_json()),



    go.layout.Annotation(
        text=voltage_name,
        font=go.layout.annotation.Font(
            size=axis_label_font_size,),
        x=0.5,
        y=0.0,
        xref="paper",
        yref="paper",
        yanchor="top",

        # yshift=-20,
        # yshift=-18,
        # yshift=-30,

        # COMBAK
        # yshift=-25,
        yshift=-21,

        name="x_axis_title",

        ).update(**shared_annot_props.to_plotly_json()),



    #__|

    ]

#| - subplot bulk label in righ-top corner
subplot_label_dict = go.layout.Annotation(
    x=layout.xaxis.range[1],
    y=layout.yaxis.range[1],
    xanchor="right",
    yanchor="top",
    yshift=0.50,
    xshift=0.40,
    font=dict(
        family=font_family,
        size=axis_tick_labels_font_size,
        color="black",
        ),
    showarrow=False,
    )

for axis_i, (i_cnt, system_i) in zip(axis_num_list, enumerate(main_systems)):

    annot_i_tmp = go.layout.Annotation(
        text=system_names_dict[system_i],
        xref="x" + str(axis_i),
        yref="y" + str(axis_i),

        # xref="x" + str(i_cnt + 1),
        # yref="y" + str(i_cnt + 1),
        bgcolor=irox_bulk_color_map[system_i],
        )

    annot_i = go.layout.Annotation(
        **subplot_label_dict.to_plotly_json(),
        ).update(annot_i_tmp.to_plotly_json())

    annotations.append(annot_i)

#__|

layout.annotations = annotations
#__|
