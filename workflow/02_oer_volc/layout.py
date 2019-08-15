#!/usr/bin/env python

"""TEMP.

Author(s): Raul A. Flores
"""

#| - IMPORT MODULES
import os
import sys

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

import plotly.graph_objs as go

from proj_data_irox import (
    axis_label_font_size,
    axis_tick_labels_font_size,
    font_family,
    base_font_color,
    exp_irox_lim_pot,
    )

#__|


layout = go.Layout({

    # 'annotations': [{'showarrow': False,
    # 'text': 'IrO<sub>3</sub> (@1mA/cm<sup>2</sup>)',
    # 'x': 1,
    # 'xanchor': 'left',
    # 'xref': 'x',
    # 'y': 1.45,
    # 'yref': 'y',
    # 'yshift': 9},
    # {'showarrow': False,
    # 'text': 'IrO<sub>2</sub> (@1mA/cm<sup>2</sup>)',
    # 'x': 1,
    # 'xanchor': 'left',
    # 'xref': 'x',
    # 'y': 1.8,
    # 'yref': 'y',
    # 'yshift': 9}],

    'font': {'color': base_font_color, 'family': font_family},

    'width': 436.45228346974983,
    'height': 317.4916535470773,

    'legend': {
        'font': {'size': 18},
        'traceorder': 'normal',
        'x': 0.0, 'y': -0.1,
        'yanchor': 'top'},

    'margin': {'b': 50.0, 'l': 50.0, 'r': 50.0, 't': 50.0},
    'paper_bgcolor': 'rgba(250,250,250,0.9)',
    'plot_bgcolor': 'rgba(0,0,0,0)',
    'showlegend': False,

    'xaxis': {
        'linecolor': 'black',
        'mirror': True,
        'range': [1.0, 2.0],
        'showgrid': False,
        'showline': True,
        'showticklabels': True,
        'tickcolor': 'black',
        'tickfont': {'size': 10.666666666666666},
        'ticklen': 2,
        'ticks': 'inside',
        'tickwidth': 1,
        'title': {
            'font': {'size': 12.0},
            # 'text': 'ΔG<sub>O</sub> - ΔG<sub>OH</sub> (eV)',
            'text': '',
            },
        'zeroline': False},

    'yaxis': {
        'linecolor': 'black',
        'mirror': 'ticks',
        'range': [2.0, 1.4],
        'showgrid': False,
        'showline': True,
        'tickcolor': 'black',
        'tickfont': {'size': 10.666666666666666},
        'ticklen': 2,
        'ticks': 'inside',
        'tickwidth': 1,
        'title': {
            'font': {'size': 12.0},
            'text': 'Theoretical Limiting Potential (V)'},
        'zeroline': False,
        },

    })

#| - Annotations
annotations = []


annotations = [
    dict(
        x=1,
        y=exp_irox_lim_pot["iro3"]["lim_pot"],
        xref='x',
        yref='y',
        text='IrO<sub>3</sub> (@1mA/cm<sup>2</sup>)',
        showarrow=False,
        xanchor="left",
        yshift=9,
        ),

    dict(
        x=1,
        y=exp_irox_lim_pot["iro2"]["lim_pot"],
        xref='x',
        yref='y',
        text='IrO<sub>2</sub> (@1mA/cm<sup>2</sup>)',
        showarrow=False,
        xanchor="left",
        yshift=9,
        ),

    # dict(
    #     x=1,
    #     y=exp_irox_lim_pot["irox"]["lim_pot"],
    #     xref='x',
    #     yref='y',
    #     text='IrO<sub>x</sub>',
    #     showarrow=False,
    #     xanchor="left",
    #     yshift=7,
    #     ),

    ]


#| - Scaling equations on volcano
annotations.append(
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

    )

annotations.append(
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
    )
#__|

annotations.append(
    go.layout.Annotation(
        x=0.5,
        # y=-0.15,
        y=-0.2,
        showarrow=False,
        text="ΔG<sub>O</sub> - ΔG<sub>OH</sub> (eV)",

        font=dict(
            family=font_family,
            size=axis_label_font_size,
            color=base_font_color,
            ),

        xref="paper",
        yref="paper"
        )
    )



layout["annotations"] = annotations
#__|


#| - __old__
# Add system paths
# from an_data_processing import load_df
#
# ###########################################################
#
# # Python Modules
# import numpy as np
# import pandas as pd
#
# import chart_studio.plotly as py
# # import plotly.plotly as py
# import plotly.graph_objs as go
#
# # import colorlover as cl
# from IPython.display import HTML
#
# # My Modules
# from oxr_reaction.oxr_rxn import ORR_Free_E_Plot
# from oxr_reaction.oxr_plotting_classes.oxr_plot_scaling import Scaling_Relations_Plot
# from oxr_reaction.oxr_plotting_classes.oxr_plot_volcano import Volcano_Plot

# Project Data
#__|
