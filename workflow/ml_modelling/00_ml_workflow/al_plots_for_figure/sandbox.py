# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python [conda env:research-new]
#     language: python
#     name: conda-env-research-new-py
# ---

# +
color_list = [
    "0c007fff",
    "0000acff",
    "0001d7ff",
    "0033e9ff",
    "0053ffff",
    "0073ffff",
    "008df3ff",
    "00b5f6ff",
    "00dcf5ff",
    "00fff3ff",
    ]

color_list

# +
# h = input('Enter hex: ').lstrip('#')

h = '0c007fff'


for color_i in color_list:
    h = color_i
    print(tuple(int(h[i:i+2], 16) for i in (0, 2, 4)))

# +
import plotly.graph_objs as go

# # go.Figure?

# go.Frame?


# # go.Layout?
# fig.update_layout(
#     annotations=[
#         go.layout.Annotation(
#             x=50,
#             y=2,
#             xref="x",
#             yref="y",
#             text="dict Text",
#             showarrow=True,
#             arrowhead=7,
#             ax=0,
#             ay=-40
#             )
#         ],

#     )

# +
# go.Scatter.textposition?
