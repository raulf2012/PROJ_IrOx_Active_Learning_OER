# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import sys
import os

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
from proj_data_irox import proj_dir_name

import plotly.plotly as py
import plotly.graph_objs as go

# Create random data with numpy
import numpy as np

# +
# smart_format_dict = [
# #     [{"facet": "A"}, {"dash": "10px,5px,2px,5px"}],
#     [{"facet": "A"}, {"dash": "4px,2px,4px,2px"}],
#     [{"facet": "B"}, {"dash": "12px,2px,12px,2px"}],
#     [{"facet": "C"}, {"dash": "20px,2px,20px,2px"}],
#     [{"facet": "D"}, {"dash": "28px,2px,28px,2px"}],
#     [{"facet": "D"}, {"dash": "36px,2px,36px,2px"}],
#     ]

# smart_format_dict = [
#     [{"facet": "001"}, {"dash": "4px,2px,4px,2px"}],
#     [{"facet": "010"}, {"dash": "12px,2px,12px,2px"}],
#     [{"facet": "100"}, {"dash": "20px,2px,20px,2px"}],
#     [{"facet": "110"}, {"dash": "28px,2px,28px,2px"}],
#     [{"facet": "111"}, {"dash": "36px,2px,36px,2px"}],
#     [{"facet": "211"}, {"dash": "44px,2px,44px,2px"}],
#     ]

#     [{"facet": "010"},  {"dash": "10px,5px,2px,5px"}],
#     [{"facet": "001"},  {"dash": "16px,3px,16px,3px"}],
#     [{"facet": "100"},  {"dash": "solid"}],
#     [{"facet": "110"},  {"dash": "dot"}],
#     [{"facet": "111"},  {"dash": "dash"}],
#     [{"facet": "211"},  {"dash": "dashdot"}],

# smart_format_dict = [
#     [{"facet": "001"}, {"dash": "4px,2px,4px,2px"}],
#     [{"facet": "010"}, {"dash": "12px,2px,12px,2px"}],
#     [{"facet": "100"}, {"dash": "20px,2px,20px,2px"}],
#     [{"facet": "110"}, {"dash": "28px,2px,28px,2px"}],
#     [{"facet": "111"}, {"dash": "36px,2px,36px,2px"}],
#     [{"facet": "211"}, {"dash": "44px,2px,44px,2px"}],
#     ]

# smart_format_dict = [
#     [{"facet": "001"}, {"dash": "2px,2px,2px,2px"}],
#     [{"facet": "010"}, {"dash": "4px,2px,4px,2px"}],
#     [{"facet": "100"}, {"dash": "8px,2px,8px,2px"}],
    
#     [{"facet": "110"}, {"dash": "16px,2px,16px,2px"}],
#     [{"facet": "111"}, {"dash": "32px,2px,32px,2px"}],
#     [{"facet": "211"}, {"dash": "solid"}],
#     ]

smart_format_dict = [
    [{"facet": "001"}, {"dash": "32px,0px,32px,0px"}],
    [{"facet": "010"}, {"dash": "32px,2px,32px,2px"}],
    [{"facet": "100"}, {"dash": "16px,2px,16px,2px"}],
    [{"facet": "110"}, {"dash": "8px,2px,8px,2px"}],
    [{"facet": "111"}, {"dash": "4px,2px,4px,2px"}],
    [{"facet": "211"}, {"dash": "2px,2px,2px,2px"}],
    ]
# -

['001', '010', '100', '110', '111', '211']


# +
format_dict = {}
for i in smart_format_dict:
    format_dict[i[0]["facet"]] = i[1]["dash"]

data = []
for ind, i in enumerate(format_dict.items()):
    trace_i = go.Scatter(
        x = [0, 0.2], y = [ind + 1., ind + 1.],
        mode = 'lines',
        name = i[0],
        line = dict(
            dash=i[1],
            width=4.,
            color="black",
            ),
        )
    data.append(trace_i)
# -

trace_i = go.Scatter(
    x = [0, 0.2], y = [0., 0.],
    mode = 'lines',
    name = i[0],
    line = dict(
#         dash="solid",
        width=1.,
        color="black",
        ),
    )
data.append(trace_i)

# +
layout = {
    "width": 100 * 2.5,
    "height": 100 * 3.,
    "showlegend": False,
    }

fig = go.Figure(data=data, layout=layout)
# -

py.iplot(fig,
    filename=os.path.join(
        proj_dir_name,
        "surface_pourbaix",
        "line_type_legend",
        )
    )

# + active=""
#
#

# +
# data = []
# i=1

# trace_i = go.Scatter(
#     x = [0, 1], y = [i, i],
#     mode = 'lines',
#     name = '100',
#     line = dict(
#         dash="solid",
#         width=2.,
#         color="black",
#         ),
#     )
# data.append(trace_i)
# i += 1

# trace_i = go.Scatter(
#     x = [0, 1],
#     y = [i, i],
#     mode = 'lines',
#     name = '110',
#     line = dict(
#         dash="dot",
#         width=2.,
#         color="black",
#         ),
#     )
# data.append(trace_i)
# i += 1

# trace_i = go.Scatter(
#     x = [0, 1],
#     y = [i, i],
#     mode = 'lines',
#     name = '111',
#     line = dict(
#         dash="dash",
#         width=2.,
#         color="black",
#         ),
#     )
# data.append(trace_i)
# i += 1

# trace_i = go.Scatter(
#     x = [0, 1],
#     y = [i, i],
#     mode = 'lines',
#     name = '211',
#     line = dict(
#         dash="dashdot",
#         width=2.,
#         color="black",
#         ),
#     )
# data.append(trace_i)
# i += 1

# trace_i = go.Scatter(
#     x = [0, 1],
#     y = [i, i],
#     mode = 'lines',
#     name = '001',
#     line = dict(
#         dash="dashdot",
#         width=2.,
#         color="black",
#         ),
#     )
# data.append(trace_i)
# i += 1
