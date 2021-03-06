# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_IrOx_Active_Learning_OER]
#     language: python
#     name: conda-env-PROJ_IrOx_Active_Learning_OER-py
# ---

# # Import Modules

# +
import os
print(os.getcwd())
import sys

import chart_studio.plotly as py
import plotly.graph_objs as go

# #########################################################
# Local Imports
from layout import layout
# -

# # Read Data

# #########################################################
import pickle; import os
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/opt_mae_err_gp_model",
    "out_data/AB3_data.pickle")
    # "out_data/data_0.pickle")
with open(path_i, "rb") as fle:
    df = pickle.load(fle)
# #########################################################

df.head()

# +
x_array = df.pca_comp
# y_array = df.mae_cv
y_array = df.mae_ave


trace = go.Scatter(
    x=x_array,
    y=y_array,
    mode="markers+lines",
    opacity=1.,
    marker=dict(
        symbol="circle",
        color='grey',
        opacity=1.,
        size=12,
        line=dict(
            color="black",
            width=1,
            )
        ),
    line=dict(
        color="black",
        width=2,
        dash="solid",
        ),
    )

data = [trace]

fig = go.Figure(data=data, layout=layout)
# -

# # Saving Plot

# +
from plotting.my_plotly import my_plotly_plot

my_plotly_plot(
    figure=fig,
    plot_name="mae_vs_pca",
    write_html=True,
    write_png=False,
    png_scale=6.0,
    write_pdf=True,
    write_svg=False,
    try_orca_write=True,
    )
# -

print(20 * "# # ")
print("All done!")
assert False

fig.show()
