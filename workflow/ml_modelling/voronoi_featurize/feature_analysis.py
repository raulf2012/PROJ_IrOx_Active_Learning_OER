# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_IrOx_Active_Learning_OER]
#     language: python
#     name: conda-env-PROJ_IrOx_Active_Learning_OER-py
# ---

# # Import Modules

# +
import os
import sys

import pickle

import plotly.express as px


from plotting.my_plotly import my_plotly_plot

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import (
    voronoi_features_data_path,
    voronoi_features_all_data_path,
    )
# -

# # Import Data

# +
with open(voronoi_features_all_data_path, "rb") as fle:
    df_voronoi = pickle.load(fle)

with open(voronoi_features_data_path, "rb") as fle:
    df_voronoi_pca = pickle.load(fle)

# +
df_voronoi.describe()

df_voronoi.head()
# -

df_voronoi_pca.describe()

# +
# df_i = df_voronoi_pca

# fig = px.scatter_matrix(df_i)
# fig.show()

# fig.layout["height"] = 1200
# fig = my_plotly_plot(figure=fig, plot_name="scatter_matrix_tmp")

# +
# df_i = df_voronoi

# fig = px.scatter_matrix(df_i)
# fig.show()

# fig.layout["height"] = 1200
# fig = my_plotly_plot(figure=fig, plot_name="scatter_matrix_voronoi_full")
