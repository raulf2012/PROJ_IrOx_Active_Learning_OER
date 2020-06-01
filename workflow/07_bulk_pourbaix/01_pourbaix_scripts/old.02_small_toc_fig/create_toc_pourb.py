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

from plotting.my_plotly import my_plotly_plot

# +
# /mnt/f/Dropbox/01_norskov/00_git_repos/PROJ_IrOx_Active_Learning_OER/workflow/07_bulk_pourbaix/01_pourbaix_scripts

# /mnt/f/Dropbox/01_norskov/00_git_repos/PROJ_IrOx_Active_Learning_OER
# workflow/07_bulk_pourbaix/01_pourbaix_scripts
# -

# #########################################################
import pickle; import os
path_i = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/07_bulk_pourbaix/01_pourbaix_scripts",
    "out_data/pourb_fig_0.pickle")
with open(path_i, "rb") as fle:
    fig = pickle.load(fle)
# #########################################################

# +
fig.layout.annotations = None

fig.layout.xaxis.range = [0, 7]
# -

fig.layout.xaxis.tickfont.size = 6 * (4/3)
fig.layout.yaxis.tickfont.size = 6 * (4/3)

import plotly.graph_objs as go

# +
fig.layout.xaxis.title.font.size = 8 * (4/3)
fig.layout.yaxis.title.font.size = 8 * (4/3)
# fig.layout.xaxis.title.font.size = 6 * (4/3)

# # go.layout.xaxis.Title?

# +
# fig.layout.width = 120
fig.layout.width = 110

# fig.layout.height = 180
# fig.layout.height = 140
# fig.layout.height = 120
fig.layout.height = 110
# -

fig.show()

my_plotly_plot(
    figure=fig,
    plot_name="bulk_pourb_small_toc",
    write_html=True,
    write_png=False,
    png_scale=6.0,
    write_pdf=True,
    write_svg=False,
    try_orca_write=False,
    )
