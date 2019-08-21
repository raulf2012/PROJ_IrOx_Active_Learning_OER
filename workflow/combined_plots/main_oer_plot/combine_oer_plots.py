# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.4'
#       jupytext_version: 1.2.1
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Import Modules

# +
import os

# #############################################################################
from svgutils.compose import Figure, Panel, SVG, Text, Line

# #############################################################################
from misc_modules.image_processing import convert_pdf_to_svg

# #############################################################################
from IPython.display import SVG as iSVG
from IPython.display import display
# -

# ## File Paths

# +
figure_file_paths = []

# #############################################################################
surf_pourb_plot_path = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/01_surface_energies/02_surface_e_pourb_plot/out_plot",
    "surf_e_pourbaix_irox__regular.pdf")
figure_file_paths.append(surf_pourb_plot_path)

# #############################################################################
oer_volc_path = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/02_oer_volc/out_plot",
    "pl_irox_volcano_plotly_default_ooh.pdf",
#     "TEMP_878.pdf",
    )
figure_file_paths.append(oer_volc_path)

# #############################################################################
subplot_c = os.path.join(
    "/home/raulf2012/Dropbox/01_norskov/01_projects",
    "04_irox_oer_orr/00_figures/00_combined_1d_volcano_surf_pourb/subplot_c_structure_vis",
    "00_subplot_c__v2.pdf",
#     "00_subplot_c__v1.svg",
#     "00_subplot_c__v0.svg",
    )
figure_file_paths.append(subplot_c)
# -

# # Converting pdf to svg 

# +
directory = "out_data/converted_svgs"
if not os.path.exists(directory):
    os.makedirs(directory)

out_dir = directory
pdf_figure_file_paths = []
for figure_path_i in figure_file_paths:
    pdf_figure_path_i = convert_pdf_to_svg(figure_path_i, out_dir)
    pdf_figure_file_paths.append(pdf_figure_path_i)
# -

# # Creating Figure Object

fig = Figure(
    "20cm", "15cm",
    Panel(
        SVG(pdf_figure_file_paths[0]),
        ).move(0, 0),
    Panel(
        SVG(pdf_figure_file_paths[1]),
        ).move(300, -39.609448819368 + 0.007999999999981355),
    Panel(
        SVG(pdf_figure_file_paths[2]),
#         ).move(330, 290),
#         ).move(330, 290),
        ).move(310 - 2.948999882 + 0.09991455078, 270 + 2.344163895),
    )

fig.save("out_data/fig_final_compose.svg")
# display(iSVG("fig_final_compose.svg"))

os.system("google-chrome out_data/fig_final_compose.svg")
