# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernel_info:
#     name: python3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Investigating OOH structures and scaling

# +
# %%capture
# %load_ext autoreload
# %autoreload 2

# from IPython.display import display

import sys
import os

an_dir = os.path.join(
    os.environ["PROJ_irox"],
    "workflow")
sys.path.insert(0, an_dir)

# +
import sys
import os

an_dir = os.path.join(
    os.environ["PROJ_irox"],
    "data")
sys.path.insert(0, an_dir)

# -

from proj_data_irox import (
    corrections_dict_tmp,
    corrections_dict,
    )

corrections_dict_tmp

corrections_dict

# + active=""
#
#
#
#
#
#
#
#
#
#
#

# +
# ads_spec = "oh"

# corr_i = (0. + \

#     + ( 0. + \
#         + tmp_dict[ads_spec]["zpe"] +  \
#         + tmp_dict[ads_spec]["cv"] + \
#         - tmp_dict[ads_spec]["ts"]  
#         ) + \

#     - (0. + \
#         + h2o_corr + \
#         + ((1. - 1. * 2.) / 2.) * h2_corr
#         )
#     )


# print(tmp)

# +
# dataframe_dir = os.path.join(
#     os.environ["PROJ_irox"],
#     "data",
#     )

# from an_data_processing import load_df
# df_pourbaix, df_ads, df_surf = load_df(
#     from_file=False,
#     root_dir=dataframe_dir,
#     data_dir=dataframe_dir,  
#     file_name="df_master.pickle",
#     process_df=True,
#     )
# df_m = df_ads

# +
# import pandas as pd
# pd.set_option("display.max_columns", None)
# pd.set_option('display.max_rows', None)
# from ase import io
# from orr_reaction.orr_fed_plot import ORR_Free_E_Plot
# from orr_reaction.orr_fed_plot import Scaling_Relations_Plot
# import plotly.plotly as py
# import copy
