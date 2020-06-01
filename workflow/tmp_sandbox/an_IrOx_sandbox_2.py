# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernel_info:
#     name: python3
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# # Import Modules

# +
# %%capture
# %load_ext autoreload
# %autoreload 2

import sys
import os

sys.path.insert(
    0,
    os.path.join(
        os.environ["PROJ_irox"],
        "scripts",
        ),
    )

sys.path.insert(
    0,
    os.path.join(
        os.environ["PROJ_irox"],
        "workflow/an_analysis_dir",
        ),
    )

sys.path.insert(
    0,
    os.path.join(
        os.environ["PROJ_irox"],
        "workflow/data",
        ),
    )

from Pourbaix_simple_generic import plot_Pourbaix
from an_data_processing import load_df
from ase_modules.ase_methods import create_species_element_dict

################################################################
################################################################
################################################################

import pandas as pd
pd.set_option("display.max_columns", None)
pd.set_option('display.max_rows', None)

from ase import io
from ase.visualize import view

import numpy as np
# -

# # Read/Process DataFrame

df_pourbaix, df_ads, df_surf = load_df(
    from_file=False,
    root_dir=os.path.join(os.environ["PROJ_irox"], "workflow/data"),
    data_dir=os.path.join(os.environ["PROJ_irox"], "workflow/data"),

    file_name="df_master.pickle",
    process_df=True,
    )

df_surf

from vasp.vasp_methods import parse_incar

parse_incar(df_ads.iloc[0].incar)

df_ads.iloc[0].incar_parsed

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
#             
#             
