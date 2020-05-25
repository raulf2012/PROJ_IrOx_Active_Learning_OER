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

# # Figuring out how many computational resources I've used
#
# Computational resources are counted by parsing the OUTCAR file and searching for the
# line with 'Total CPU time used'.

# # Import Modules

# %%capture
# %load_ext autoreload
# %autoreload 2

# +
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

import pandas as pd
pd.set_option("display.max_columns", None)
pd.set_option('display.max_rows', None)

from ase.visualize import view

from Pourbaix_simple_generic import plot_Pourbaix
from an_data_processing import load_df
from ase_modules.ase_methods import create_species_element_dict
# -

# # Read/Process DataFrame

# +
    # /mnt/c/Users/raulf/Dropbox/01_acad_folder/01_grad_school/01_norskov/04_comp_clusters/02_DATA/04_IrOx_surfaces_OER

data_dir = os.path.join(
os.environ["dropbox"],
"01_acad_folder/01_grad_school/01_norskov/04_comp_clusters/02_DATA/04_IrOx_surfaces_OER",
)
# -

df_master = load_df(
    from_file=False,
    
    root_dir="../data",
    data_dir="../../data",
    file_name="df_master.pickle",
    process_df=False,
    filter_early_revisions=False,
    )

df_master

len(df_master)

2 + 2


def parse_cpu_time(row_i):
    outcar_list = row_i.outcar
    search_lines = [i for i in outcar_list if "Total CPU time used" in i]
    if len(search_lines) == 1:
        time_i = float(search_lines[0].split()[-1])
        return(time_i)    
    else:
        return(0.)


df_master = df_master[df_master["outcar"].notnull()]

# +
df_master["cpu_time"] = df_master.apply(
    parse_cpu_time,
    axis=1,
    )

df_master["cpu_hours"] = (df_master["cpu_time"] * 240.) / 3600.

# +
df_m = df_master

print(df_m["cpu_hours"].mean())
print(df_m["cpu_hours"].sum())
# -

110700 * 2. * 1.4
