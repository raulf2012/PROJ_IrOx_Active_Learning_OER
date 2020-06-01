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

# + [markdown] Collapsed="false"
# # TEMP
# ---

# + [markdown] Collapsed="false"
# # Import Modules

# + Collapsed="false" attributes={"classes": [], "id": "", "n": "1"}
import os

import pandas as pd

from ase.db import connect

# + [markdown] Collapsed="false"
# # Read Data

# + Collapsed="false" jupyter={}
# Structural Analysis db file
FinalStructuresdb_file = os.path.join(
    os.environ["PROJ_irox"],
    "workflow/ml_modelling/energy_vs_volume/kirsten_E_vs_V_analysis/scripts/out_data",
    "FinalStructures_1.db")
db = connect(FinalStructuresdb_file)


# + [markdown] Collapsed="false"
# # Construct DataFrame
#
#


# +
data_list = []
for row in db.select():
    row_dict = dict(
        **row.key_value_pairs)
    data_list.append(row_dict)

df = pd.DataFrame(data_list)
df = df[~df["stoich"].isna()]

df = df.drop(columns=["stoich", "id_old"])
# -

df

# Pickling data ###########################################
import os; import pickle
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)
with open(os.path.join(directory, "df_coord_env.pickle"), "wb") as fle:
    pickle.dump(df, fle)
# #########################################################

# + active=""
#
#
#

# + Collapsed="false" attributes={"classes": [], "id": "", "n": "1"} jupyter={}
# structure_id_map = {
#     '64cg6j9any': 'i (rutile)',  # rutile
#     'cg8p7fxq65': 'anatase', # anatase
#     'm2bs8w82x5': 'brookite', # Brookite
#     'n36axdbw65': 'ii', # 2nd stable columbite like?
#     '85z4msnl6o': 'iii (pyrite)', # Pyrite                    
#     #'myc4ng73xh': 'v', # Fm3m
#     'zizr7rvpxs': 'vi', # Porous
#     'b49kx4c19q': 'v (columbite)', # Columbite
#     'nscdbpmdct': 'iv',  # P63 (layered)                    
#     #'m2bs8w82x5': 'vi',
#     # IrO3
#     'mp6lno9jzr': 'i', # 482_2d
#     '9i6ixublcr': 'iii', # porous
#     'v2blxebixh': 'ii', # sg=2
#     'nrml6dms9l': 'iv',   # 472_mplowest _63
#     #'xozr8f7p7g': 'iv',  # Mp 2nd sg=38                    
#     '6tmjv4myvg': 'v',  # 1D sg=1
#     #'9lmkmh8s8r': '', # 489_alpha
#     #'zimixdvdxd': '', #492_alpha_like
#     'b5cgvsb16w': '(3)', #'rutile-like',
#     '8p8evt9pcg': '(1)', #'alpha',
#     'zimixdvdxd': '(2)', #'P6_322',
#     'mj7wbfb5nt': '(4)', #'sg=52, battery?',
#     '949rnem5z2': '(5)'   #'sg=53',
#     }

# dx = 0.2

# + Collapsed="false" jupyter={}
# #############################################################################
# Duplicates list
# duplicates = pickle.load(open("../duplicates.pickle", "rb"))


# #############################################################################
# # Bulk DFT Dataframe
# sys.path.insert(0, os.path.join(
#     os.environ["PROJ_irox"], "workflow/ml_modelling"))
# from ml_methods import get_data_for_al

# data_dict = get_data_for_al(stoich="AB2", drop_too_many_atoms=True)
# df_bulk_dft_ab2 = data_dict["df_bulk_dft"]

# data_dict = get_data_for_al(stoich="AB3", drop_too_many_atoms=True)
# df_bulk_dft_ab3 = data_dict["df_bulk_dft"]

# # Combine AB2/3 Dataframes
# df_bulk_dft = pd.concat([df_bulk_dft_ab2, df_bulk_dft_ab3])
# df_bulk_dft = df_bulk_dft[df_bulk_dft.source == "raul"]


# + Collapsed="false" attributes={"classes": [], "id": "", "n": "1"} jupyter={}
# import sys
# import copy
# import pickle

# import numpy as np

# from plotly.subplots import make_subplots

# import chart_studio.plotly as py
# import plotly.graph_objs as go
# import plotly.express as px
# import plotly.io as plio


# from inputs import (
#     Ir_ref,
#     O_ref,
#     coord_env_style)

# from layout import layout
