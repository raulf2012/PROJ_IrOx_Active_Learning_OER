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

# +
import os
print(os.getcwd())
import sys

import numpy as np

# +
# %%capture

sys.path.insert(
    0, os.path.join(
        os.environ["PROJ_irox"],
        "workflow"))

sys.path.insert(
    0, os.path.join(
        os.environ["PROJ_irox"],
        "data"))

from an_data_processing import load_df

from proj_data_irox import data_dir
# -

df_pourbaix, df_ads, df_surf = load_df(
    from_file=True,
    root_dir=data_dir,
    data_dir=data_dir,
    file_name="df_master.pickle",
    process_df=False)
df_m = df_ads

# + active=""
#
#
#
# -

adsorbates = ["o", "oh", "ooh"]
# adsorbates = ["adsorbate"]

# +
df_i = df_m[(df_m.bulk_system == "IrO2")]

mean_ads_dict = dict()
group = df_i.groupby("adsorbate")
for ads, df_j in group:
    ads_e_ave = df_j.ads_e.mean()
    mean_ads_dict[ads] = ads_e_ave

mean_ads_dict_ab2 = mean_ads_dict

# +
df_i = df_m[(df_m.bulk_system == "IrO3")]

mean_ads_dict = dict()
group = df_i.groupby("adsorbate")
for ads, df_j in group:
    ads_e_ave = df_j.ads_e.mean()
    mean_ads_dict[ads] = ads_e_ave
    
mean_ads_dict_ab3 = mean_ads_dict
# -

ads_diff_dict = dict()
for ads_i in adsorbates:
    ads_ab2 = mean_ads_dict_ab2[ads_i]
    ads_ab3 = mean_ads_dict_ab3[ads_i]
    
    ads_diff = ads_ab3 - ads_ab2
    
    print(ads_i, ":", ads_diff)
    ads_diff_dict[ads_i] = ads_diff

# +
# for i, j in ads_diff_dict.items():
#     print(j)

np.mean(list(ads_diff_dict.values()))
# -

ads_diff_dict
