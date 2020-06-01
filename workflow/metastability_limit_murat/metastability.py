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

sys.path.insert(0, os.path.join(os.environ["PROJ_irox"], "data"))
from proj_data_irox import calc_dH
# -

raw_dft_most_stable_amorph = dict(
    AB2=-6.542,
    AB3=-6.163,
    )

# +
ab2_meta_lim = calc_dH(
    raw_dft_most_stable_amorph["AB2"],
    stoich="AB2")

ab3_meta_lim = calc_dH(
    raw_dft_most_stable_amorph["AB3"],
    stoich="AB3")

# calc_dH(stoich="AB2")
# -

print("AB2:", ab2_meta_lim)
print("AB3:", ab3_meta_lim)
