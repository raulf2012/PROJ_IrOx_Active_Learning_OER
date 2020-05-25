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

import ase.db
db = ase.db.connect("FloresActive2020/FloresActive2020.db")

for row in db.select():
    tmp = 42
    # print(row.forces[0, 2], row.relaxed)

row.key_value_pairs
