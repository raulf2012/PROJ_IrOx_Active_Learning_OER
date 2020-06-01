# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_irox]
#     language: python
#     name: conda-env-PROJ_irox-py
# ---

# +
import os
print(os.getcwd())
import sys

from json import dump, load
from shutil import copyfile

# #########################################################
from methods import clean_ipynb
from methods import get_ipynb_notebook_paths
# -

dirs_list = get_ipynb_notebook_paths()
for file_i in dirs_list:
    clean_ipynb(file_i, True)
