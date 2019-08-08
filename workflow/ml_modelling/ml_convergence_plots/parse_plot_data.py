#!/usr/bin/env python

"""TEMP.

Author(s): Raul A. Flores
"""

#| - IMPORT MODULES
import pandas as pd
import csv
#__|


def get_plot_data():
    """
    """
    #| - get_plot_data ********************************************************

    #| - IrO2

    #| - IrO2 | round 1
    with open("chris_data_for_figures/iro2/IrO2round1.csv", "r") as fle:
        data = csv.reader(fle)

        data_list = []
        for i in data:
            row_i = {
                "ind": int(i[0]),
                "form_e_0": float(i[1]),
                "uncert_0": float(i[2]),
                "uncert_1": float(i[3]),
                }

            data_list.append(row_i)

    df0 = pd.DataFrame(data_list)
    #__|

    #| - IrO2 | round 2
    with open("chris_data_for_figures/iro2/IrO2round2.csv", "r") as fle:
        data = csv.reader(fle)

        data_list = []
        for i in data:
            row_i = {
                "ind": int(i[0]),
                "form_e_0": float(i[1]),
                "uncert_0": float(i[2]),
                "uncert_1": float(i[3]),
                }

            data_list.append(row_i)

    df1 = pd.DataFrame(data_list)
    #__|

    #| - IrO2 | round 3
    with open("chris_data_for_figures/iro2/IrO2round3.csv", "r") as fle:
        data = csv.reader(fle)

        data_list = []
        for i in data:
            row_i = {
                "ind": int(i[0]),
                "form_e_0": float(i[1]),
                "uncert_0": float(i[2]),
                "uncert_1": float(i[3]),
                }

            data_list.append(row_i)

    df2 = pd.DataFrame(data_list)
    #__|

    d = {
        "round_1": df0,
        "round_2": df1,
        "round_3": df2,
        }

    df_iro2 = pd.concat(d, axis=1)
    #__|

    #| - IrO3

    #| - IrO3 | round 1
    with open("chris_data_for_figures/iro3/IrO3Model1.csv", "r") as fle:
        data = csv.reader(fle)

        data_list = []
        for i, row_i in enumerate(data):

            if i == 0:
                continue

            row_i = {
                "id": row_i[0],
                "form_e_0": float(row_i[1]),
                "uncert_0": float(row_i[2]),
                "uncert_1": float(row_i[2]),
                }

            data_list.append(row_i)

    df0 = pd.DataFrame(data_list)
    #__|

    #| - IrO3 | round 2
    with open("chris_data_for_figures/iro3/IrO3Last.csv", "r") as fle:
        data = csv.reader(fle)

        data_list = []
        for i, row_i in enumerate(data):

            if i == 0:
                continue

            row_i = {
                "id": row_i[0],
                "form_e_0": float(row_i[1]),
                "uncert_0": float(row_i[2]),
                "uncert_1": float(row_i[2]),
                }

            data_list.append(row_i)

    df1 = pd.DataFrame(data_list)

    #__|

    df0["ind"] = [i + 1 for i in range(len(df0))]
    df1["ind"] = [i + 1 for i in range(len(df1))]

    d = {
        "round_1": df0,
        "round_2": df1,
        }


    df_iro3 = pd.concat(d, axis=1)
    #__|form_e_0

    out_dict = {
        "iro2": df_iro2,
        "iro3": df_iro3}

    return(out_dict)

    #__| **********************************************************************
