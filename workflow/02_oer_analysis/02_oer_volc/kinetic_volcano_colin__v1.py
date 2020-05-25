# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python [conda env:PROJ_irox] *
#     language: python
#     name: conda-env-PROJ_irox-py
# ---

# # Kinetic OER Volcano (From Colin's microkinetic OER paper)
# ## Insights into the Electrochemical Oxygen Evolution Reaction with ab Initio Calculations and Microkinetic Modeling: Beyond the Limiting Potential Volcano
# https://pubs.acs.org/doi/10.1021/acs.jpcc.9b03830

# # Import Modules

# +
import os
print(os.getcwd())

import sys

import numpy as np
import pandas as pd

from scipy.optimize import minimize


from decimal import Decimal

import chart_studio.plotly as py
import plotly.graph_objs as go
# -

# # Data From Colin's Plot Extracted

# +
df_A = pd.read_csv(
    os.path.join(
        # os.environ["PROJ_DATA"],
        # "04_IrOx_surfaces_OER/kinetic_oer_volcano_data",
        "in_data",
        "kinetic_volc_1A.csv"))

df_mA = pd.read_csv(
    os.path.join(
        # os.environ["PROJ_DATA"],
        # "04_IrOx_surfaces_OER/kinetic_oer_volcano_data",
        "in_data",
        "kinetic_volc_1mA.csv"))

data = []

df = df_A
trace = go.Scatter(
    x=df["x"],
    y=df["y"],
    mode="lines",
    line=dict(
        color="firebrick",
        width=2,
        dash="dot",
        ),
    )
data.append(trace)

df = df_mA
trace = go.Scatter(
    x=df["x"],
    y=df["y"],
    mode="lines",
    line=dict(
        color="firebrick",
        width=2,
        dash="dot",
        ),
    )
data.append(trace)

layout = go.Layout(yaxis=dict(range=[2.4, 1.5]))

fig = go.Figure(data=data, layout=layout)
# fig.show()
# -

# ## Saving traces

# Pickling data ######################################################
import os; import pickle
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)
with open(os.path.join(directory, "kinetic_volcano_trace.pickle"), "wb") as fle:
    pickle.dump(data, fle)
# #####################################################################

# # Model Imputs

# +
# #############################################################################
# Constants fit to data OER data in paper
xi = 2.41
alpha = 0.48
beta = 0.58

# Extent of reaction, (close to 0)
# gamma = 0.01
gamma = 0.

# #############################################################################
# Physical constants
T = 300
h = 6.62607015E-34

kb = 1.38064852E-23  # m2 kg s-2 K-1
kb_1 = 8.617333262145E-5

# #############################################################################
# Initial guess for numerical solver
G__O_OH = 0.5

# +
data_kin_volc = data

def rate_eqn(
    xi=None,
    alpha=None,
    beta=None,
    gamma=None,
    kb=None,
    T=None,
    h=None,
    G__O_OH=None,
    U_SHE=None,
    verbose=False,
    ):
    

    # #########################################################################
    # K__OH_O Calculation #####################################################
    K__OH_O = np.exp(-(-G__O_OH + 1 * U_SHE) / (kb_1 * T))
    # K__OH_O = np.exp(-(G__O_OH - 1 * U_SHE) / (kb_1 * T))
    # K__OH_O = np.exp(-(G__O_OH) / (kb_1 * T))

    if verbose:
        print("K__OH_O:", K__OH_O)

    # #########################################################################
    # Constant Term ###########################################################
    term_const_0 = (1 + K__OH_O) ** (-1.)
    term_const_1 = (kb * T / h)
    term_const_2 = (1 - gamma)
    term_const = term_const_0 * term_const_1 * term_const_2
    
    if verbose:
        term_const_tmp = '%.4E' % Decimal(float(term_const))
        print("term_const:", term_const_tmp)

    # #########################################################################
    # Exponential Term ########################################################
    term_exp_0 = -((xi) - alpha * (G__O_OH) - beta * U_SHE) / (kb_1 * T)
    term_exp_1 = np.exp(term_exp_0)

    if verbose:
        print("term_exp_0:", term_exp_0)

        term_exp_tmp = '%.4E' % Decimal(float(term_exp_1))
        print("term_exp:", term_exp_tmp)

    # #########################################################################
    # Putting terms together to calc rate #####################################
    r_1 = term_const * term_exp_1
    r_out = 0.64 * r_1

    return(r_out)


# +

# U_SHE = 2.4483

# rate_eqn(
#     xi=xi,
#     alpha=alpha,
#     beta=beta,
#     gamma=gamma,
#     kb=kb,
#     T=T,
#     h=h,
#     G__O_OH=G__O_OH,
#     U_SHE=U_SHE,
#     verbose=True,
#     )
# -

def test_meth(
    U,
    xi=xi,
    alpha=alpha,
    beta=beta,
    gamma=gamma,
    kb=kb,
    T=T,
    h=h,
    G__O_OH=G__O_OH,
    ):

    rate_i = rate_eqn(
        xi=xi,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        kb=kb,
        T=T,
        h=h,
        G__O_OH=G__O_OH,
        U_SHE=U,
        )
    
    desired_current_density = 10  # mA/cm2
    out = np.abs(rate_i - desired_current_density)

    return(out)


# # Load data from previous run

# +
# #############################################################################
import pickle; import os
path_i = os.path.join(
    "out_data",
    "df_10mA.pickle")

from pathlib import Path
my_file = Path(path_i)
if my_file.is_file():
    with open(path_i, "rb") as fle:
        df_old = pickle.load(fle)
else:
    df_old = pd.DataFrame()
# #############################################################################
# -

# # Solve equation numerically across a range of dG_O-OH values
#
# Change the np.linspace arguments to sample different regions

# +
data_list = []
# for i in np.linspace(1.5, 1.7, 100, endpoint=True):
# for i in np.linspace(0.5, 2.5, 100, endpoint=True):

# for i in np.linspace(1.5, 1.7, 193, endpoint=True):
for i in np.linspace(0.5, 2.5, 100, endpoint=True):
    data_dict_i = dict()
    
    data_dict_i["descriptor"] = i
    print("DG_O - DG_OH", i)
   
    for U_i in np.linspace(-0.5, 3.0, 50, endpoint=True):
        fun = test_meth
        min_out = minimize(
            fun,
            U_i,
            args=(xi, alpha, beta, gamma, kb, T, h, i),
            method="L-BFGS-B",
            tol=0.1,
            options={
                'disp': None, 'maxcor': 10, 'ftol': 2.220446049250313e-09,
                'gtol': 1e-05, 'eps': 1e-09, 'maxfun': 15000, 'maxiter': 15000,
                'iprint': -1, 'maxls': 20})

        if min_out.fun[0] < 1e-8:
            U_lim = min_out.x[0]
            data_dict_i["U_lim"] = U_lim
            data_list.append(data_dict_i)
            print("Found solution")
            break


    print(1 * "\n")
# -

# # Analytical expressions for the left and right legs of volcano
#
# This works but doesn't resolve the "hump" of the volcano, i.e. the top will be pointy

# +
data_kin_volc = data


def leg_0(
    xi=None,
    alpha=None,
    beta=None,
    gamma=None,
    kb=None,
    T=None,
    h=None,
    G__O_OH=None,
    U_SHE=None,
    ):
    verbose = True

    # 1 mA/cm2
    r = 1 / 0.64

    U = ((kb_1 * T * np.log(r * h / (kb * T))) + xi - (G__O_OH) * alpha) / beta

    return(U)


def leg_1(
    xi=None,
    alpha=None,
    beta=None,
    gamma=None,
    kb=None,
    T=None,
    h=None,
    G__O_OH=None,
    U_SHE=None,
    ):
    verbose = True

    # 1 mA/cm2
    r = 1 / 0.64

    U = ((kb_1 * T * np.log(r * h / (kb * T))) + xi + (G__O_OH) * (1 - alpha)) / (1 + beta)

    return(U)

# G__O_OH = 2

leg_args = dict(
    xi=xi,
    alpha=alpha,
    beta=beta,
    gamma=gamma,
    kb=kb,
    T=T,
    h=h)

x_array = np.linspace(0.5, 2.5, 10, endpoint=True)
leg_0_array = [leg_0(G__O_OH=i, **leg_args) for i in x_array]
leg_1_array = [leg_1(G__O_OH=i, **leg_args) for i in x_array]

trace_leg_0 = go.Scatter(
    x=x_array,
    y=leg_0_array,
    mode="lines",
    line=dict(
        color="black",
        width=1,
        ),
    )

trace_leg_1 = go.Scatter(
    x=x_array,
    y=leg_1_array,
    mode="lines",
    line=dict(
        color="black",
        width=1,
        ),
    )

# +
# Combining new data with old
df = pd.DataFrame(data_list)

df = pd.concat([df_old, df])

# Sort by descriptor value
df = df.sort_values("descriptor")
# -

# # Plotting data

# +
x_array = df["descriptor"]
y_array = df["U_lim"]

trace = go.Scatter(
    x=x_array,
    y=y_array,
    mode="markers+lines",
    opacity=0.8,
    marker=dict(
        symbol="circle",
        color='LightSkyBlue',

        opacity=0.8,
        size=8,
        line=dict(
            color='MediumPurple',
            width=2
            )
        ),

    line=dict(
        color="black",
        width=1,
        ),
    )


data = []
data.append(trace)
data.append(trace_leg_0)
data.append(trace_leg_1)

for trace_i in data_kin_volc:
    data.append(trace_i)

# #############################################################################
layout = go.Layout(yaxis=dict(autorange="reversed"))

fig = go.Figure(data=data, layout=layout)
fig.show()
# -

# Pickling data ######################################################
import os; import pickle
directory = "out_data"
if not os.path.exists(directory): os.makedirs(directory)
with open(os.path.join(directory, "df_10mA.pickle"), "wb") as fle:
    pickle.dump(df, fle)
# #####################################################################
