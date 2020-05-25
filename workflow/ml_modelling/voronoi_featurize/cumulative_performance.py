# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.2
#   kernelspec:
#     display_name: Python [conda env:research-new]
#     language: python
#     name: conda-env-research-new-py
# ---



# # Cumaltive plots

# + jupyter={}
y_tot_sum = AL.CandidateSpace.Y_data_series["y_real"].sum()

Y_data_series = AL.CandidateSpace.Y_data_series

# + jupyter={}
y_min = Y_data_series["y_real"].min()
y_max = Y_data_series["y_real"].max()

m = (1 - 0) / (y_min - y_max)
b = y_max / (y_max - y_min)


Y_data_series["y_real_scaled"] = m * Y_data_series["y_real"] + b

# + jupyter={}
# #############################################################################
# #############################################################################
Y_data_series_sorted = copy.deepcopy(Y_data_series).sort_values("y_real_scaled", ascending=False)

y_gen_summed_ideal = []
for i_cnt, (gen_i, AL_i) in enumerate(AL.al_gen_dict.items()):
    acq_size_i = len(AL_i.prev_acquisition)

    Y_data_acq_i = Y_data_series_sorted.iloc[0:acq_size_i]
    
    y_gen_sum_i = Y_data_acq_i["y_real_scaled"].sum()
    y_gen_summed_ideal.append(y_gen_sum_i)

    # Drop acquired rows
    Y_data_series_sorted = Y_data_series_sorted.drop(
        labels=Y_data_acq_i.index.tolist(), axis=0)

# #############################################################################
# #############################################################################
Y_data_series_sorted = copy.deepcopy(Y_data_series).sort_values("y_real_scaled", ascending=True)

y_gen_summed_nonideal = []
for i_cnt, (gen_i, AL_i) in enumerate(AL.al_gen_dict.items()):
    acq_size_i = len(AL_i.prev_acquisition)

    Y_data_acq_i = Y_data_series_sorted.iloc[0:acq_size_i]
    
    y_gen_sum_i = Y_data_acq_i["y_real_scaled"].sum()
    y_gen_summed_nonideal.append(y_gen_sum_i)

    # Drop acquired rows
    Y_data_series_sorted = Y_data_series_sorted.drop(
        labels=Y_data_acq_i.index.tolist(), axis=0)

# #############################################################################
# #############################################################################
y_gen_summed_al = []
for i_cnt, (gen_i, AL_i) in enumerate(AL.al_gen_dict.items()):
    # #########################################################################
    y_gen_sum_i = 0.
    for id_i in AL_i.prev_acquisition:
        y_i = Y_data_series.loc[id_i]["y_real_scaled"]
        y_gen_sum_i += y_i

    y_gen_summed_al.append(y_gen_sum_i)

# + jupyter={}
# print("y_tot_sum:", -1. * y_tot_sum)

print(np.sum(y_gen_summed_al))
print(np.sum(y_gen_summed_ideal))

# + jupyter={}
# #############################################################################
y_cumulative_al = []
for i_cnt, y_i in enumerate(y_gen_summed_al):
    if i_cnt == 0: y_cumulative_al.append(y_i)
    else: y_cumulative_al.append(y_i + y_cumulative_al[-1])
y_cumulative_al = np.array(y_cumulative_al)


# #############################################################################
y_cumulative_ideal = []
for i_cnt, y_i in enumerate(y_gen_summed_ideal):
    if i_cnt == 0: y_cumulative_ideal.append(y_i)
    else: y_cumulative_ideal.append(y_i + y_cumulative_ideal[-1])
y_cumulative_ideal = np.array(y_cumulative_ideal)


# #############################################################################
y_cumulative_nonideal = []
for i_cnt, y_i in enumerate(y_gen_summed_nonideal):
    if i_cnt == 0: y_cumulative_nonideal.append(y_i)
    else: y_cumulative_nonideal.append(y_i + y_cumulative_nonideal[-1])
y_cumulative_nonideal = np.array(y_cumulative_nonideal)

# + jupyter={}
y_cumulative_al = y_cumulative_al / y_cumulative_al.max()
y_cumulative_ideal = y_cumulative_ideal / y_cumulative_ideal.max()
y_cumulative_nonideal = y_cumulative_nonideal / y_cumulative_nonideal.max()

# + jupyter={}
import chart_studio.plotly as py
import plotly.graph_objs as go
import os


trace = go.Scatter(
    y=y_cumulative_al,
    mode="markers",
    name="AL",
    )

trace_ideal = go.Scatter(
    y=y_cumulative_ideal,
    mode="markers",
    name="ideal",
    )


trace_nonideal = go.Scatter(
    y=y_cumulative_nonideal,
    mode="markers",
    name="non-ideal",
    )

data = [trace, trace_ideal, trace_nonideal]

fig = go.Figure(data=data)
fig.show()

# + jupyter={}
from scipy.integrate import simps
from numpy import trapz

a = len(y_cumulative_al)
b = len(y_cumulative_ideal)
c = len(y_cumulative_nonideal)
assert a == b == c, "IJFIDJS"

area_al = []
area_ideal = []
area_nonideal = []
for i in range(len(y_cumulative_al)):
    area_al_i = trapz(y_cumulative_al[0:i], x=None, dx=1.0)
    area_ideal_i = trapz(y_cumulative_ideal[0:i], x=None, dx=1.0)
    area_nonideal_i = trapz(y_cumulative_nonideal[0:i], x=None, dx=1.0)
    
    area_al.append(area_al_i)
    area_ideal.append(area_ideal_i)
    area_nonideal.append(area_nonideal_i)

# + jupyter={}
data = []

trace_i = go.Scatter(
    y=(np.array(area_al) - np.array(area_nonideal)) / (np.array(area_ideal) - np.array(area_nonideal)),
    mode="markers",
    name="ratio",
    )
data.append(trace_i)


trace_i = go.Scatter(
    y=(np.array(area_al) - np.array(area_nonideal)),
    mode="markers",
    name="AL",
    )
data.append(trace_i)

trace_i = go.Scatter(
    y=(np.array(area_ideal) - np.array(area_nonideal)),
    mode="markers",
    name="ideal",
    )
data.append(trace_i)


fig = go.Figure(data=data)
fig.show()

# + jupyter={}
np.array(area_al) - np.array(area_nonideal)
