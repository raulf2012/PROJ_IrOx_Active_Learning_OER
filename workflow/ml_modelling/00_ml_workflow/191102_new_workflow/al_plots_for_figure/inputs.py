stoich_i = "AB3"
# stoich_i = "AB2"

lowest_N_sys_to_track = 10

# | - Generations to plot
# generations to plot

# AB3
if stoich_i == "AB3":
#     gen_0 = 0
#     gen_1 = 2
#     gen_2 = 9  # 5/10 top structures found
#     gen_3 = 14 # 10/10 top structures found
#     gen_4 = "last"

#     gen_0 = 0
#     gen_1 = 3
#     gen_2 = 6  # 7/10 top structures found
#     gen_3 = 20 # 10/10 top structures found
#     gen_4 = -10

    gen_0 = 0   # 00/10
    gen_1 = 3   # 00/10
    gen_2 = 5   # 04/10
    gen_3 = 14  # 09/10
    gen_4 = 20 # 10/10


# AB2  # 92 total gens
# 01 found | gen 3
# 05 found | gen 5
# 09 found | gen 30
# 10 found | gen 75

if stoich_i == "AB2":
    gen_0 = 0
    gen_1 = 3
    # Gen 5 looks really bad, maybe use a slightly better trained version
    gen_2 = 12  # 5/10 top structures found
    # gen_2 = 5  # 5/10 top structures found
    gen_3 = 30 # 9/10 top structures found
    gen_4 = "last"


gens_to_plot = [gen_0, gen_1, gen_2, gen_3, gen_4]

main_gen = gen_2
#__|


#| - Minor tick axis settings
minor_ticks_data = {

#| - Main AL Plots
    "main_al_plots": {
        "x": 100,
        "y": 0.25,
        },

# # #############################################################################
# add_duplicate_axes(
#     fig,
#     axis_type='x',
#     axis_data={"dtick": 100, "tickcolor": "black", "ticklen": 3},
#     axis_num_list=[7, 8, 9, 10, 11],
#     # tmp_define_both_axis_types=False,
#     tmp_define_both_axis_types=True,
#     )
# add_duplicate_axes(
#     fig,
#     axis_type='y',
#     axis_data={"dtick": 0.25, "tickcolor": "black", "ticklen": 3},
#     axis_num_list=[7, 8, 9, 10, 11],
#     # tmp_define_both_axis_types=False,
#     tmp_define_both_axis_types=True,
#     )

#__|

    }

# # #############################################################################
# add_duplicate_axes(
#     fig,
#     axis_type='x',
#     axis_data={"dtick": 1, "tickcolor": "black", "ticklen": 3},
#     # axis_num_list=[12, 13, 14, 15, 16],
#     axis_num_list=[12,],
#     # tmp_define_both_axis_types=False,
#     tmp_define_both_axis_types=True,
#     )
# add_duplicate_axes(
#     fig,
#     axis_type='y',
#     axis_data={"dtick": 0.1, "tickcolor": "black", "ticklen": 3},
#     # axis_num_list=[12, 13, 14, 15, 16],
#     axis_num_list=[12, ],
#     # tmp_define_both_axis_types=False,
#     tmp_define_both_axis_types=True,
#     )

# add_duplicate_axes(
#     fig,
#     axis_type='x',
#     axis_data={"dtick": 1, "tickcolor": "black", "ticklen": 3},
#     # axis_num_list=[12, 13, 14, 15, 16],
#     axis_num_list=[13,],
#     # tmp_define_both_axis_types=False,
#     tmp_define_both_axis_types=True,
#     )
# add_duplicate_axes(
#     fig,
#     axis_type='y',
#     axis_data={"dtick": 0.05, "tickcolor": "black", "ticklen": 3},
#     # axis_num_list=[12, 13, 14, 15, 16],
#     axis_num_list=[13, ],
#     # tmp_define_both_axis_types=False,
#     tmp_define_both_axis_types=True,
#     )


# add_duplicate_axes(
#     fig,
#     axis_type='x',
#     axis_data={"dtick": 1, "tickcolor": "black", "ticklen": 3},
#     # axis_num_list=[12, 13, 14, 15, 16],
#     axis_num_list=[14,],
#     # tmp_define_both_axis_types=False,
#     tmp_define_both_axis_types=True,
#     )
# add_duplicate_axes(
#     fig,
#     axis_type='y',
#     axis_data={"dtick": 0.05, "tickcolor": "black", "ticklen": 3},
#     # axis_num_list=[12, 13, 14, 15, 16],
#     axis_num_list=[14, ],
#     # tmp_define_both_axis_types=False,
#     tmp_define_both_axis_types=True,
#     )

# add_duplicate_axes(
#     fig,
#     axis_type='x',
#     axis_data={"dtick": 1., "tickcolor": "black", "ticklen": 3},
#     # axis_num_list=[12, 13, 14, 15, 16],
#     axis_num_list=[15, 16],
#     # tmp_define_both_axis_types=False,
#     tmp_define_both_axis_types=True,
#     )
# add_duplicate_axes(
#     fig,
#     axis_type='y',
#     axis_data={"dtick": 0.025, "tickcolor": "black", "ticklen": 3},
#     # axis_num_list=[12, 13, 14, 15, 16],
#     axis_num_list=[15, 16],
#     # tmp_define_both_axis_types=False,
#     tmp_define_both_axis_types=True,
#     )

# # #############################################################################

# add_duplicate_axes(
#     fig,
#     axis_type='x',
#     axis_data={"dtick": 50., "tickcolor": "black", "ticklen": 3},
#     axis_num_list=[20, ],
#     tmp_define_both_axis_types=True,
#     )

# add_duplicate_axes(
#     fig,
#     axis_type='y',
#     axis_data={"dtick": 1., "tickcolor": "black", "ticklen": 3},
#     axis_num_list=[20, ],
#     tmp_define_both_axis_types=True,
#     )

#__|

# #############################################################################
# Variablers to export
stoich_i,
lowest_N_sys_to_track,
gens_to_plot,
main_gen,
