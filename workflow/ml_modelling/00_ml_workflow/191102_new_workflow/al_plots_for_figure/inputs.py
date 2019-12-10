stoich_i = "AB3"

lowest_N_sys_to_track = 10

#| - Generations to plot
# generations to plot

# AB3
if stoich_i == "AB3":
#     gen_0 = 0
#     gen_1 = 2
#     gen_2 = 9  # 5/10 top structures found
#     gen_3 = 14 # 10/10 top structures found
#     gen_4 = "last"

    gen_0 = 0
    gen_1 = 3
    gen_2 = 6  # 7/10 top structures found
    gen_3 = 20 # 10/10 top structures found
    gen_4 = -10


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
# __|



# #############################################################################
# Variablers to export
stoich_i,
lowest_N_sys_to_track,
gens_to_plot,
main_gen,
