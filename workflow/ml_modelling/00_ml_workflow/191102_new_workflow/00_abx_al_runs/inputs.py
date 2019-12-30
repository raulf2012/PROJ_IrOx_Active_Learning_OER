# TEMP

stoich_i = "AB2"
verbose = False

# num_gen_stop = 5
# num_gen_stop = 10

gp_settings = {
    "noise": 0.02542,
    "sigma_l": 1.0049,
    "sigma_f": 5.19,
    "alpha": 0.018,
    }


# name_i = "AL_test_0"
# acquisition_method = "gp_ucb"
# acquisition_method = "random"


# #############################################################################
runs_list = list(range(5))
acquisition_methods = ["gp_ucb", "random"]
duplicate_analysis = [True, False]

# TEST SETTINGS # #############################################################
runs_list = list(range(1))
acquisition_methods = ["gp_ucb"]
duplicate_analysis = [True]