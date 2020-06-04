# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
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

import subprocess

import pandas as pd

# +
import time

t0 = time.time()
# -

orig_dir = os.getcwd()

# Don't change the order of the `notebooks_to_run_list` list! The order is criticial for this to work

# +
# file_i = "workflow/ml_modelling/processing_bulk_dft/static_prototypes_structures/create_atoms_df.ipynb"

notebooks_to_run_list = [
    # #####################################################
    # Figure S1 (MAE vs PCA comp)
                        "workflow/ml_modelling/processing_bulk_dft/static_prototypes_structures/create_atoms_df.ipynb",
                    "workflow/ml_modelling/processing_bulk_dft/parse_oqmd_data/prepare_oqmd_data.ipynb ",
                    "workflow/ml_modelling/processing_bulk_dft/parse_my_oer_bulk_dft/an_parse_oer_bulk_systems.ipynb",
                    "workflow/ml_modelling/processing_bulk_dft/parse_my_bulk_dft/parse_my_bulk_dft.ipynb",
                "workflow/ml_modelling/processing_bulk_dft/collect_all_bulk_dft_data.ipynb",
                        "workflow/ml_modelling/ccf_similarity_analysis/compute_ccf_and_dij_matrix/compute_ccf__v1.ipynb",
                    "workflow/ml_modelling/ccf_similarity_analysis/compute_ccf_and_dij_matrix/get_d_ij_matrix__v1.ipynb",
                    "workflow/ml_modelling/processing_bulk_dft/static_prototypes_structures/elim_high_num_atom_cnt.ipynb",
                "workflow/ml_modelling/00_ml_workflow/get_duplicates_from_al/get_duplicates_list.ipynb",
            "workflow/ml_modelling/processing_bulk_dft/creating_final_dataset_for_upload/create_final_dft.ipynb",
            "workflow/ml_modelling/voronoi_featurize/02_fingerprints_pre_opt.ipynb",
        "workflow/ml_modelling/opt_mae_err_gp_model/cv_error__v3.ipynb",
    "workflow/ml_modelling/opt_mae_err_gp_model/plotting/plot_cv_error.ipynb",

    # #####################################################
    # Figure S2 (IrO2 AL Figure)
        "workflow/ml_modelling/00_ml_workflow/performance_comp/top_10_disc_vs_dft/comp_perf_plots__v2.ipynb",
        "workflow/ml_modelling/00_ml_workflow/al_plots_for_figure/ml_plots__v5.ipynb",
    "workflow/ml_modelling/00_ml_workflow/combined_al_plot/create_subplots__v5.ipynb",

    # #####################################################
    # Figure S3 (All AL performance traces)
        "workflow/ml_modelling/00_ml_workflow/performance_comp/top_10_disc_vs_dft/comp_perf_plots__v2.ipynb",
    "workflow/ml_modelling/00_ml_workflow/performance_comp/top_10_disc_vs_dft/01_plot_all_runs/plot_all_runs.ipynb",

    # #####################################################
    # Figure S4 (Bulk Pourbaix no IrO3)
            "workflow/energy_treatment_deriv/calc_references/calc_O_Ir_refs__H_G.ipynb",
        "workflow/07_bulk_pourbaix/01_pourbaix_scripts/sc_create_all_entries.py",
    "workflow/07_bulk_pourbaix/01_pourbaix_scripts/an_pourbaix_plot.ipynb",

    # #####################################################
    # Figure S5 (OER Scaling)
        "parse_dft_data/parse_all_data_new.ipynb",
    "workflow/02_oer_analysis/03_ads_e_scaling/an_irox_scaling__v2.ipynb",

    # #####################################################
    # Figure 2 (IrO3 AL Figure)
    # Taken care of by previous run of IrO2

    # #####################################################
    # Figure 3 (GP parity plot)
        "workflow/ml_modelling/00_ml_workflow/parity_plots/parity_pre_dft.ipynb",
        "workflow/ml_modelling/00_ml_workflow/parity_plots/parity_post_dft.ipynb",
    "workflow/ml_modelling/00_ml_workflow/parity_plots/plotting_results.ipynb",

    # #####################################################
    # Figure 4 (Bulk Pourbaix w/ IrO3)
    # Taken care of above

    # #####################################################
    # Figure 5 (OER Plot)
        "workflow/07_bulk_pourbaix/01_pourbaix_scripts/sandbox_pourb_transitions.ipynb",
            "workflow/02_oer_analysis/02_oer_volc/kinetic_volcano_colin__v1.ipynb",
        "workflow/02_oer_analysis/02_oer_volc/an_irox_volcano_plotly__v1.ipynb",
            "workflow/energy_treatment_deriv/energy_derivation.ipynb",
        "workflow/01_surface_energies/02_surface_e_pourb_plot/an_surface-energy_pourbaix__v4.ipynb",
    "workflow/03_combined_oer_surf_e/an_surface-energy_pourbaix__v4.ipynb",
    ]
# -

notebooks_that_take_stoich_args = [
    "comp_perf_plots__v2.ipynb",

    # Bulk Pourbaix
    "sc_create_all_entries.py",
    "an_pourbaix_plot.ipynb",

    # AL IrOx Plots
    "create_subplots__v5.ipynb",

    # Parity Plot
    "parity_post_dft.ipynb",
    "parity_pre_dft.ipynb",
    ]


data_dict_list = []
for file_i in notebooks_to_run_list:
    data_dict_i = dict()

    t0_i = time.time()

    file_name = file_i.split("/")[-1]
    file_ext = file_name.split(".")[-1]

    print("8356", 55 * "#")
    print(file_i)
    print(file_name)
    print("8356", 55 * "#")

    data_dict_i["full_path"] = file_i
    data_dict_i["notebook_name"] = file_name

    # #########################################################
    file_path = "/".join(file_i.split("/")[0:-1])


    python_file_name = file_name.split(".")[0] + ".py"
    python_file_full_path = os.path.join(
        os.environ["PROJ_irox"],
        file_path, python_file_name)



    # #########################################################
    # Change directory
    os.chdir(
        os.path.join(
            os.environ["PROJ_irox"],
            file_path))

    # #########################################################
    # Convert .ipynb file to .py ##############################
    bash_comm = "jupytext --to py " + file_name

    output = None
    if not file_ext == "py":
        try:
                output = subprocess.check_output(
                    bash_comm,
                    stderr=subprocess.STDOUT,
                    shell=True)
        except subprocess.CalledProcessError as e:
            output = e.output

        output = output.decode("utf-8")
        print(output)

    # #########################################################
    # Run .py file ############################################

    if file_name in notebooks_that_take_stoich_args:
        bash_comm = "python " + python_file_full_path + " AB2"

        output = None
        try:
            output = subprocess.check_output(
                bash_comm, stderr=subprocess.STDOUT, shell=True)
        except subprocess.CalledProcessError as e:
            output = e.output
        output = output.decode("utf-8")
        print(output)

        # I'll run the next one outside of the if statement as normal by replacing the bash_comm
        bash_comm = "python " + python_file_full_path + " AB3"

    else:
        bash_comm = "python " + python_file_full_path

    output = None
    try:
        output = subprocess.check_output(
            bash_comm,
            stderr=subprocess.STDOUT,
            shell=True)
    except subprocess.CalledProcessError as e:
        output = e.output

    output = output.decode("utf-8")
    print(output)

    # #########################################################
    completed_correctly = "All done!" in output
    data_dict_i["completed"] = completed_correctly

    # #########################################################
    # Change directory
    os.chdir(orig_dir)


    tf_i = time.time()
    notebook_run_time = tf_i - t0_i
    data_dict_i["run_time_s"] = notebook_run_time
    data_dict_i["run_time_min"] = notebook_run_time / 60

    data_dict_list.append(data_dict_i)
    print(3 * "\n")

# +
df = pd.DataFrame(data_dict_list)

# df.iloc[4].full_path
# df.sort_values("run_time_min")

df
# -

os.getcwd()

# +
run_time = time.time() - t0
run_time = run_time / 60

print(run_time, "min")

# + active=""
#
#
#
#
#

# +

# # #########################################################
# #
# file_name = file_i.split("/")[-1]

# file_path = "/".join(file_i.split("/")[0:-1])

# python_file_name = file_name.split(".")[0] + ".py"

# python_file_full_path = os.path.join(
#     os.environ["PROJ_irox"],
#     file_path, python_file_name)



# # #########################################################
# # Change directory
# os.chdir(
#     os.path.join(
#         os.environ["PROJ_irox"],
#         file_path))

# # #########################################################
# # Convert .ipynb file to .py ##############################
# bash_comm = "jupytext --to py " + file_name

# output = None
# try:
#     output = subprocess.check_output(
#         bash_comm,
#         stderr=subprocess.STDOUT,
#         shell=True)
# except subprocess.CalledProcessError as e:
#     output = e.output

# output = output.decode("utf-8")
# print(output)

# # #########################################################
# # Run .py file ############################################
# bash_comm = "python " + python_file_full_path

# output = None
# try:
#     output = subprocess.check_output(
#         bash_comm,
#         stderr=subprocess.STDOUT,
#         shell=True)
# except subprocess.CalledProcessError as e:
#     output = e.output

# output = output.decode("utf-8")
# print(output)

# # #########################################################
# #


# # #########################################################
# # Change directory
# os.chdir(orig_dir)

# +
# os.getcwd()

# os.system("ls")

# +
# bash_comm = "run_jupy " + file_name

# os.system(bash_comm)

# +
# bash_comm = "python " + file_name

# # os.system(bash_comm)

# +
# import subprocess

# correct = subprocess.run(
#     ["python", file_name],
#     # check=True,
#     # text=True,
#     shell=True
#     )

# +
# file_path

# # file_i

# # file_path =
# # "/".join(file_i.split("/")[0:-1])
# file_i


# df.iloc[15].full_path


# +
# import subprocess

# bash_comm = "python " + python_file_name
# # bash_comm = "python " + python_file_full_path
# print(bash_comm)

# output = subprocess.check_output(
#     bash_comm,
#     stderr=subprocess.STDOUT,
#     shell=True,
#     )
# output

# +
# type(str(output))
# print(str(output))

# str(output)[0:10]

# bash_comm
# correct
