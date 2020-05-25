#!/usr/bin/env python

"""Job maintance methods.

Author: Raul A. Flores
"""

# | - Import Modules
import os
import shutil
from ase.build import add_adsorbate

# My Modules
from ase_modules.dft_params import Espresso_Params
from ase_modules.ase_methods import find_diff_between_atoms_objects
#__|

def dir_setup(step_i, path_i, job_i_params, wf_vars):
    """Everything needed to specify a job completely.

    Args:
        step_i:
        path_i:
        job_i_params:
        wf_vars:
    """
    # | - dir_setup ************************************************************
    from ase import io
    # from ase_modules.adsorbates import Adsorbate

    # If the .READY file exists in a dir, don't do anything, the job has
    # already been setup
    if os.path.exists(os.path.join(path_i, ".READY")):
        # print("Job already setup, skipping setup:")
        # print(path_i)
        return(None)

    # | - Parsing wf_vars
    master_root_dir = wf_vars["root_dir"]
    mod_dir = wf_vars["mod_dir"]
    model_names = wf_vars["model_names"]
    # JobsAn = wf_vars["jobs_an_list"][step_i]
    atoms_list_names = wf_vars["atoms_list_names"]
    atoms_ext = wf_vars["atoms_ext"]
    atoms_dir = wf_vars["atoms_dir"]
    #__|

    # | - Copying Model Scripts
    model_dir = master_root_dir + "/" + mod_dir + "/" + model_names[step_i]
    shutil.copyfile(model_dir, path_i + "/model.py")
    #__|

    # | - Atoms Object
    for atoms_info_i in atoms_list_names:
        params_i = atoms_info_i["job_params"]

        if job_i_params == params_i:
            atoms_path_i = atoms_info_i["atoms_path"]
            break

    slab = io.read(os.path.join(
        # master_root_dir,
        ".",
        atoms_dir,
        atoms_path_i,
        ))

    slab.write(path_i + "/init.traj")

    # ads_i = job_i_params["adsorbate"]
    # site_i = job_i_params["site"]
    #
    # z_dist = 1.3
    # site_coords_dict = {
    #     "ring-center":  [3.201, 1.915, z_dist],
    #     "C1-ontop":     [2.493, 0.613, z_dist],
    #     "C2-ontop":     [1.774, 1.843, z_dist],
    #     "N-ontop":      [3.901, 0.615, z_dist],
    #     "C-N-bridged":  [3.201, 0.615, z_dist],
    #     "C-C bridged":  [2.201, 1.215, z_dist],
    #     }
    #
    # Ads = Adsorbate()
    # ads = Ads.get_adsorbate(ads_i)
    #
    # atoms_file = atoms_list_names[step_i] + atoms_ext
    # slab = io.read(master_root_dir + "/" + atoms_dir + "/" + atoms_file)
    #
    # z_space = site_coords_dict[site_i][2]
    #
    # slab_before_adsorbate = slab.copy()
    #
    # add_adsorbate(slab, ads, z_space, position=site_coords_dict[site_i][0:2])
    #
    # tmp1, tmp2 = find_diff_between_atoms_objects(slab_before_adsorbate, slab)
    #
    # # print("Adsorbate indices - klsjfksjkfjdkk")
    # # print(tmp2)
    #
    # slab.info["adsorbates"] = tmp2
    #
    # slab.write(path_i + "/init.traj")
    #__|

    # | - Writing Ready File for First Step
    if step_i == 0:
        open(path_i + "/.READY", "w")
    #__|

    #__| **********************************************************************

def job_maint(
    step_i,
    job_i,
    job_i_params,
    wf_vars,
    tally,
    file_ops=False,
    ):
    """Job maintance, rerunning jobs, continuning jobs, etc.

    Args:
        step_i:
        job_i:
        job_i_params:
        wf_vars:
        tally:
        file_ops:
            If False will not create new revision folders and copy files into
            them.
    """
    # | - job_maint ************************************************************

    # | - Parsing wf_vars
    # master_root_dir = wf_vars["root_dir"]
    # mod_dir = wf_vars["mod_dir"]
    # model_names = wf_vars["model_names"]

    Jobs = wf_vars["jobs_man_list"][step_i]

    # atoms_list_names = wf_vars["atoms_list_names"]
    # atoms_ext = wf_vars["atoms_ext"]
    # atoms_dir = wf_vars["atoms_dir"]
    #__|

    import copy
    # path_i = Jobs.var_lst_to_path(job_i, job_rev="Auto", relative_path=False)
    path_i = copy.deepcopy(job_i)

    # | - Not READY | Job Dir Not Set Up
    # if Jobs._job_setup(job_i):
        # print("SUBMITTING" + " | " + path_i)
    #__|

    # | - READY  | Not SUBMITTED
    if Jobs._job_ready(job_i):
        print("SUBMITTING" + " | " + path_i)

        # Jobs.submit_job(
        #     path_i=path_i,
        #     wall_time="400",
        #     queue="regular",
        #     )
    #__|

    # | - PENDING | Leave Alone
    # print("Is job pending:")
    # print(Jobs._job_pending(job_i))

    if Jobs._job_pending(job_i):
        print("PENDING" + " | " + path_i)  # PERM_PRINT

        tally["pending"] += 1
    #__|

    # | - RUNNING
    if Jobs._job_running(job_i):
        print("RUNNING" + " | " + path_i)  # PERM_PRINT

        tally["running"] += 1
    #__|

    # | - SUCCEEDED | Copy Files to Next Step
    if Jobs._job_succeeded(job_i):
        print("SUCCESS" + " | " + path_i)  # PERM_PRINT

        tally["successes"] += 1

        if file_ops:
            tmp = 42
            print(tmp)

            # if step_num == 1:
            #     jd.copyfiles_onestep_up(job_i, step_num, Jobs_Inst_list,
            #         files_lst=[
            #             ["CONTCAR", "init.POSCAR"],
            #             ".READY",
            #             ]
            #         )

    # else:
    #     pass
    #__|

    # | - FAILED | Rerun Job
    if Jobs._job_failed(job_i):
        print("FAILURE" + " | " + path_i)  # PERM_PRINT
        tally["failures"] += 1

        # | - Step 1 Restart
        if step_i == 0:
            prev_files = [
                "dft-params.json",
                "out_opt.traj",
                "model.py",
                ".READY",
                ]

            sub_params = {
                "wall_time": "2000"
                }

            if file_ops:
                Jobs.restart_job(job_i, prev_files, sub_params=sub_params)
        #__|

    #__|


    return(tally)
    #__| **********************************************************************
